"""
Chemical Reaction Sandbox — Streamlit app (thin presentation layer).

Two tabs:
  🧪 Molecular Sandbox     — resolve molecules from PubChem, pack them into a box
                             at a chosen density/solvent/temperature, give them
                             Maxwell–Boltzmann velocities and propagate ab-initio
                             Langevin MD. Live Plotly charts of PE/KE/temperature.
  🏭 Reactor Engineering   — ideal Batch / CSTR / PFR models for an nth-order
                             reaction A → products. Conversion, sizing, Levenspiel
                             and Arrhenius plots.

Compute lives in the ``chemmodel`` package; this module only wires widgets,
threads the MD worker, and renders state. Run with:
  streamlit run src/chemmodel/ui/sandbox_app.py
"""

from __future__ import annotations

import queue
import threading
import time

import numpy as np
import streamlit as st

from chemmodel.chemistry import assemble_system, auto_multiplicity, xyz_to_atoms
from chemmodel.chemistry.pubchem import fetch_pubchem_xyz as _fetch_pubchem_xyz
from chemmodel.constants import EV_TO_KJMOL, L_TO_FT3
from chemmodel.engines import available_engines, make_psi4_calculator, probe
from chemmodel.md import md_worker
from chemmodel.reactors import (
    arrhenius_k,
    batch_profile,
    batch_time_for_X,
    inv_rate,
    k_units,
    tau_cstr,
    tau_pfr,
)
from chemmodel.ui.theme import LANE_COLORS, LANES, SANDBOX_CSS
from chemmodel.ui.viz import md_charts_fig, system_3d_fig

PSI4_OK = available_engines()["psi4"]
PCP_OK = probe("pubchempy")

# process-global live MD worker handle (one Psi4 worker per process)
_MD_WORKER: "threading.Thread | None" = None


@st.cache_data(show_spinner=False)
def fetch_pubchem_xyz(name: str) -> dict:
    """Cached PubChem lookup so reruns/repeat names don't re-hit the network."""
    return _fetch_pubchem_xyz(name)


# ══════════════════════════════════════════════════════════════════════════════
# Page
# ══════════════════════════════════════════════════════════════════════════════

st.set_page_config(page_title="Reaction Sandbox", page_icon="🧪", layout="wide")
st.markdown(SANDBOX_CSS, unsafe_allow_html=True)

st.title("🧪 Chemical Reaction Sandbox")
st.caption("PubChem · ASE · Psi4 (PCM) · ideal-reactor kinetics")

tab_md, tab_reactor = st.tabs(["🧪 Molecular Sandbox", "🏭 Reactor Engineering"])

# ══════════════════════════════════════════════════════════════════════════════
# TAB 1 — Molecular Sandbox (ab-initio Langevin MD)
# ══════════════════════════════════════════════════════════════════════════════

with tab_md:
    # ── availability / environment banners ─────────────────────────────────────
    cols = st.columns(2)
    cols[0].markdown(
        f'<span class="engine-badge">{"✅ Psi4 ready" if PSI4_OK else "⚠️ Psi4 missing — demo mode"}</span>',
        unsafe_allow_html=True)
    cols[1].markdown(
        f'<span class="engine-badge">{"✅ PubChemPy ready" if PCP_OK else "⚠️ PubChemPy missing"}</span>',
        unsafe_allow_html=True)

    if not PSI4_OK:
        st.warning("Psi4 is not importable in this interpreter. Run the app with the "
                   "`chemmodel` conda environment (see README) for real simulations; "
                   "otherwise the MD loop falls back to a cheap analytic force field.")

    # ── input lanes ────────────────────────────────────────────────────────────
    st.subheader("1 · Molecule input lanes")
    st.caption("Type a chemical name (e.g. *ethylene*, *acetic acid*, *palladium*) and "
               "press **Resolve** to fetch 3-D coordinates from PubChem. "
               "At least *Reactant A* is required; other lanes are optional.")

    lane_cols = st.columns(4)
    resolved: dict[str, dict] = {}

    for lane, col in zip(LANES, lane_cols):
        with col:
            st.markdown(
                f'<span class="lane-badge" style="background:{LANE_COLORS[lane]};">{lane}</span>',
                unsafe_allow_html=True)
            name = st.text_input("PubChem name", key=f"name_{lane}",
                                 placeholder="e.g. water", label_visibility="collapsed")
            if st.button("🔍 Resolve", key=f"btn_{lane}", use_container_width=True):
                if name.strip():
                    with st.spinner(f"Querying PubChem for '{name}'…"):
                        st.session_state[f"res_{lane}"] = fetch_pubchem_xyz(name)
                else:
                    st.session_state.pop(f"res_{lane}", None)

            with st.expander("…or paste XYZ"):
                manual = st.text_area("XYZ", key=f"manual_{lane}", height=90,
                                      label_visibility="collapsed")
                if manual.strip():
                    st.session_state[f"res_{lane}"] = {
                        "ok": True, "xyz": manual, "formula": "(manual)",
                        "cid": None, "flat": False, "error": None}

            res = st.session_state.get(f"res_{lane}")
            if res:
                if res["ok"]:
                    resolved[lane] = res
                    tag = f"✅ {res['formula']}"
                    if res.get("cid"):
                        tag += f" · CID {res['cid']}"
                    st.success(tag)
                    if res.get("flat"):
                        st.warning("Only 2-D coords available — geometry is flat.")
                else:
                    st.error(res["error"])

    active = [lane for lane in LANES if lane in resolved]

    # ── build ASE molecules from resolved lanes ────────────────────────────────
    mol_atoms, lane_used, parse_err = [], [], None
    for lane in active:
        try:
            mol_atoms.append(xyz_to_atoms(resolved[lane]["xyz"]))
            lane_used.append(lane)
        except Exception as e:
            parse_err = f"{lane}: {e}"

    if parse_err:
        st.error(f"XYZ parse error — {parse_err}")

    st.divider()

    # ── simulation controls ─────────────────────────────────────────────────────
    st.subheader("2 · Simulation conditions")
    st.caption(
        "Think of this as a tiny **batch reactor**: the box is the vessel (its volume "
        "fixes reactant concentration), the solvent is the charge medium, the "
        "temperature is the operating setpoint, and steps × timestep is how long the "
        "reactor is run. Hover the ⓘ on each control for details."
    )
    c1, c2, c3, c4 = st.columns(4)
    with c1:
        box_L = st.slider(
            "Reactor volume — box edge (Å)", 8.0, 40.0, 16.0, 0.5,
            help="The cubic vessel that holds the charge. Its edge length sets the "
                 "reactor volume (shown below), which fixes reactant **concentration / "
                 "packing density**: a smaller box crowds molecules together (high "
                 "concentration, more frequent collisions, higher effective pressure); "
                 "a larger box dilutes them. Set it large enough that the starting "
                 "molecules don't overlap.")
        _conc_M = 1.0 / (6.022e23 * (box_L * 1e-10) ** 3 * 1000.0)   # mol/L per molecule
        st.caption(f"Volume ≈ {box_L**3:,.0f} Å³  ·  ~{_conc_M:.2f} M per molecule")
    with c2:
        solvent = st.selectbox(
            "Reaction medium (solvent charge)",
            ["Vacuum", "Water", "Ethanol", "Benzene"],
            help="What fills the reactor around the reactants. **Vacuum** = gas-phase / "
                 "neat (no medium). Choosing **Water, Ethanol or Benzene** immerses the "
                 "solutes in a Psi4 PCM continuum — an implicit solvent with that "
                 "fluid's dielectric constant — so polarity and solvation stabilization "
                 "affect the energetics, as a real solvent charge would. "
                 "Cost: PCM is ~50× slower per step than vacuum.")
    with c3:
        temperature = st.slider(
            "Operating temperature (K)", 50, 1500, 350, 10,
            help="The reactor setpoint. It seeds the initial Maxwell–Boltzmann velocity "
                 "distribution and is the temperature the Langevin thermostat holds the "
                 "system at throughout the run. Higher T = faster molecular motion and "
                 "more collision energy to climb reaction barriers (the Arrhenius "
                 "effect): hotter reactors react faster.")
        st.caption(f"≈ {temperature - 273.15:.0f} °C  ·  {temperature * 9/5 - 459.67:.0f} °F")
    with c4:
        n_steps = st.slider(
            "Reaction run length — steps", 5, 200, 30, 5,
            help="Number of molecular-dynamics integration steps = how long the reactor "
                 "is run. Total simulated time ≈ steps × timestep (timestep is in "
                 "Advanced). More steps follow the reaction longer but cost "
                 "proportionally more compute. Note: ab-initio MD reaches only "
                 "femtoseconds–picoseconds of real time — enough to watch bonds form / "
                 "break and molecules collide, not minutes-scale reactor residence.")

    with st.expander("⚙️ Advanced (level of theory & integrator)"):
        a1, a2, a3, a4 = st.columns(4)
        method = a1.selectbox(
            "Method", ["hf", "b3lyp", "pbe", "mp2"], index=0,
            help="Quantum-chemistry model that computes the forces driving the reaction. "
                 "**hf** (Hartree–Fock) is cheapest; **b3lyp / pbe** are DFT — better "
                 "chemistry at ~2× cost; **mp2** adds electron correlation, most "
                 "accurate and most expensive.")
        basis = a2.selectbox(
            "Basis set", ["sto-3g", "6-31g", "6-31g*", "cc-pvdz"], index=0,
            help="Size of the wavefunction expansion per atom. **sto-3g** is minimal and "
                 "fast (qualitative); **6-31g* / cc-pvdz** are larger and more accurate "
                 "but several× slower. Bigger basis = sharper energetics, longer runs.")
        timestep = a3.number_input(
            "Timestep (fs)", 0.1, 2.0, 0.5, 0.1,
            help="How much real time each MD step advances, in femtoseconds. Smaller is "
                 "more stable and accurate — important for light H atoms and high "
                 "temperature — while larger covers more time per step but can blow up "
                 "the integration. 0.5 fs is a safe default.")
        friction = a4.number_input(
            "Langevin friction (1/fs)", 0.001, 0.2, 0.02, 0.001, format="%.3f",
            help="How tightly the thermostat couples the system to the heat bath "
                 "(the reactor walls/jacket). Higher friction locks temperature to the "
                 "setpoint faster but makes motion more diffusive/viscous; lower "
                 "friction is more ballistic, closer to constant-energy dynamics.")
    st.caption(
        f"⏱️ Total simulated time ≈ {n_steps} steps × {timestep:g} fs = "
        f"**{n_steps * timestep:g} fs** ({n_steps * timestep / 1000:.3f} ps)"
    )

    # ── live system preview + runtime estimate ──────────────────────────────────
    combined = None
    group_sizes, charge, mult = None, 0, 1
    if mol_atoms:
        combined, group_sizes, min_center_d = assemble_system(mol_atoms, lane_used, box_L)
        n_atoms = len(combined)
        charge, mult = auto_multiplicity(combined)

        st.divider()
        st.subheader("3 · Assembled system")
        pcol, icol = st.columns([3, 1])
        with pcol:
            st.plotly_chart(
                system_3d_fig(combined, group_sizes,
                              title=f"{n_atoms} atoms · {combined.get_chemical_formula()}"),
                use_container_width=True)
            legend = "  ·  ".join(
                f'<span style="color:{LANE_COLORS[L]};">●</span> {L} ({resolved[L]["formula"]})'
                for L in lane_used)
            st.markdown(legend, unsafe_allow_html=True)
        with icol:
            st.metric("Total atoms", n_atoms)
            st.metric("Net charge", f"{charge:+d}")
            st.metric("Multiplicity", mult)
            if min_center_d < 3.0 and len(mol_atoms) > 1:
                st.warning(f"Molecule centers only {min_center_d:.1f} Å apart — "
                           "increase box size to avoid overlap.")

        # crude runtime heuristic from the calibration runs
        per_step = max(0.15, 0.02 * n_atoms ** 1.6)          # vacuum HF/STO-3G
        if solvent != "Vacuum":
            per_step *= 55                                   # PCM penalty
        if basis != "sto-3g":
            per_step *= 4
        if method != "hf":
            per_step *= 2
        est = per_step * n_steps
        est_txt = f"{est:.0f} s" if est < 90 else f"{est/60:.1f} min"
        (st.info if est < 120 else st.warning)(
            f"Rough estimate: ~{per_step:.1f} s/step × {n_steps} steps ≈ **{est_txt}** "
            f"of ab-initio MD." + (" PCM solvation dominates the cost."
                                   if solvent != "Vacuum" else ""))
    else:
        st.info("Resolve at least one molecule to assemble a system.")

    st.divider()

    # ── run ──────────────────────────────────────────────────────────────────────
    st.subheader("4 · Run")

    # MD runs in a background thread (see chemmodel.md.langevin). The script thread
    # only starts/stops it and renders state from session_state — so Stop is a
    # normal button on a responsive thread and a 22 s PCM step never freezes the UI.
    ss = st.session_state
    ss.setdefault("md", dict(
        phase="idle", k=0, n=0, t_start=None, current_T=None, s_per_step=0.0,
        engine_note="", stopped_at=None, error_text=None, target_T=None,
        t=[], pe=[], ke=[], T=[], final=None))
    ss.setdefault("md_queue", None)
    ss.setdefault("md_cancel", None)
    ss.setdefault("md_worker", None)

    # snapshot the current inputs every run (cheap, no copy) so the Run callback
    # — which fires BEFORE the next script body — sees the latest widget values
    ss["_md_inputs"] = dict(
        combined=combined, group_sizes=group_sizes, method=method, basis=basis,
        solvent=solvent, charge=charge, mult=mult, T=temperature, dt=timestep,
        fric=friction, n_steps=n_steps, psi4_ok=PSI4_OK)

    def _start_md():
        global _MD_WORKER
        inp = ss.get("_md_inputs")
        if inp is None or inp["combined"] is None:
            return
        # cross-session safety: cancel/reap any worker still alive process-wide
        if _MD_WORKER is not None and _MD_WORKER.is_alive():
            if ss.md_cancel is not None:
                ss.md_cancel.set()
            _MD_WORKER.join(timeout=2.0)
        if ss.md_worker is not None and not ss.md_worker.is_alive():
            ss.md_worker = None

        m_, b_, solv_ = inp["method"], inp["basis"], inp["solvent"]
        chg_, mul_, ok_ = inp["charge"], inp["mult"], inp["psi4_ok"]

        def make_calc():
            if ok_:
                return make_psi4_calculator(m_, b_, solv_, chg_, mul_)
            from ase.calculators.lj import LennardJones
            return LennardJones(epsilon=0.0103, sigma=3.4, rc=10.0)

        if ok_:
            engine_note = f"Psi4 {m_}/{b_}" + (
                f" + PCM({solv_})" if solv_ != "Vacuum" else " (vacuum)")
        else:
            engine_note = "Lennard-Jones demo force field (Psi4 unavailable)"

        params = dict(
            atoms=inp["combined"].copy(), group_sizes=inp["group_sizes"],
            make_calc=make_calc, psi4_ok=ok_, T=inp["T"], dt=inp["dt"],
            fric=inp["fric"], n_steps=inp["n_steps"])

        md = ss.md
        md.update(phase="running", k=0, n=inp["n_steps"], t_start=time.time(),
                  current_T=None, s_per_step=0.0, stopped_at=None,
                  error_text=None, t=[], pe=[], ke=[], T=[], final=None,
                  engine_note=engine_note, target_T=inp["T"])
        ss.md_queue = queue.Queue()
        ss.md_cancel = threading.Event()
        w = threading.Thread(target=md_worker,
                             args=(params, ss.md_queue, ss.md_cancel), daemon=True)
        ss.md_worker = w
        _MD_WORKER = w
        w.start()

    def _request_stop():
        if ss.md_cancel is not None:
            ss.md_cancel.set()
        ss.md["phase"] = "stopping"

    md = ss.md
    c_run, c_stop = st.columns([2, 1])
    c_run.button("🚀 Run Advanced Simulation", type="primary",
                 use_container_width=True,
                 disabled=(combined is None) or md["phase"] in ("running", "stopping"),
                 on_click=_start_md)
    c_stop.button("⏹ Stop", use_container_width=True,
                  disabled=md["phase"] not in ("running", "stopping"),
                  on_click=_request_stop)

    # ── live heartbeat: drains the worker queue and repaints state every 0.5 s ──
    run_every = "0.5s" if md["phase"] in ("running", "stopping") else None

    @st.fragment(run_every=run_every)
    def md_live():
        md = ss.md
        q = ss.md_queue
        just_finished = False
        if q is not None:
            while True:
                try:
                    kind, payload = q.get_nowait()
                except queue.Empty:
                    break
                if kind == "step":
                    md["k"] = payload["i"] + 1
                    md["current_T"] = payload["T"]
                    md["t"].append(payload["t"])
                    md["pe"].append(payload["pe"])
                    md["ke"].append(payload["ke"])
                    md["T"].append(payload["T"])
                else:
                    md["phase"] = {"done": "done", "stopped": "stopped",
                                   "error": "error"}[kind]
                    if kind == "stopped":
                        md["stopped_at"] = payload
                    elif kind == "error":
                        md["error_text"] = payload
                        md["stopped_at"] = md["k"]
                    elif kind == "done":
                        md["final"] = payload
                    if ss.md_worker is not None:
                        ss.md_worker.join(timeout=1.0)
                    ss.md_worker = None
                    ss.md_queue = None
                    ss.md_cancel = None
                    just_finished = True

        phase = md["phase"]
        if phase == "idle":
            return

        elapsed = (time.time() - md["t_start"]) if md["t_start"] else 0.0
        md["s_per_step"] = elapsed / max(md["k"], 1)
        label, state = {
            "running":  (f"Running ab-initio MD — step {md['k']}/{md['n']}", "running"),
            "stopping": ("Stopping after the current step (a PCM step can take ~20 s)…", "running"),
            "done":     (f"Completed {md['k']}/{md['n']} steps", "complete"),
            "stopped":  (f"Stopped at step {md['stopped_at']}/{md['n']}", "complete"),
            "error":    (f"Failed at step {md['stopped_at']}/{md['n']}", "error"),
        }[phase]
        with st.status(label, state=state,
                       expanded=phase in ("running", "stopping", "error")):
            st.caption(f"Engine: **{md['engine_note']}**  ·  "
                       f"Langevin thermostat @ {md['target_T']} K")
            st.progress(min(md["k"] / max(md["n"], 1), 1.0),
                        text=f"Step {md['k']}/{md['n']}")
            a, b, c, d = st.columns(4)
            a.metric("Step", f"{md['k']}/{md['n']}")
            b.metric("Elapsed", f"{elapsed:.0f} s")
            c.metric("s / step", f"{md['s_per_step']:.1f}")
            d.metric("Current T",
                     f"{md['current_T']:.0f} K" if md["current_T"] is not None else "—")
            if md["s_per_step"] > 5 and phase in ("running", "stopping"):
                st.caption("⏳ A PCM gradient can take ~20 s — the run is alive "
                           "between updates.")
            if phase == "error" and md["error_text"]:
                st.error(md["error_text"])

        if md["t"]:
            st.plotly_chart(
                md_charts_fig(md["t"], md["pe"], md["ke"], md["T"], md["target_T"]),
                use_container_width=True, key="md_live_chart")

        # one full rerun on completion so run_every re-evaluates to None and the
        # 0.5 s polling actually stops (a fragment auto-rerun cannot re-read it)
        if just_finished:
            st.rerun()

    md_live()

    # ── partial / final results (rendered from session_state, survives reruns) ──
    md = ss.md
    if md["k"] > 0 and md["phase"] in ("done", "stopped", "error"):
        head = {
            "done":    f"Completed {md['k']} MD steps.",
            "stopped": f"Stopped at step {md['stopped_at']} of {md['n']} — partial results below.",
            "error":   f"Failed at step {md['stopped_at']} of {md['n']} — partial results below.",
        }[md["phase"]]
        (st.success if md["phase"] == "done" else st.warning)(head)

        if md["pe"]:
            m1, m2, m3, m4 = st.columns(4)
            m1.metric("Final ΔPE (eV)", f"{md['pe'][-1]:+.4f}")
            m2.metric("Mean T (K)", f"{np.mean(md['T']):.0f}")
            m3.metric("Mean KE (eV)", f"{np.mean(md['ke']):.4f}")
            m4.metric("T drift (K)", f"{md['T'][-1] - md['T'][0]:+.0f}")

        if md["final"]:
            from ase import Atoms as _Atoms
            fin = md["final"]
            final_atoms = _Atoms(symbols=fin["symbols"], positions=fin["positions"])
            st.subheader("Final geometry")
            st.plotly_chart(
                system_3d_fig(final_atoms, fin["group_sizes"], height=360,
                              title=f"After {md['k']} steps"),
                use_container_width=True)

        with st.expander("Trajectory data"):
            import pandas as pd
            st.dataframe(pd.DataFrame({
                "t (fs)": md["t"], "ΔPE (eV)": md["pe"],
                "KE (eV)": md["ke"], "T (K)": md["T"],
            }), use_container_width=True)


# ══════════════════════════════════════════════════════════════════════════════
# TAB 2 — Reactor Engineering (ideal Batch / CSTR / PFR)
# ══════════════════════════════════════════════════════════════════════════════

with tab_reactor:
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots

    st.subheader("Ideal reactor design — A → products")
    st.caption(
        "Closed-form kinetics for a single **nth-order** reaction, rate law "
        "−rₐ = k·Cₐⁿ, in an ideal **Batch**, **CSTR** (continuous stirred tank), or "
        "**PFR** (plug flow) reactor. This bridges the molecular tab to reactor "
        "scale: a barrier you measure on a potential-energy surface (in eV) becomes "
        "an Arrhenius activation energy, which sets the rate constant k, which sizes "
        "the reactor. Runs instantly — no Psi4 needed."
    )

    rc1, rc2, rc3 = st.columns(3)
    with rc1:
        reactor_type = st.selectbox(
            "Reactor type", ["Batch", "CSTR", "PFR"],
            help="**Batch**: closed vessel, composition evolves in time (conversion vs "
                 "time). **CSTR**: continuous, perfectly mixed, operates at the exit "
                 "composition (needs the most volume for a given conversion at positive "
                 "order). **PFR**: continuous, no back-mixing, composition changes down "
                 "the tube (most volume-efficient). For CSTR/PFR the x-axis is space-time "
                 "τ = V/v₀.")
        order = st.selectbox(
            "Reaction order n", [0, 1, 2], index=1,
            help="Exponent in −rₐ = k·Cₐⁿ. **0**: rate independent of concentration "
                 "(e.g. saturated catalyst). **1**: rate ∝ Cₐ (most common, e.g. "
                 "radioactive-style decay, many decompositions). **2**: rate ∝ Cₐ² "
                 "(bimolecular A+A). The order fixes the units of k.")
    with rc2:
        CA0 = st.number_input(
            "Inlet concentration Cₐ₀ (mol/L)", 0.001, 100.0, 1.0, 0.1,
            help="Initial (batch) or feed (continuous) concentration of reactant A. "
                 "Sets the scale of the rate and, for orders ≠ 1, enters the design "
                 "equations directly.")
        k_mode = st.radio(
            "Rate-constant source", ["Direct k", "Arrhenius k(T)"], horizontal=True,
            help="**Direct k**: type the rate constant outright. **Arrhenius k(T)**: "
                 "compute k = A·exp(−Eₐ/RT) from a pre-exponential factor A and an "
                 "activation energy Eₐ at the operating temperature — this is where a "
                 "PES barrier in eV plugs in.")
    with rc3:
        tempF = st.number_input(
            "Operating temperature (°F)", -50.0, 1500.0, 140.0, 5.0,
            help="Reactor operating temperature (US-customary °F to match the plant "
                 "reactor spec sheet; default 140 °F mirrors reactor_specs_template.csv). "
                 "Used by the Arrhenius equation to evaluate k(T).")
        T_K = (tempF - 32.0) * 5.0 / 9.0 + 273.15
        st.caption(f"= {(tempF - 32.0) * 5.0 / 9.0:.1f} °C = **{T_K:.1f} K**")

    # ── rate constant ─────────────────────────────────────────────────────────
    st.markdown("##### Kinetics")
    if k_mode.startswith("Direct"):
        default_k = {0: 0.10, 1: 0.05, 2: 0.02}[order]
        k = st.number_input(
            f"Rate constant k  ({k_units(order)})",
            min_value=1e-9, value=float(default_k), format="%.4g",
            help=f"For order {order}, k carries units of {k_units(order)}. Larger k = "
                 "faster reaction = smaller reactor / shorter batch time.")
        Ea_kJ = A_pre = None
    else:
        ac1, ac2, ac3 = st.columns(3)
        A_pre = ac1.number_input(
            f"Pre-exponential A  ({k_units(order)})", min_value=1e-3,
            value=1.0e13, format="%.3e",
            help="Arrhenius frequency factor — the rate constant in the high-temperature "
                 "limit (collision/attempt frequency). Typical first-order values are "
                 "~10¹²–10¹⁴ s⁻¹.")
        ea_unit = ac2.selectbox(
            "Eₐ units", ["kJ/mol", "kcal/mol", "eV"],
            help="Pick eV to enter an activation barrier straight from a potential-energy "
                 "surface scan; it is converted internally (1 eV = 96.49 kJ/mol).")
        ea_val = ac3.number_input(
            f"Activation energy Eₐ  ({ea_unit})", min_value=0.0,
            value={"kJ/mol": 50.0, "kcal/mol": 12.0, "eV": 0.5}[ea_unit], step=0.5,
            help="Energy barrier the reaction must climb. From the molecular tab / a PES "
                 "scan this is the height of the transition state above the reactants.")
        Ea_kJ = {"kJ/mol": ea_val, "kcal/mol": ea_val * 4.184,
                 "eV": ea_val * EV_TO_KJMOL}[ea_unit]
        k = float(arrhenius_k(A_pre, Ea_kJ, T_K))
        st.caption(
            f"Eₐ = **{Ea_kJ:.1f} kJ/mol** = {Ea_kJ / EV_TO_KJMOL:.3f} eV "
            f"= {Ea_kJ / 4.184:.1f} kcal/mol  →  k(T) = **{k:.4e} {k_units(order)}**"
        )

    st.markdown(f"**Active rate constant:** k = `{k:.4e}` {k_units(order)} at {T_K:.0f} K")
    st.divider()

    # ── operating / sizing controls ────────────────────────────────────────────
    oc1, oc2 = st.columns(2)
    with oc1:
        Xtarget = st.slider(
            "Target conversion X", 0.01, 0.99, 0.80, 0.01,
            help="Fraction of A consumed you want the reactor to achieve. Drives the "
                 "required batch time (Batch) or space-time/volume (CSTR, PFR).")
    with oc2:
        if reactor_type == "Batch":
            t_required = batch_time_for_X(CA0, k, order, Xtarget)
            t_max = st.slider(
                "Batch time horizon (s)", 1.0, float(max(10.0, 3.0 * t_required)),
                float(max(5.0, 1.5 * t_required)), step=1.0,
                help="How long to plot the batch trajectory. Defaults to ~1.5× the time "
                     "needed to hit the target conversion.")
        else:
            v0 = st.number_input(
                "Volumetric feed rate v₀ (L/s)", min_value=1e-4, value=1.0, format="%.4g",
                help="Feed flow into the continuous reactor. With the computed space-time "
                     "τ it sets the required reactor volume V = τ·v₀.")

    # ════════════════════════════════════════════════════════════════════════════
    # Batch reactor
    # ════════════════════════════════════════════════════════════════════════════
    if reactor_type == "Batch":
        t_grid = np.linspace(0.0, t_max, 400)
        CA, X = batch_profile(CA0, k, order, t_grid)
        t_req = batch_time_for_X(CA0, k, order, Xtarget)

        m1, m2, m3 = st.columns(3)
        m1.metric("Time to target X", f"{t_req:.3g} s",
                  help=f"Batch time to reach X = {Xtarget:.2f}")
        m2.metric("Conversion at horizon", f"{X[-1]*100:.1f} %")
        m3.metric("C_A at horizon", f"{CA[-1]:.4g} mol/L")

        fig = make_subplots(specs=[[{"secondary_y": True}]])
        fig.add_trace(go.Scatter(x=t_grid, y=X, name="Conversion X",
                                 line=dict(color="#5b9bd5", width=3)),
                      secondary_y=False)
        fig.add_trace(go.Scatter(x=t_grid, y=CA, name="C_A (mol/L)",
                                 line=dict(color="#70ad47", width=3, dash="dot")),
                      secondary_y=True)
        if t_req <= t_max:
            fig.add_vline(x=t_req, line_dash="dash", line_color="#e74c3c",
                          annotation_text=f"X={Xtarget:.2f} @ {t_req:.3g}s")
        fig.update_xaxes(title_text="Time (s)")
        fig.update_yaxes(title_text="Conversion X", range=[0, 1], secondary_y=False)
        fig.update_yaxes(title_text="C_A (mol/L)", secondary_y=True)
        fig.update_layout(height=460, template="plotly_dark",
                          paper_bgcolor="#1a1a2e", plot_bgcolor="#1a1a2e",
                          font_color="#ecf0f1", title="Batch reactor — conversion & concentration vs time",
                          legend=dict(orientation="h", y=1.08))
        st.plotly_chart(fig, use_container_width=True)

    # ════════════════════════════════════════════════════════════════════════════
    # CSTR / PFR
    # ════════════════════════════════════════════════════════════════════════════
    else:
        Xg = np.linspace(0.0, 0.99, 400)
        tau_p = tau_pfr(CA0, k, order, Xg)
        tau_c = tau_cstr(CA0, k, order, Xg)
        tau_p[0] = tau_c[0] = 0.0

        tau_sel = (tau_pfr if reactor_type == "PFR" else tau_cstr)(CA0, k, order, Xtarget)
        tau_pfr_t = float(tau_pfr(CA0, k, order, Xtarget))
        tau_cstr_t = float(tau_cstr(CA0, k, order, Xtarget))
        V_sel = tau_sel * v0
        V_pfr = tau_pfr_t * v0
        V_cstr = tau_cstr_t * v0

        m1, m2, m3, m4 = st.columns(4)
        m1.metric(f"Space-time τ ({reactor_type})", f"{tau_sel:.3g} s",
                  help="τ = V/v₀ required to reach the target conversion.")
        m2.metric("Required volume V", f"{V_sel:.3g} L",
                  help=f"V = τ·v₀ = {tau_sel:.3g}·{v0:.3g}")
        m3.metric("…in ft³", f"{V_sel * L_TO_FT3:.3g} ft³")
        m4.metric("CSTR/PFR volume ratio", f"{(V_cstr / V_pfr):.2f}×",
                  help="How much larger a CSTR must be than a PFR for the same "
                       "conversion (≥1 for positive order).")

        # X vs τ, both reactors
        fig1 = go.Figure()
        fig1.add_trace(go.Scatter(x=tau_p, y=Xg, name="PFR",
                                  line=dict(color="#5b9bd5", width=3)))
        fig1.add_trace(go.Scatter(x=tau_c, y=Xg, name="CSTR",
                                  line=dict(color="#e67e22", width=3)))
        fig1.add_hline(y=Xtarget, line_dash="dot", line_color="#888",
                       annotation_text=f"target X={Xtarget:.2f}")
        fig1.add_vline(x=tau_sel, line_dash="dash", line_color="#e74c3c",
                       annotation_text=f"τ={tau_sel:.3g}s ({reactor_type})")
        fig1.update_xaxes(title_text="Space-time τ = V/v₀ (s)",
                          range=[0, float(min(np.nanmax(tau_c[Xg <= 0.99]),
                                              5 * max(tau_pfr_t, 1e-9)))])
        fig1.update_yaxes(title_text="Conversion X", range=[0, 1])
        fig1.update_layout(height=420, template="plotly_dark",
                           paper_bgcolor="#1a1a2e", plot_bgcolor="#1a1a2e",
                           font_color="#ecf0f1", legend=dict(orientation="h", y=1.1),
                           title="Conversion vs space-time — PFR is more volume-efficient")
        st.plotly_chart(fig1, use_container_width=True)

        # Levenspiel plot: 1/(-rA) vs X with the reactor's "area" shaded
        Xl = np.linspace(0.0, Xtarget, 300)
        inv = inv_rate(CA0, k, order, Xl)
        fig2 = go.Figure()
        fig2.add_trace(go.Scatter(x=Xl, y=inv, name="1/(−rₐ)",
                                  line=dict(color="#9b59b6", width=3)))
        if reactor_type == "PFR":
            fig2.add_trace(go.Scatter(
                x=np.r_[Xl, Xl[::-1]],
                y=np.r_[inv, np.zeros_like(inv)],
                fill="toself", fillcolor="rgba(91,155,213,0.25)",
                line=dict(width=0), name="τ = ∫ dX/(−rₐ)  (area)"))
        else:  # CSTR rectangle: height = 1/(-rA at exit), width = Xtarget
            inv_exit = float(inv_rate(CA0, k, order, Xtarget))
            fig2.add_trace(go.Scatter(
                x=[0, Xtarget, Xtarget, 0, 0],
                y=[inv_exit, inv_exit, 0, 0, inv_exit],
                fill="toself", fillcolor="rgba(230,126,34,0.25)",
                line=dict(color="#e67e22", width=1), name="τ = X·[1/(−rₐ)]_exit  (rectangle)"))
        fig2.update_xaxes(title_text="Conversion X", range=[0, 1])
        fig2.update_yaxes(title_text="1 / (−rₐ)   (L·s/mol)")
        fig2.update_layout(height=400, template="plotly_dark",
                           paper_bgcolor="#1a1a2e", plot_bgcolor="#1a1a2e",
                           font_color="#ecf0f1", legend=dict(orientation="h", y=1.12),
                           title=f"Levenspiel plot — shaded area = τ for the {reactor_type}")
        st.plotly_chart(fig2, use_container_width=True)

    # ── Arrhenius plot (only meaningful in Arrhenius mode) ─────────────────────
    if not k_mode.startswith("Direct"):
        with st.expander("📈 Arrhenius behaviour — k(T)", expanded=False):
            T_sweep = np.linspace(max(250.0, T_K - 150.0), T_K + 200.0, 200)
            k_sweep = arrhenius_k(A_pre, Ea_kJ, T_sweep)
            af = make_subplots(rows=1, cols=2,
                               subplot_titles=("k vs T", "ln k vs 1000/T (linear)"))
            af.add_trace(go.Scatter(x=T_sweep, y=k_sweep, line=dict(color="#1abc9c", width=3),
                                    name="k(T)"), row=1, col=1)
            af.add_trace(go.Scatter(x=1000.0 / T_sweep, y=np.log(k_sweep),
                                    line=dict(color="#f1c40f", width=3), name="ln k"),
                         row=1, col=2)
            af.add_vline(x=T_K, line_dash="dash", line_color="#e74c3c", row=1, col=1)
            af.update_xaxes(title_text="T (K)", row=1, col=1)
            af.update_yaxes(title_text=f"k ({k_units(order)})", type="log", row=1, col=1)
            af.update_xaxes(title_text="1000/T (1/K)", row=1, col=2)
            af.update_yaxes(title_text="ln k", row=1, col=2)
            af.update_layout(height=360, template="plotly_dark", showlegend=False,
                             paper_bgcolor="#1a1a2e", plot_bgcolor="#1a1a2e",
                             font_color="#ecf0f1",
                             title=f"Slope of ln k vs 1/T = −Eₐ/R  (Eₐ = {Ea_kJ:.1f} kJ/mol)")
            st.plotly_chart(af, use_container_width=True)

    # ── optional adiabatic temperature rise ────────────────────────────────────
    with st.expander("🔥 Adiabatic temperature rise (energy balance)", expanded=False):
        st.caption("For an adiabatic reactor, the heat released by reaction raises the "
                   "mixture temperature: ΔT_ad = (−ΔH_rxn · Cₐ₀ · X) / (ρ · C_p).")
        h1, h2, h3 = st.columns(3)
        dH = h1.number_input("ΔH_rxn (kJ/mol)", value=-80.0, step=5.0,
                             help="Heat of reaction per mole of A. Negative = exothermic "
                                  "(releases heat, temperature rises).")
        rho = h2.number_input("Mixture density ρ (kg/L)", min_value=0.01, value=1.0, step=0.05)
        Cp = h3.number_input("Heat capacity C_p (kJ/kg·K)", min_value=0.01, value=4.18, step=0.1,
                             help="Default 4.18 ≈ water.")
        dT_ad = (-dH * CA0 * Xtarget) / (rho * Cp)
        T_out = T_K + dT_ad
        e1, e2 = st.columns(2)
        e1.metric("Adiabatic ΔT at target X", f"{dT_ad:+.1f} K")
        e2.metric("Outlet temperature", f"{T_out:.1f} K  ({(T_out-273.15)*9/5+32:.0f} °F)")
        if dT_ad > 50:
            st.warning("Large exotherm — a real reactor would need cooling / jacket duty "
                       "to stay near setpoint and avoid runaway.")
