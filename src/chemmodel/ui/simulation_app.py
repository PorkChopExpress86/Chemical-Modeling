"""
Chemical & Polymer Simulation Web App — Streamlit frontend.

Thin presentation layer over the ``chemmodel`` domain package:
- ASE as the molecular data manager (``chemmodel.chemistry``)
- Quantum-chemistry engines via the registry (``chemmodel.engines``)
- Catalytic-reaction potential-energy-surface scans (``chemmodel.reactions``)

Run with:  streamlit run src/chemmodel/ui/simulation_app.py
"""

from __future__ import annotations

import traceback

import numpy as np
import streamlit as st

from chemmodel.chemistry import build_reaction_system, xyz_to_atoms
from chemmodel.chemistry.presets import CATALYSTS, PRESET_REACTIONS, PRESETS
from chemmodel.engines import available_engines, get_engine
from chemmodel.reactions import decompose_energies, run_reaction_scan
from chemmodel.ui.theme import SIMULATION_CSS
from chemmodel.ui.viz import (
    mol_3d_fig,
    render_reaction_dashboard,
    render_single_dashboard,
)

_AVAIL = available_engines()
PSI4_OK = _AVAIL["psi4"]
PYSCF_OK = _AVAIL["pyscf"]
CP2K_OK = _AVAIL["cp2k"]


def run_single_molecule(engine_name: str, atoms, **params) -> dict:
    """Single-point energy + forces via the chosen engine, demo fallback if absent."""
    eng = get_engine(engine_name)
    if eng.available:
        er = eng.single_point(atoms, with_forces=True, **params)
    else:
        short = eng.label.split(" —")[0]
        er = get_engine("demo").single_point(
            atoms, with_forces=True, label=f"{short} (demo — not installed)", **params
        )
    return {"energy_eV": er.energy_eV, "forces": er.forces, "note": er.note}


# ═══════════════════════════════════════════════════════════════════════════════
# Streamlit page
# ═══════════════════════════════════════════════════════════════════════════════

st.set_page_config(
    page_title="Chemical & Polymer Simulator",
    page_icon="⚗️",
    layout="wide",
)

st.markdown(SIMULATION_CSS, unsafe_allow_html=True)

st.title("⚗️ Chemical & Polymer Simulation Suite")
st.caption("ASE · Psi4 · PySCF · CP2K  |  Single Molecule & Catalytic Reaction")

tab_single, tab_reaction = st.tabs(["🔬 Single Molecule", "⚗️ Catalytic Reaction"])

# ═══════════════════════════════════════════════════════════════════════════════
# TAB 1 — Single Molecule
# ═══════════════════════════════════════════════════════════════════════════════

with tab_single:
    with st.sidebar:
        st.header("Single Molecule")
        input_mode = st.radio("Input mode", ["Preset", "Raw XYZ"], horizontal=True)
        if input_mode == "Preset":
            preset_name = st.selectbox("Molecule", list(PRESETS.keys()))
            xyz_text = PRESETS[preset_name]
            st.code(xyz_text, language="text")
        else:
            xyz_text = st.text_area("XYZ coordinates", height=200,
                                    value=PRESETS["Water (H₂O)"])

        st.divider()
        engine_s = st.selectbox(
            "Simulation Engine",
            ["Psi4 — Quantum Chemistry", "PySCF — QC + Periodic BC", "CP2K — MD Forces"],
        )

        st.divider()
        if engine_s.startswith("Psi4"):
            st.subheader("Psi4 Parameters")
            p4_method = st.selectbox("Method", ["B3LYP", "HF", "MP2", "CCSD", "PBE"])
            p4_basis = st.selectbox("Basis Set", ["6-31G*", "cc-pVDZ", "cc-pVTZ", "def2-TZVP"])
            p4_maxiter = st.slider("Max SCF iterations", 50, 500, 150, 10)
            bc = "#5d3a9b" if PSI4_OK else "#7a3b00"
            st.markdown(f'<span class="engine-badge" style="background:{bc};">'
                        f'{"✅ psi4" if PSI4_OK else "⚠️ demo mode"}</span>',
                        unsafe_allow_html=True)
        elif engine_s.startswith("PySCF"):
            st.subheader("PySCF Parameters")
            py_method = st.selectbox("Method", ["HF", "RHF", "UHF", "B3LYP", "PBE"])
            py_basis = st.selectbox("Basis Set", ["sto-3g", "6-31g", "6-31g*", "cc-pvdz", "def2-svp"])
            use_pbc = st.toggle("Periodic Boundary Conditions", value=False)
            cell_size = st.slider("Cell size (Å)", 5.0, 30.0, 15.0, 0.5, disabled=not use_pbc)
            kpts = st.select_slider("k-grid", [1, 2, 3, 4], value=1, disabled=not use_pbc)
            st.markdown(f'<span class="engine-badge">{"✅ pyscf" if PYSCF_OK else "⚠️ demo"}</span>',
                        unsafe_allow_html=True)
        else:
            st.subheader("CP2K Parameters")
            cp2k_ts = st.number_input("Timestep (fs)", 0.1, 5.0, 0.5, 0.1)
            cp2k_st = st.number_input("MD Steps", 1, 10_000, 100, 10)
            cp2k_cut = st.slider("PW Cutoff (Ry)", 100, 1000, 400, 50)
            cp2k_rel = st.slider("Rel Cutoff (Ry)", 20, 100, 60, 10)
            st.markdown(f'<span class="engine-badge">{"✅ cp2k" if CP2K_OK else "⚠️ demo"}</span>',
                        unsafe_allow_html=True)

    # ── molecule preview ─────────────────────────────────────────────────────
    col_mol, col_run = st.columns([3, 1])
    with col_mol:
        st.subheader("Molecule Preview")
        try:
            atoms = xyz_to_atoms(xyz_text)
            st.plotly_chart(mol_3d_fig(atoms, height=320), use_container_width=True)
            st.caption(f"{len(atoms)} atoms · {atoms.get_chemical_formula()}")
        except Exception as exc:
            st.error(f"XYZ parse error: {exc}")
            atoms = None

    with col_run:
        st.subheader("Run")
        run_btn = st.button("▶ Run Simulation", type="primary", use_container_width=True)
        if atoms:
            st.metric("Atoms", len(atoms))
            st.metric("Formula", atoms.get_chemical_formula())

    st.divider()

    if run_btn:
        if atoms is None:
            st.error("Fix the XYZ input first.")
        else:
            status = st.status("Running…", expanded=True)
            with status:
                try:
                    if engine_s.startswith("Psi4"):
                        result = run_single_molecule(
                            "psi4", atoms, basis=p4_basis, method=p4_method, max_iter=p4_maxiter)
                    elif engine_s.startswith("PySCF"):
                        result = run_single_molecule(
                            "pyscf", atoms, basis=py_basis, method=py_method,
                            use_pbc=use_pbc, cell_size=cell_size, kpts=kpts)
                    else:
                        result = run_single_molecule(
                            "cp2k", atoms, basis="", method="",
                            cutoff=cp2k_cut, rel_cutoff=cp2k_rel)
                    st.write(f"✅ Energy: **{result['energy_eV']:.6f} eV**")
                    if result.get("note"):
                        st.info(result["note"])
                except Exception as exc:
                    status.update(label="Failed", state="error")
                    st.error(f"{type(exc).__name__}: {exc}")
                    with st.expander("Traceback"):
                        st.code(traceback.format_exc())
                    st.stop()
            status.update(label="Complete", state="complete")
            render_single_dashboard(result, atoms)
            with st.expander("Raw JSON"):
                st.json({"energy_eV": result["energy_eV"],
                         "forces": result["forces"].tolist(),
                         "note": result.get("note", "")})


# ═══════════════════════════════════════════════════════════════════════════════
# TAB 2 — Catalytic Reaction
# ═══════════════════════════════════════════════════════════════════════════════

with tab_reaction:
    st.subheader("Catalytic Reaction — Potential Energy Surface Scan")
    st.caption(
        "Select two reactants and an optional catalyst. "
        "The app builds the combined system, scans molecule B from "
        "a far separation toward A, and plots the energy profile."
    )

    # ── preset reaction loader ────────────────────────────────────────────────
    with st.expander("⚡ Load a preset reaction scenario", expanded=True):
        preset_rxn = st.selectbox("Preset reactions", ["— custom —"] + list(PRESET_REACTIONS.keys()))
        if preset_rxn != "— custom —":
            rxn = PRESET_REACTIONS[preset_rxn]
            st.info(f"Loaded: A={rxn['A']}, B={rxn['B']}, catalyst={rxn['cat']}, "
                    f"scan {rxn['start']}→{rxn['end']} Å, {rxn['n']} points")

    use_preset = preset_rxn != "— custom —"
    rxn_data = PRESET_REACTIONS.get(preset_rxn, {})

    # ── reactant / catalyst selectors ────────────────────────────────────────
    col_a, col_b, col_c = st.columns(3)

    with col_a:
        st.markdown("**Reactant A**")
        mol_names = list(PRESETS.keys())
        default_a = mol_names.index(rxn_data.get("A", mol_names[0])) if use_preset else 0
        rxn_a_name = st.selectbox("Molecule A", mol_names, index=default_a, key="rxn_a")
        rxn_a_custom = st.text_area("or paste XYZ for A", height=80, key="xyz_a",
                                    placeholder="Leave blank to use preset above")
        try:
            atoms_a = xyz_to_atoms(rxn_a_custom) if rxn_a_custom.strip() else xyz_to_atoms(PRESETS[rxn_a_name])
            st.caption(f"✅ {atoms_a.get_chemical_formula()} ({len(atoms_a)} atoms)")
        except Exception as e:
            st.error(str(e))
            atoms_a = None

    with col_b:
        st.markdown("**Reactant B**")
        default_b = mol_names.index(rxn_data.get("B", mol_names[1])) if use_preset and rxn_data.get("B") in mol_names else min(1, len(mol_names) - 1)
        rxn_b_name = st.selectbox("Molecule B", mol_names, index=default_b, key="rxn_b")
        rxn_b_custom = st.text_area("or paste XYZ for B", height=80, key="xyz_b",
                                    placeholder="Leave blank to use preset above")
        try:
            atoms_b = xyz_to_atoms(rxn_b_custom) if rxn_b_custom.strip() else xyz_to_atoms(PRESETS[rxn_b_name])
            st.caption(f"✅ {atoms_b.get_chemical_formula()} ({len(atoms_b)} atoms)")
        except Exception as e:
            st.error(str(e))
            atoms_b = None

    with col_c:
        st.markdown("**Catalyst**")
        cat_names = list(CATALYSTS.keys())
        default_cat = cat_names.index(rxn_data.get("cat", "None")) if use_preset else 0
        cat_name = st.selectbox("Catalyst", cat_names, index=default_cat, key="cat_sel")
        cat_info = CATALYSTS[cat_name]
        st.markdown(f'<span class="engine-badge" style="border-color:{cat_info["color"]};color:{cat_info["color"]};">'
                    f'{cat_info["label"]}</span>', unsafe_allow_html=True)
        if cat_info["xyz"]:
            try:
                catalyst_atoms = xyz_to_atoms(cat_info["xyz"])
                st.caption(f"✅ {catalyst_atoms.get_chemical_formula()}, charge={cat_info['charge']:+d}")
            except Exception as e:
                st.error(str(e))
                catalyst_atoms = None
        else:
            catalyst_atoms = None

    catalyst_charge = cat_info["charge"]

    # ── preview of initial arrangement ───────────────────────────────────────
    if atoms_a and atoms_b:
        prev_dist = rxn_data.get("start", 4.0) if use_preset else 4.0
        try:
            preview_sys = build_reaction_system(atoms_a, atoms_b, catalyst_atoms, prev_dist)
            nA = len(atoms_a)
            nB = len(atoms_b)
            nC = len(catalyst_atoms) if catalyst_atoms else 0
            st.subheader("Initial arrangement preview")
            st.plotly_chart(
                mol_3d_fig(preview_sys,
                           title=f"Starting configuration — separation {prev_dist:.1f} Å",
                           height=360, group_sizes=(nA, nB, nC)),
                use_container_width=True,
            )
            st.caption("🔵 A  ·  🟢 B  ·  🟡 Catalyst  |  B approaches A along z-axis")
        except Exception as e:
            st.warning(f"Preview failed: {e}")

    st.divider()

    # ── engine + scan parameters ──────────────────────────────────────────────
    col_eng, col_scan = st.columns([1, 1])

    with col_eng:
        st.subheader("Engine")
        rxn_engine = st.selectbox("Engine", ["PySCF", "Psi4"], key="rxn_engine",
                                  help="CP2K uses force fields, not recommended for bond-formation scans")
        if rxn_engine == "PySCF":
            rxn_basis = st.selectbox("Basis Set", ["sto-3g", "6-31g", "6-31g*", "cc-pvdz"], key="rxn_basis")
            rxn_method = st.selectbox("Method", ["HF", "B3LYP", "PBE"], key="rxn_method")
            ok_tag = "✅ pyscf" if PYSCF_OK else "⚠️ demo mode"
        else:
            rxn_basis = st.selectbox("Basis Set", ["sto-3g", "6-31G*", "cc-pVDZ"], key="rxn_basis2")
            rxn_method = st.selectbox("Method", ["HF", "B3LYP", "MP2"], key="rxn_method2")
            ok_tag = "✅ psi4" if PSI4_OK else "⚠️ demo mode"
        st.markdown(f'<span class="engine-badge">{ok_tag}</span>', unsafe_allow_html=True)

    with col_scan:
        st.subheader("Scan Parameters")
        scan_start = st.slider("Start separation (Å)", 2.5, 10.0,
                               float(rxn_data.get("start", 5.0)), 0.25, key="scan_start")
        scan_end = st.slider("End separation (Å)", 0.8, 4.0,
                             float(rxn_data.get("end", 1.5)), 0.1, key="scan_end")
        scan_n = st.slider("Number of points", 4, 20,
                           int(rxn_data.get("n", 10)), 1, key="scan_n")
        st.caption(f"Step size: {(scan_start - scan_end) / max(scan_n - 1, 1):.3f} Å")

    # ── run button ────────────────────────────────────────────────────────────
    rxn_run = st.button("▶ Run Reaction Scan", type="primary",
                        use_container_width=False, key="rxn_run")

    if rxn_run:
        if not (atoms_a and atoms_b):
            st.error("Fix molecule inputs before running.")
        elif scan_end >= scan_start:
            st.error("End separation must be less than start separation.")
        else:
            distances = np.linspace(scan_start, scan_end, scan_n)

            prog_bar = st.progress(0, text="Starting scan…")
            scan_status = st.status("Scanning reaction coordinate…", expanded=True)

            scan_results: list[dict] = []

            def _progress(done, total, d, e):
                pct = int(done / total * 100)
                prog_bar.progress(pct,
                    text=f"Point {done}/{total}  d={d:.2f} Å"
                         + (f"  ΔE={e:.4f} eV" if e is not None else ""))
                with scan_status:
                    tag = f"✅ d={d:.2f} Å  ΔE={e:.4f} eV" if e is not None else f"❌ d={d:.2f} Å  failed"
                    st.write(tag)

            with scan_status:
                try:
                    scan_results = run_reaction_scan(
                        atoms_a, atoms_b, catalyst_atoms, catalyst_charge,
                        distances, rxn_engine, rxn_basis, rxn_method,
                        progress_cb=_progress,
                    )
                    decomp = decompose_energies(
                        atoms_a, atoms_b, catalyst_atoms, catalyst_charge,
                        rxn_engine, rxn_basis, rxn_method,
                    )
                except Exception as exc:
                    scan_status.update(label="Scan failed", state="error")
                    st.error(f"{type(exc).__name__}: {exc}")
                    with st.expander("Traceback"):
                        st.code(traceback.format_exc())
                    st.stop()

            prog_bar.progress(100, text="Complete")
            scan_status.update(label="Scan complete", state="complete")

            st.subheader("Reaction Analysis")
            render_reaction_dashboard(
                scan_results, atoms_a, atoms_b, catalyst_atoms,
                decomp, cat_name, rxn_engine, rxn_basis, rxn_method,
            )
