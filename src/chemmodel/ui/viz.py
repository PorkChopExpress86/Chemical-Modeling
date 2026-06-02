"""Plotly figures and Streamlit result dashboards for the simulation apps."""

from __future__ import annotations

from typing import Any

import numpy as np
import streamlit as st
from ase import Atoms

from chemmodel.ui.theme import ATOM_RADII, GROUP_COLORS, LANE_COLORS, LANES


# ── 3-D molecule viewer (simulation app) ──────────────────────────────────────
def mol_3d_fig(atoms: Atoms, title: str = "", height: int = 320,
               group_sizes: tuple[int, ...] = ()):
    import plotly.graph_objects as go

    syms = atoms.get_chemical_symbols()
    pos = atoms.get_positions()
    n = len(syms)

    # colour atoms by group (A=blue, B=green, cat=gold)
    atom_group = [0] * n
    cursor = 0
    for gi, gs in enumerate(group_sizes):
        for k in range(cursor, cursor + gs):
            atom_group[k] = gi
        cursor += gs

    colors = [GROUP_COLORS[atom_group[i]] for i in range(n)]
    sizes = [ATOM_RADII.get(s, 12) for s in syms]

    fig = go.Figure()
    fig.add_trace(go.Scatter3d(
        x=pos[:, 0], y=pos[:, 1], z=pos[:, 2],
        mode="markers+text", text=syms, textposition="top center",
        marker=dict(size=sizes, color=colors, line=dict(width=1, color="#ffffff")),
        hovertemplate="<b>%{text}</b><br>x=%{x:.2f} y=%{y:.2f} z=%{z:.2f}<extra></extra>",
    ))
    fig.update_layout(
        height=height, title=title,
        margin=dict(l=0, r=0, t=30, b=0),
        scene=dict(bgcolor="#0f0f1a",
                   xaxis=dict(color="#888"), yaxis=dict(color="#888"),
                   zaxis=dict(color="#888")),
        paper_bgcolor="#0f0f1a", font_color="#ecf0f1",
        template="plotly_dark",
    )
    return fig


# ── 3-D system viewer (sandbox app) ───────────────────────────────────────────
def system_3d_fig(atoms, group_sizes, height=380, title=""):
    import plotly.graph_objects as go

    syms = atoms.get_chemical_symbols()
    pos = atoms.get_positions()
    n = len(syms)

    group_of = [0] * n
    cur = 0
    for gi, gs in enumerate(group_sizes):
        for k in range(cur, cur + gs):
            group_of[k] = gi
        cur += gs
    palette = [LANE_COLORS[L] for L in LANES]
    colors = [palette[group_of[i]] for i in range(n)]
    sizes = [ATOM_RADII.get(s, 12) for s in syms]

    fig = go.Figure(go.Scatter3d(
        x=pos[:, 0], y=pos[:, 1], z=pos[:, 2],
        mode="markers+text", text=syms, textposition="top center",
        marker=dict(size=sizes, color=colors, line=dict(width=1, color="#fff")),
        hovertemplate="<b>%{text}</b><br>%{x:.2f}, %{y:.2f}, %{z:.2f}<extra></extra>",
    ))
    fig.update_layout(
        height=height, title=title, margin=dict(l=0, r=0, t=30, b=0),
        scene=dict(bgcolor="#0f0f1a",
                   xaxis=dict(color="#888", title="X (Å)"),
                   yaxis=dict(color="#888", title="Y (Å)"),
                   zaxis=dict(color="#888", title="Z (Å)")),
        paper_bgcolor="#0f0f1a", font_color="#ecf0f1", template="plotly_dark",
    )
    return fig


# ── MD live charts (sandbox app) ──────────────────────────────────────────────
def md_charts_fig(t, pe_rel, ke, temp, target_T):
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots

    fig = make_subplots(
        rows=3, cols=1, shared_xaxes=True, vertical_spacing=0.07,
        subplot_titles=("Potential Energy ΔE (eV)",
                        "Kinetic Energy (eV)",
                        "Temperature (K)"),
    )
    fig.add_trace(go.Scatter(x=t, y=pe_rel, mode="lines+markers",
                             line=dict(color="#5b9bd5", width=2), name="ΔPE"),
                  row=1, col=1)
    fig.add_trace(go.Scatter(x=t, y=ke, mode="lines+markers",
                             line=dict(color="#70ad47", width=2), name="KE"),
                  row=2, col=1)
    fig.add_trace(go.Scatter(x=t, y=temp, mode="lines+markers",
                             line=dict(color="#e67e22", width=2), name="T"),
                  row=3, col=1)
    if target_T:
        fig.add_hline(y=target_T, line_dash="dot", line_color="#888",
                      row=3, col=1, annotation_text=f"target {target_T} K")
    fig.update_xaxes(title_text="Time (fs)", row=3, col=1)
    fig.update_layout(
        height=620, template="plotly_dark", showlegend=False,
        paper_bgcolor="#1a1a2e", plot_bgcolor="#1a1a2e", font_color="#ecf0f1",
        margin=dict(l=10, r=10, t=40, b=10),
    )
    return fig


# ── single-molecule results dashboard (simulation app) ────────────────────────
def render_single_dashboard(result: dict[str, Any], atoms: Atoms) -> None:
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots

    forces = result["forces"]
    syms = atoms.get_chemical_symbols()
    n = len(syms)
    labels = [f"{s}{i}" for i, s in enumerate(syms)]
    fx, fy, fz = forces[:, 0], forces[:, 1], forces[:, 2]
    fmag = np.linalg.norm(forces, axis=1)

    fig = make_subplots(
        rows=2, cols=2,
        subplot_titles=("Force Magnitude (eV/Å)", "Force Components (eV/Å)",
                        "Force Vectors 3-D", "Per-Atom Summary"),
        specs=[[{"type": "bar"}, {"type": "bar"}],
               [{"type": "scatter3d"}, {"type": "table"}]],
    )
    fig.add_trace(go.Bar(x=labels, y=fmag, name="|F|",
                         marker_color="steelblue", showlegend=False), row=1, col=1)
    for vals, name, color in [(fx, "Fx", "#e74c3c"), (fy, "Fy", "#2ecc71"), (fz, "Fz", "#f39c12")]:
        fig.add_trace(go.Bar(x=labels, y=vals, name=name, marker_color=color), row=1, col=2)

    pos = atoms.get_positions()
    scale = 0.5 / (fmag.max() + 1e-10)
    for i in range(n):
        x0, y0, z0 = pos[i]
        dx, dy, dz = forces[i] * scale
        fig.add_trace(go.Scatter3d(
            x=[x0, x0 + dx], y=[y0, y0 + dy], z=[z0, z0 + dz],
            mode="lines+markers", line=dict(width=4),
            marker=dict(size=[4, 8]), name=labels[i], showlegend=False,
        ), row=2, col=1)

    fig.add_trace(go.Table(
        header=dict(values=["Atom", "Fx", "Fy", "Fz", "|F|"],
                    fill_color="#2c3e50", font=dict(color="white")),
        cells=dict(values=[labels,
                           [f"{v:.4f}" for v in fx],
                           [f"{v:.4f}" for v in fy],
                           [f"{v:.4f}" for v in fz],
                           [f"{v:.4f}" for v in fmag]],
                   fill_color=[["#f0f4f8", "#ffffff"] * (n // 2 + 1)]),
    ), row=2, col=2)

    fig.update_layout(
        height=700, barmode="group", template="plotly_dark",
        paper_bgcolor="#1a1a2e", plot_bgcolor="#1a1a2e", font_color="#ecf0f1",
        title_text=(f"Potential Energy: {result['energy_eV']:.6f} eV"
                    + (f"  |  {result['note']}" if result.get("note") else "")),
    )
    st.plotly_chart(fig, use_container_width=True)


# ── reaction results dashboard (simulation app) ───────────────────────────────
def render_reaction_dashboard(scan: list[dict], atoms_a: Atoms, atoms_b: Atoms,
                              catalyst_atoms: Atoms | None,
                              decomp: dict, cat_name: str,
                              engine: str, basis: str, method: str) -> None:
    import plotly.graph_objects as go

    good = [r for r in scan if r["error"] is None]
    if not good:
        st.error("All scan points failed — check traceback above.")
        return

    dists = np.array([r["distance"] for r in good])
    e_rel = np.array([r["energy_rel"] for r in good])  # eV relative
    e_abs = np.array([r["energy_abs"] for r in good])

    # ── key metrics ────────────────────────────────────────────────────────────
    i_min = int(np.argmin(e_rel))
    i_max = int(np.argmax(e_rel))
    E_min = e_rel[i_min]
    E_max = e_rel[i_max]
    d_min = dists[i_min]
    has_barrier = i_max < i_min or (E_max > 0.005)
    E_act = E_max if has_barrier else 0.0
    delta_E = e_rel[-1] - e_rel[0]

    col1, col2, col3, col4 = st.columns(4)
    col1.metric("Interaction well (eV)", f"{E_min:.4f}", help="Min relative energy along scan")
    col2.metric("Activation barrier (eV)", f"{E_act:.4f}" if has_barrier else "None")
    col3.metric("ΔE (far→close, eV)", f"{delta_E:.4f}")
    demo_tag = " [demo]" if not decomp.get("demo") is False else ""
    col4.metric("Engine", f"{engine}/{method}/{basis}{demo_tag}")

    # ── energy profile ─────────────────────────────────────────────────────────
    pes_fig = go.Figure()
    pes_fig.add_trace(go.Scatter(
        x=dists, y=e_rel, mode="lines+markers",
        line=dict(color="#5b9bd5", width=2),
        marker=dict(size=8, color=e_rel, colorscale="RdYlGn_r",
                    showscale=True, colorbar=dict(title="eV")),
        name="E(d) – E(d_max)",
        hovertemplate="d=%{x:.2f} Å<br>ΔE=%{y:.4f} eV<extra></extra>",
    ))
    if has_barrier:
        pes_fig.add_vline(x=dists[i_max], line_dash="dash", line_color="#e74c3c",
                          annotation_text=f"‡ {E_act:.3f} eV")
    pes_fig.add_vline(x=d_min, line_dash="dot", line_color="#2ecc71",
                      annotation_text=f"min {E_min:.3f} eV")
    pes_fig.add_hline(y=0, line_dash="dot", line_color="#888888")

    pes_fig.update_layout(
        height=380, template="plotly_dark",
        paper_bgcolor="#1a1a2e", plot_bgcolor="#1a1a2e", font_color="#ecf0f1",
        xaxis_title="Separation A···B (Å)", yaxis_title="Relative Energy (eV)",
        title=f"Potential Energy Surface — {atoms_a.get_chemical_formula()} + "
              f"{atoms_b.get_chemical_formula()}"
              + (f" + {cat_name}" if cat_name != "None" else ""),
        showlegend=False,
    )
    st.plotly_chart(pes_fig, use_container_width=True)

    # ── energy decomposition bar chart ─────────────────────────────────────────
    decomp_labels = ["Reactant A", "Reactant B", "Catalyst", "Sum fragments", "Complex (closest)"]
    sum_frags = decomp["E_A"] + decomp["E_B"] + decomp["E_cat"]
    complex_e = e_abs[-1]  # closest approach
    interaction = complex_e - sum_frags if not decomp["demo"] else e_rel[-1]

    dec_fig = go.Figure()
    dec_fig.add_trace(go.Bar(
        x=decomp_labels,
        y=[decomp["E_A"], decomp["E_B"], decomp["E_cat"], sum_frags, complex_e],
        marker_color=["#5b9bd5", "#70ad47", "#ffc000", "#9b59b6", "#e74c3c"],
        text=[f"{v:.2f}" for v in [decomp["E_A"], decomp["E_B"], decomp["E_cat"],
                                   sum_frags, complex_e]],
        textposition="outside",
    ))
    dec_fig.update_layout(
        height=340, template="plotly_dark",
        paper_bgcolor="#1a1a2e", plot_bgcolor="#1a1a2e", font_color="#ecf0f1",
        yaxis_title="Energy (eV)",
        title=f"Energy Decomposition  |  Interaction energy: {interaction:.4f} eV",
        showlegend=False,
    )
    st.plotly_chart(dec_fig, use_container_width=True)

    # ── 3-D geometry viewer with scan slider ───────────────────────────────────
    st.subheader("Geometry along the reaction coordinate")
    n_good = len(good)
    pt = st.slider("Scan point", 0, n_good - 1, 0,
                   format=f"d=%(value)s of {n_good - 1}")
    selected = good[pt]
    nA = len(atoms_a)
    nB = len(atoms_b)
    nC = len(catalyst_atoms) if catalyst_atoms else 0
    viewer = mol_3d_fig(
        selected["atoms"],
        title=f"d(A···B) = {selected['distance']:.2f} Å  |  ΔE = {selected['energy_rel']:.4f} eV",
        height=380,
        group_sizes=(nA, nB, nC),
    )
    st.plotly_chart(viewer, use_container_width=True)
    st.caption("🔵 Molecule A  ·  🟢 Molecule B  ·  🟡 Catalyst")

    # ── data table ─────────────────────────────────────────────────────────────
    with st.expander("Full scan data"):
        import pandas as pd
        df = pd.DataFrame([
            {"d (Å)": r["distance"],
             "E_abs (eV)": r["energy_abs"],
             "ΔE (eV)": r["energy_rel"],
             "error": r["error"] or ""}
            for r in scan
        ])
        st.dataframe(df, use_container_width=True)
