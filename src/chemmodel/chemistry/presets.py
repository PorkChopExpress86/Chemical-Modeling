"""Static preset libraries: molecules, catalysts, and reaction scenarios.

Pure data (XYZ strings + metadata) consumed by the simulation UI. Kept out of the
UI module so it can be reused and tested independently.
"""

from __future__ import annotations

import textwrap

PRESETS: dict[str, str] = {
    "Ethylene (C₂H₄)": textwrap.dedent("""\
        6
        Ethylene
        C   0.000000   0.000000   0.000000
        C   1.339000   0.000000   0.000000
        H  -0.540000   0.927000   0.000000
        H  -0.540000  -0.927000   0.000000
        H   1.879000   0.927000   0.000000
        H   1.879000  -0.927000   0.000000"""),
    "Water (H₂O)": textwrap.dedent("""\
        3
        Water
        O   0.000000   0.000000   0.119748
        H   0.000000   0.756950  -0.478993
        H   0.000000  -0.756950  -0.478993"""),
    "Methane (CH₄)": textwrap.dedent("""\
        5
        Methane
        C   0.000000   0.000000   0.000000
        H   0.629118   0.629118   0.629118
        H  -0.629118  -0.629118   0.629118
        H  -0.629118   0.629118  -0.629118
        H   0.629118  -0.629118  -0.629118"""),
    "CO₂": textwrap.dedent("""\
        3
        CO2
        C   0.000000   0.000000   0.000000
        O   1.162000   0.000000   0.000000
        O  -1.162000   0.000000   0.000000"""),
    "H₂": textwrap.dedent("""\
        2
        H2
        H   0.000000   0.000000   0.000000
        H   0.000000   0.000000   0.741000"""),
    "H (atom)": textwrap.dedent("""\
        1
        H-atom
        H   0.000000   0.000000   0.000000"""),
    "CO": textwrap.dedent("""\
        2
        CO
        C   0.000000   0.000000   0.000000
        O   0.000000   0.000000   1.128000"""),
    "NH₃": textwrap.dedent("""\
        4
        Ammonia
        N   0.000000   0.000000   0.000000
        H   0.000000   0.939692  -0.333333
        H   0.813493  -0.469846  -0.333333
        H  -0.813493  -0.469846  -0.333333"""),
    "HF": textwrap.dedent("""\
        2
        HF
        H   0.000000   0.000000   0.000000
        F   0.000000   0.000000   0.917000"""),
    "N₂": textwrap.dedent("""\
        2
        N2
        N   0.000000   0.000000   0.000000
        N   0.000000   0.000000   1.098000"""),
}

CATALYSTS: dict[str, dict] = {
    "None": {
        "xyz": None,
        "charge": 0,
        "label": "No catalyst — direct bimolecular",
        "color": "#555555",
    },
    "H⁺ (Brønsted acid)": {
        "xyz": "1\nProton\nH 0.0 0.0 0.0",
        "charge": +1,
        "label": "Proton — acid-catalysed reactions",
        "color": "#e74c3c",
    },
    "Pt₄ cluster": {
        "xyz": textwrap.dedent("""\
            4
            Pt4
            Pt  0.000000   0.000000   0.000000
            Pt  2.770000   0.000000   0.000000
            Pt  1.385000   2.399000   0.000000
            Pt  1.385000   0.799667   2.263000"""),
        "charge": 0,
        "label": "Pt₄ cluster — metallic heterogeneous catalyst",
        "color": "#95a5a6",
    },
    "Pd₄ cluster": {
        "xyz": textwrap.dedent("""\
            4
            Pd4
            Pd  0.000000   0.000000   0.000000
            Pd  2.750000   0.000000   0.000000
            Pd  1.375000   2.381000   0.000000
            Pd  1.375000   0.793667   2.245000"""),
        "charge": 0,
        "label": "Pd₄ cluster — hydrogenation catalyst",
        "color": "#bdc3c7",
    },
    "Al³⁺ (Lewis acid)": {
        "xyz": "1\nAlLewis\nAl 0.0 0.0 0.0",
        "charge": +3,
        "label": "Al³⁺ — Lewis acid catalysis",
        "color": "#f39c12",
    },
}

PRESET_REACTIONS: dict[str, dict] = {
    "H₂ + H → H abstraction": {
        "A": "H₂", "B": "H (atom)", "cat": "None",
        "start": 4.5, "end": 0.9, "n": 12,
    },
    "NH₃ + H⁺ → NH₄⁺ (protonation)": {
        "A": "NH₃", "B": "H (atom)", "cat": "H⁺ (Brønsted acid)",
        "start": 4.0, "end": 1.0, "n": 10,
    },
    "H₂ adsorption on Pt₄": {
        "A": "H₂", "B": "H₂", "cat": "Pt₄ cluster",
        "start": 5.0, "end": 1.5, "n": 10,
    },
    "CO + H₂O (water-gas shift)": {
        "A": "CO", "B": "Water (H₂O)", "cat": "None",
        "start": 4.5, "end": 1.5, "n": 10,
    },
    "H₂O dimer formation": {
        "A": "Water (H₂O)", "B": "Water (H₂O)", "cat": "None",
        "start": 5.0, "end": 2.0, "n": 12,
    },
    "CO₂ + H₂O (carbonic acid)": {
        "A": "CO₂", "B": "Water (H₂O)", "cat": "None",
        "start": 5.0, "end": 1.8, "n": 10,
    },
    "N₂ + H₂ on Pd₄ (Haber-Bosch step)": {
        "A": "N₂", "B": "H₂", "cat": "Pd₄ cluster",
        "start": 5.5, "end": 1.5, "n": 10,
    },
}
