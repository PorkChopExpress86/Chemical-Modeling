"""Command-line interface for the nitrogen-solubility model.

Presentation only: argument parsing and table printing. All physics lives in
:mod:`chemmodel.solubility.model`.
"""

from __future__ import annotations

import argparse
from typing import Optional, Sequence

from chemmodel.solubility.model import (
    SOLVENTS,
    MissingHenryParameterError,
    SolubilityResult,
    solvent_parameters,
    temperature_sweep,
)


def _format_float(value: float, digits: int = 4) -> str:
    return f"{value:.{digits}g}"


def print_sweep_table(rows: Sequence[SolubilityResult]) -> None:
    if not rows:
        print("No rows to display.")
        return

    solvent = rows[0].solvent
    solvent_name = rows[0].solvent_name
    print(f"\nN2 solubility in {solvent_name} ({solvent})")
    print("-" * 88)
    print(
        f"{'Temp (F)':>10} {'P (psig)':>10} {'P (psia)':>10} "
        f"{'Hcc (psia*ft3/mol)':>22} {'c_N2 (mol/ft3)':>18} {'x_N2':>12}"
    )
    print("-" * 88)
    for row in rows:
        print(
            f"{row.temperature_f:10.1f} "
            f"{row.pressure_psig:10.1f} "
            f"{row.pressure_psia:10.2f} "
            f"{row.henry_cc_psia_ft3_per_mol:22.2f} "
            f"{row.dissolved_mol_per_ft3:18.6f} "
            f"{_format_float(row.liquid_mole_fraction, 6):>12}"
        )


def print_parameter_notes(solvents: Sequence[str]) -> None:
    print("\nParameter notes")
    print("-" * 88)
    for solvent in solvents:
        params = solvent_parameters(solvent)
        print(f"{params.key}: {params.henry_reference}")
        print(f"{params.key}: {params.density_reference}")
        if params.notes:
            print(f"{params.key}: {params.notes}")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Nitrogen solubility in EO/PO at fixed psig.")
    parser.add_argument("--mode", choices=("sweep",), default="sweep")
    parser.add_argument("--solvent", choices=sorted(SOLVENTS), action="append")
    parser.add_argument("--pressure-psig", type=float, default=50.0)
    parser.add_argument("--t-min-f", type=float, default=68.0)
    parser.add_argument("--t-max-f", type=float, default=140.0)
    parser.add_argument("--t-step-f", type=float, default=12.0)
    parser.add_argument("--no-notes", action="store_true", help="Do not print parameter source notes.")
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    solvents = args.solvent or ["EO", "PO"]

    for solvent in solvents:
        try:
            rows = temperature_sweep(
                solvent=solvent,
                pressure_psig=args.pressure_psig,
                t_min_f=args.t_min_f,
                t_max_f=args.t_max_f,
                t_step_f=args.t_step_f,
            )
        except MissingHenryParameterError as exc:
            print(f"\nN2 solubility in {solvent_parameters(solvent).name} ({solvent})")
            print("-" * 88)
            print(f"Unavailable: {exc}")
            continue
        print_sweep_table(rows)

    if not args.no_notes:
        print_parameter_notes(solvents)
    return 0
