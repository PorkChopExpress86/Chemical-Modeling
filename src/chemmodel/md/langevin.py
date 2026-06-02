"""Ab-initio Langevin molecular-dynamics worker.

Runs off the Streamlit script thread: the worker NEVER touches ``st.*`` — it only
``put()``s plain data onto a queue the main thread drains. Cancellation is
cooperative via a ``threading.Event`` checked between steps (a Psi4 gradient
cannot be aborted mid-call, so stop takes effect after the current step).

Psi4 is a non-thread-safe process-global singleton, so MD must be single-flight
across the WHOLE process — :data:`MD_PROCESS_LOCK` enforces that exactly one
worker ever touches Psi4, even across browser refreshes / multiple sessions.
"""

from __future__ import annotations

import contextlib
import threading
import traceback

#: process-global guard — only one MD worker may touch Psi4 at a time
MD_PROCESS_LOCK = threading.Lock()


def md_worker(params, q, cancel):
    """Run Langevin MD, streaming step records to ``q``. No ``st.*`` calls here.

    ``params`` keys: ``atoms``, ``psi4_ok``, ``make_calc`` (→ ASE calculator),
    ``T``, ``dt`` (fs), ``fric`` (1/fs), ``n_steps``, ``group_sizes``.
    Queue messages: ``("step", {...})``, ``("done", {...})``, ``("stopped", n)``,
    ``("error", text)``.
    """
    if not MD_PROCESS_LOCK.acquire(blocking=False):
        q.put(("error", "Another molecular-dynamics run is already active in "
                         "this process. Wait for it to finish or stop it first."))
        return

    from ase import units
    from ase.md.langevin import Langevin
    from ase.md.velocitydistribution import (
        MaxwellBoltzmannDistribution, Stationary, ZeroRotation)

    atoms = params["atoms"]
    try:
        if params["psi4_ok"]:
            import psi4  # reset the singleton before we start

            with contextlib.suppress(Exception):
                psi4.core.clean()
                psi4.core.clean_options()
                psi4.core.clean_variables()

        atoms.calc = params["make_calc"]()
        MaxwellBoltzmannDistribution(atoms, temperature_K=params["T"])
        Stationary(atoms)
        ZeroRotation(atoms)
        dyn = Langevin(atoms, params["dt"] * units.fs,
                       temperature_K=params["T"],
                       friction=params["fric"] / units.fs)

        pe0 = None
        for step in range(params["n_steps"]):
            if cancel.is_set():
                q.put(("stopped", step))
                return
            dyn.run(1)
            epot = atoms.get_potential_energy()
            ekin = atoms.get_kinetic_energy()
            temp = atoms.get_temperature()
            if pe0 is None:
                pe0 = epot
            q.put(("step", {"i": step, "t": (step + 1) * params["dt"],
                            "pe": epot - pe0, "ke": ekin, "T": temp}))
        q.put(("done", {"positions": atoms.get_positions().copy(),
                        "symbols": list(atoms.get_chemical_symbols()),
                        "group_sizes": params["group_sizes"]}))
    except BaseException as exc:  # never raise on a worker thread
        q.put(("error", f"{type(exc).__name__}: {exc}\n{traceback.format_exc()}"))
    finally:
        with contextlib.suppress(Exception):
            calc = getattr(atoms, "calc", None)
            if calc is not None and hasattr(calc, "cleanup"):
                calc.cleanup()
        MD_PROCESS_LOCK.release()
