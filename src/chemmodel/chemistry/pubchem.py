"""PubChem name → 3-D XYZ resolution. Pure (no Streamlit); the UI wraps this with
``st.cache_data`` so repeated look-ups and reruns don't re-hit the network.
"""

from __future__ import annotations


def fetch_pubchem_xyz(name: str) -> dict:
    """Resolve a chemical name to a 3-D XYZ string via PubChem.

    Returned dict: ``{ok, xyz, formula, cid, flat, error}``. Falls back to 2-D
    coordinates (``flat=True``) when no 3-D conformer is available.
    """
    import pubchempy as pcp

    name = name.strip()
    if not name:
        return {"ok": False, "error": "empty name"}
    try:
        comps = pcp.get_compounds(name, "name", record_type="3d")
        if not comps:
            # fall back to 2-D then let the caller know coords are flat
            comps = pcp.get_compounds(name, "name")
            if not comps:
                return {"ok": False, "error": f"'{name}' not found on PubChem"}
            c = comps[0]
            if not c.atoms:
                return {"ok": False, "error": f"'{name}' has no atom records"}
            flat = True
        else:
            flat = False

        c = comps[0]
        atoms = c.atoms
        lines = [str(len(atoms)), f"{name} (PubChem CID {c.cid})"]
        for a in atoms:
            z = a.z if a.z is not None else 0.0
            lines.append(f"{a.element}  {a.x:.6f}  {a.y:.6f}  {z:.6f}")
        return {
            "ok": True,
            "xyz": "\n".join(lines),
            "formula": c.molecular_formula or "?",
            "cid": c.cid,
            "flat": flat,
            "error": None,
        }
    except Exception as exc:  # network errors, timeouts, bad responses
        return {"ok": False, "error": f"{type(exc).__name__}: {exc}"}
