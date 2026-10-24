"""
Import-time guard that PDSim's compiled extensions match the installed CoolProp.

CoolProp renumbers its property/input enums between releases.  Because PDSim's
Cython code compiles those enum members to integer literals, a PDSim binary
built against one CoolProp silently reads the wrong properties from another --
see ``PDSim/misc/_coolprop_abi.pyx`` for the full explanation.

This module compares the enum values frozen into the binary at build time
against the values reported by the *installed* CoolProp, and raises a clear,
actionable error at ``import PDSim`` instead of letting the mismatch surface as
a cryptic ``ValueError: p is not a valid number`` deep inside ``solve()``.
"""
from __future__ import absolute_import


class CoolPropABIMismatch(ImportError):
    """PDSim's compiled CoolProp enum values disagree with the installed
    CoolProp -- PDSim must be rebuilt against the installed CoolProp."""
    pass


def _format_message(mismatches, version):
    lines = [
        "PDSim was compiled against a different CoolProp than the one "
        "installed (CoolProp {0}).".format(version),
        "",
        "CoolProp renumbers its property/input enums between releases, and "
        "PDSim bakes those numbers into its compiled extensions, so a PDSim "
        "binary built against one CoolProp reads the wrong properties from "
        "another -- typically surfacing much later as a cryptic "
        "'ValueError: p is not a valid number' inside solve().",
        "",
        "Mismatched enum values (name: built -> installed):",
    ]
    for name in sorted(mismatches):
        built, installed = mismatches[name]
        shown = 'MISSING' if installed is None else installed
        lines.append("    {0}: {1} -> {2}".format(name, built, shown))
    lines += [
        "",
        "Resolution -- rebuild PDSim against the installed CoolProp:",
        "  * from a source checkout (rebuild the in-place extensions):",
        "        python setup.py build_ext --inplace --force",
        "  * installed with pip (force a fresh compile, no cached wheel):",
        "        pip install --force-reinstall --no-binary PDSim PDSim",
        "  * with uv (a plain `uv sync` reuses uv's cached build -- bust it):",
        "        uv cache clean pdsim && uv sync",
    ]
    return "\n".join(lines)


def check_coolprop_abi(raise_on_mismatch=True):
    """Compare the CoolProp enum values baked into PDSim against the installed
    CoolProp.

    Returns a dict ``{name: (built_value, installed_value)}`` of any mismatches
    (empty if all agree).  When ``raise_on_mismatch`` is true (the default), a
    non-empty result raises :class:`CoolPropABIMismatch` with rebuild
    instructions.

    Returns ``{}`` without raising if either the compiled fingerprint module or
    CoolProp cannot be imported -- there is then nothing to compare and the
    real import error will surface on its own.
    """
    try:
        from PDSim.misc._coolprop_abi import build_abi_fingerprint
    except Exception:
        return {}
    try:
        import CoolProp
        from CoolProp import constants
    except Exception:
        return {}

    built = build_abi_fingerprint()
    mismatches = {}
    for name, built_value in built.items():
        installed_value = getattr(constants, name, None)
        if installed_value is None:
            mismatches[name] = (built_value, None)
        elif int(installed_value) != int(built_value):
            mismatches[name] = (built_value, int(installed_value))

    if mismatches and raise_on_mismatch:
        version = getattr(CoolProp, '__version__', 'unknown')
        raise CoolPropABIMismatch(_format_message(mismatches, version))

    return mismatches
