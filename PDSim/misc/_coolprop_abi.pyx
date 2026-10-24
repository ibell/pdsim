# cython: language_level=3
"""
Build-time fingerprint of the CoolProp enum values that PDSim bakes in.

PDSim's compiled extensions ``cimport CoolProp.constants_header`` and use
members such as ``constants.iDmass`` to read fluid properties.  Cython compiles
those enum members to *bare integer literals* in the generated C, so the
numbers are frozen at PDSim build time.

CoolProp renumbers its ``parameters`` / ``input_pairs`` enums between releases
(e.g. ``iDmass`` was 36 in 7.1.0 and 39 in 7.2.0; more move again in 8.x).
Nothing in the binary records which numbering it was built with -- the
``State`` C struct is unchanged, so Cython's own cross-module size/checksum
import check passes silently.  A PDSim binary built against one CoolProp then
asks a different runtime CoolProp for the wrong parameter, returning a garbage
density that only blows up much later as a cryptic
``ValueError: p is not a valid number`` deep inside ``solve()``.

The values captured here are whatever this binary was compiled against.  At
import, :func:`PDSim._abi_check.check_coolprop_abi` compares them against the
*installed* CoolProp (read from its pure-Python ``constants`` module) and fails
fast with rebuild instructions if they disagree.

Only the constants PDSim actually uses as values are tracked -- see the
``constants.*`` / ``constants_header.*`` usages in containers.pyx,
state_flooded.pyx and flow_models.pyx.
"""
cimport CoolProp.constants_header as ch


cpdef dict build_abi_fingerprint():
    """The CoolProp enum integer values this PDSim binary was compiled against."""
    return {
        # parameters enum (property output keys)
        'iT': <long>ch.iT,
        'iP': <long>ch.iP,
        'iQ': <long>ch.iQ,
        'iDmass': <long>ch.iDmass,
        'iHmass': <long>ch.iHmass,
        'iSmass': <long>ch.iSmass,
        'iUmass': <long>ch.iUmass,
        'iCpmass': <long>ch.iCpmass,
        'iCp0mass': <long>ch.iCp0mass,
        'iCvmass': <long>ch.iCvmass,
        'imolar_mass': <long>ch.imolar_mass,
        'iviscosity': <long>ch.iviscosity,
        'iconductivity': <long>ch.iconductivity,
        'ispeed_sound': <long>ch.ispeed_sound,
        'iT_critical': <long>ch.iT_critical,
        'iP_critical': <long>ch.iP_critical,
        # input_pairs enum (update() input specifiers)
        'PT_INPUTS': <long>ch.PT_INPUTS,
    }
