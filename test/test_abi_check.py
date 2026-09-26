"""Tests for the CoolProp ABI guard (PDSim/_abi_check.py)."""
import pytest

from PDSim._abi_check import check_coolprop_abi, CoolPropABIMismatch
from PDSim.misc import _coolprop_abi


def test_matching_abi_does_not_raise():
    """The binary we are running was built against the installed CoolProp, so
    the real check must find no mismatches."""
    assert check_coolprop_abi() == {}


def test_renumbered_enum_is_detected(monkeypatch):
    """If a baked-in enum value disagrees with the installed CoolProp, the
    guard reports it and raises with rebuild advice."""
    real = _coolprop_abi.build_abi_fingerprint()
    # Simulate a binary built against an older CoolProp where iDmass was 36.
    tampered = dict(real, iDmass=real['iDmass'] + 9999)
    monkeypatch.setattr(_coolprop_abi, 'build_abi_fingerprint', lambda: tampered)

    mismatches = check_coolprop_abi(raise_on_mismatch=False)
    assert 'iDmass' in mismatches
    assert mismatches['iDmass'][0] == tampered['iDmass']

    with pytest.raises(CoolPropABIMismatch) as exc:
        check_coolprop_abi()
    msg = str(exc.value)
    assert 'iDmass' in msg
    assert 'build_ext --inplace' in msg


def test_missing_constant_is_detected(monkeypatch):
    """A constant that no longer exists in the installed CoolProp is flagged."""
    real = _coolprop_abi.build_abi_fingerprint()
    tampered = dict(real, iNonexistentParameter=12345)
    monkeypatch.setattr(_coolprop_abi, 'build_abi_fingerprint', lambda: tampered)

    mismatches = check_coolprop_abi(raise_on_mismatch=False)
    assert mismatches.get('iNonexistentParameter') == (12345, None)
