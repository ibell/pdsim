"""
Tests that setDiscGeo assigns plain Python floats (not size-1 numpy arrays) to
the geoVals ``cdef double`` attributes.

numpy >= 1.25 deprecates converting a >0-dim array to a scalar; numpy will
eventually turn this into a hard TypeError. PDSim's setDiscGeo derives the arc
geometry from coords_norm(), which returns size-1 arrays, so without explicit
scalar coercion the cdef-double setters fire that deprecation.
"""
import warnings

import numpy as np
import pytest

from PDSim.scroll.scroll_geo import set_scroll_geo, geoVals
import PDSim.scroll.symm_scroll_geo as symm


@pytest.mark.parametrize("disc", [("2Arc", 0.0), ("2Arc", 0.0001), ("ArcLineArc", 0.001)])
def test_setdiscgeo_no_array_to_scalar_deprecation(disc):
    disc_type, disc_r2 = disc
    geo = geoVals()
    set_scroll_geo(
        0.0001048,  # Vdisp [m^3/rev]
        2.2,        # Vratio [-]
        0.004,      # thickness [m]
        0.005,      # orbiting radius [m]
        phi_i0=0.0,
        phi_os=0.3,
        phi_is=3.141,
        geo=geo,
    )
    with warnings.catch_warnings():
        warnings.simplefilter("error", DeprecationWarning)
        symm.setDiscGeo(geo, Type=disc_type, r2=disc_r2)

    # The geometry attributes must be plain scalars, not size-1 arrays.
    for attr in ("xa_arc1", "ya_arc1", "ra_arc1", "t1_arc1", "t2_arc1",
                 "xa_arc2", "ya_arc2", "ra_arc2", "t1_arc2", "t2_arc2"):
        value = getattr(geo, attr)
        assert np.isscalar(value) or isinstance(value, float)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
