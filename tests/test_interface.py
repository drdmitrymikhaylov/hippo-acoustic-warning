"""Checks on the shallow-source interface calculation."""
import json
import os
import sys

sys.path.insert(0, 'src')

import numpy as np

import interface_pinn as I
from propagation import RHO_AIR, RHO_WATER, radiated_fraction_into_air


def test_deep_source_recovers_the_plane_wave_loss():
    """Many wavelengths down, the wavenumber integral must give the 35.5 dB of
    the plane-wave calculation on the project page."""
    deep = I.power_fraction_reference(3.0 * I.LAM_W, geometry="3d")
    assert abs(10 * np.log10(deep / radiated_fraction_into_air())) < 0.3


def test_transparency_depends_only_on_depth_in_wavelengths():
    a = I.power_fraction_reference(0.05 * 1480 / 50, geometry="3d", f_hz=50.0)
    b = I.power_fraction_reference(0.05 * 1480 / 200, geometry="3d", f_hz=200.0)
    assert abs(np.log10(a / b)) < 0.02


def test_shallow_source_loses_less_and_monotonically():
    ds = [0.25, 0.5, 1.0, 2.0, 5.0, 10.0]
    fr = [I.power_fraction_reference(d, geometry="3d") for d in ds]
    assert all(fr[i] > fr[i + 1] for i in range(len(fr) - 1))
    assert -10 * np.log10(fr[2]) < 15.0          # 1 m at 50 Hz: well under 15 dB
    assert -10 * np.log10(fr[-1]) > 33.0         # 10 m: back near the plane-wave value


def test_reference_field_reduces_to_the_free_field_without_an_interface():
    from scipy.special import hankel1
    orig = I.reflection_transmission
    I.reflection_transmission = lambda kx: (orig(kx)[0] * 0, orig(kx)[1] * 0 + 1, orig(kx)[2], orig(kx)[2])
    try:
        p = I.field_reference(np.array([3.0]), np.array([-25.0]), 20.0)[0]
    finally:
        I.reflection_transmission = orig
    p0 = 0.25j * hankel1(0, I.K_W * np.hypot(3.0, 5.0))
    assert abs(p - p0) / abs(p0) < 0.01


def test_reference_field_satisfies_both_interface_conditions():
    x = np.array([0.0, 3.0, 10.0]); h = 1e-3
    pw = I.field_reference(x, np.full(3, -1e-6), 1.0); pa = I.field_reference(x, np.full(3, 1e-6), 1.0)
    assert np.all(np.abs(pw - pa) < 1e-2 * np.abs(pw))
    vw = (pw - I.field_reference(x, np.full(3, -h), 1.0)) / h / RHO_WATER
    va = (I.field_reference(x, np.full(3, h), 1.0) - pa) / h / RHO_AIR
    assert np.all(np.abs(vw - va) < 2e-2 * np.abs(vw))


def test_pinn_result_is_what_the_page_says():
    """The network satisfies the Neumann data by construction and gets the
    water right; the page says it does NOT reproduce the air power.  If a
    future run does, this check should fail and the page be rewritten."""
    p = "results/interface_pinn.json"
    if not os.path.exists(p):
        return
    r = json.load(open(p))["pinn"]["flat_3.0"]
    assert r["rel_field_error_water"] < 0.05
    gap_db = 10 * np.log10(r["reference_2d"] / max(r["fraction"], 1e-12))
    assert gap_db > 3.0, gap_db

if __name__ == "__main__":
    fns = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    failed = 0
    for fn in fns:
        try:
            fn(); print(f"PASS  {fn.__name__}")
        except AssertionError as e:
            failed += 1; print(f"FAIL  {fn.__name__}: {e}")
    print(f"\n{len(fns) - failed}/{len(fns)} passed")
    sys.exit(1 if failed else 0)
