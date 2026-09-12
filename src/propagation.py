"""Acoustic detection budget for a two-channel hippopotamus warning system.

Two questions this module answers with physics rather than assertion:

1. How much of a submerged hippo's acoustic output actually reaches a
   listener standing in air?  This is fixed by the impedance contrast at the
   water/air interface and is the physical reason a person on the bank hears
   almost nothing while the animal is calling under water.

2. Given a source level and a background noise level, how far can each
   channel -- hydrophone in the water, microphone in the air -- detect the
   call?  Reported as a curve over source level, because no published source
   level for Hippopotamus amphibius exists; the curve is the honest form of
   the answer.

References for the constants are in docs/references.md.
"""

from __future__ import annotations

import numpy as np

# --- media -----------------------------------------------------------------
C_WATER = 1480.0        # m/s, fresh water at 20 C
C_AIR = 343.0           # m/s, 20 C
RHO_WATER = 998.0       # kg/m3
RHO_AIR = 1.204         # kg/m3
Z_WATER = RHO_WATER * C_WATER      # ~1.48e6 Pa*s/m
Z_AIR = RHO_AIR * C_AIR            # ~413 Pa*s/m


def refraction_cone_deg() -> float:
    """Half-angle of the cone in air into which all upward energy refracts.

    Snell's law with c_air < c_water bounds the refracted angle:
    sin(t2) = (c_air / c_water) * sin(t1) <= c_air / c_water.
    """
    return float(np.degrees(np.arcsin(C_AIR / C_WATER)))


def transmission_coefficient(theta1_rad):
    """Plane-wave power transmission, water -> air, at incidence theta1.

    Uses normal-incidence-corrected impedances Z/cos(theta).  Returns the
    fraction of incident intensity that crosses the interface.
    """
    t1 = np.asarray(theta1_rad, dtype=float)
    s2 = np.clip((C_AIR / C_WATER) * np.sin(t1), -1.0, 1.0)
    t2 = np.arcsin(s2)
    z1n = Z_WATER / np.maximum(np.cos(t1), 1e-12)
    z2n = Z_AIR / np.maximum(np.cos(t2), 1e-12)
    return 4.0 * z1n * z2n / (z1n + z2n) ** 2


def radiated_fraction_into_air(n: int = 200_001) -> float:
    """Fraction of a submerged omnidirectional source's total radiated power
    that crosses into the air, integrated over the upward hemisphere.
    """
    t1 = np.linspace(0.0, np.pi / 2, n)
    w = np.sin(t1)                       # solid-angle weight
    tc = transmission_coefficient(t1)
    upward = np.trapezoid(tc * w, t1)    # transmitted, upper hemisphere
    total = 2.0 * np.trapezoid(w, t1)    # full sphere, same weighting
    return float(upward / total)


def interface_loss_db() -> float:
    """Total insertion loss of the water/air interface, in dB."""
    return float(-10.0 * np.log10(radiated_fraction_into_air()))


# --- atmospheric absorption, ISO 9613-1 ------------------------------------
def air_absorption_db_per_m(f_hz, temp_c: float = 25.0, rel_hum: float = 60.0,
                            pressure_kpa: float = 101.325):
    """Pure-tone atmospheric attenuation coefficient, dB/m (ISO 9613-1:1993)."""
    f = np.asarray(f_hz, dtype=float)
    T = temp_c + 273.15
    T0 = 293.15
    T01 = 273.16
    pa = pressure_kpa
    pr = 101.325

    psat_over_pr = 10.0 ** (-6.8346 * (T01 / T) ** 1.261 + 4.6151)
    h = rel_hum * psat_over_pr / (pa / pr)          # molar water vapour, %

    frO = (pa / pr) * (24.0 + 4.04e4 * h * (0.02 + h) / (0.391 + h))
    frN = (pa / pr) * (T / T0) ** -0.5 * (
        9.0 + 280.0 * h * np.exp(-4.170 * ((T / T0) ** (-1.0 / 3.0) - 1.0))
    )

    term_class = 1.84e-11 * (pa / pr) ** -1 * (T / T0) ** 0.5
    term_O = 0.01275 * np.exp(-2239.1 / T) / (frO + f ** 2 / frO)
    term_N = 0.1068 * np.exp(-3352.0 / T) / (frN + f ** 2 / frN)
    return 8.686 * f ** 2 * (term_class + (T / T0) ** -2.5 * (term_O + term_N))


# --- transmission loss ------------------------------------------------------
def tl_air_db(range_m, f_hz, **kw):
    """Spherical spreading plus atmospheric absorption, re 1 m."""
    r = np.maximum(np.asarray(range_m, dtype=float), 1.0)
    return 20.0 * np.log10(r) + air_absorption_db_per_m(f_hz, **kw) * r


def tl_water_db(range_m, depth_m: float = 3.0, alpha_db_per_km: float = 0.0):
    """Shallow-water transmission loss re 1 m.

    Spherical spreading out to one water depth, cylindrical beyond it, plus an
    optional effective leakage term.  Fresh-water absorption below 5 kHz is
    < 1e-3 dB/km and is negligible; in a shallow body the real loss beyond a
    few hundred metres is bottom and surface interaction, which this term
    stands in for.  It is a free parameter, to be fixed from the recordings.
    """
    r = np.maximum(np.asarray(range_m, dtype=float), 1.0)
    r0 = max(depth_m, 1.0)
    near = 20.0 * np.log10(r)
    far = 20.0 * np.log10(r0) + 10.0 * np.log10(r / r0)
    return np.where(r <= r0, near, far) + alpha_db_per_km * r / 1000.0


def detection_range_m(source_level_db, noise_level_db, det_threshold_db,
                      medium: str, f_hz: float = 200.0, depth_m: float = 3.0,
                      alpha_db_per_km: float = 0.0, r_max: float = 20_000.0):
    """Largest range at which SL - TL - NL >= threshold.  Bisection on TL."""
    budget = float(source_level_db) - float(noise_level_db) - float(det_threshold_db)
    tl = (lambda r: tl_water_db(r, depth_m, alpha_db_per_km)) if medium == "water" \
        else (lambda r: tl_air_db(r, f_hz))
    if tl(1.0) > budget:
        return 0.0
    lo, hi = 1.0, r_max
    if tl(hi) <= budget:
        return float(hi)
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if tl(mid) <= budget:
            lo = mid
        else:
            hi = mid
    return float(lo)
