"""Figure 6: the interface for a shallow source. Run after interface_pinn.py."""
import json
import os
import sys

sys.path.insert(0, 'src')

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from propagation import C_WATER

plt.rcParams.update({
    "figure.dpi": 130, "savefig.dpi": 130, "font.size": 9,
    "axes.grid": True, "grid.alpha": 0.25, "axes.spines.top": False,
    "axes.spines.right": False, "legend.frameon": False,
})


def main():
    d = json.load(open("results/interface_pinn.json"))
    xs = np.array([float(k) for k in d["vs_lambda_3d"]])
    ys = -10 * np.log10(np.array(list(d["vs_lambda_3d"].values())))
    plane = -10 * np.log10(d["plane_wave_fraction_3d"])
    pinn = d["pinn"].get("flat_3.0")
    have_field = os.path.exists("results/interface_field_flat_3.0.npz")
    n = 1 + (2 if (pinn and have_field) else 0)
    fig, axes = plt.subplots(1, n, figsize=(4.8 * n + 0.6, 3.9), squeeze=False)
    ax = axes[0]

    a = ax[0]
    a.semilogx(xs, ys, color="#2c6fb0", lw=1.6, label="point source, exact (wavenumber integral)")
    a.axhline(plane, color="#6a6a6a", ls="--", lw=1.2, label=f"plane-wave value, {plane:.1f} dB")
    for depth, f, mk, dy in ((1.0, 50.0, "o", 6), (3.0, 50.0, "s", -11), (1.0, 200.0, "^", 6), (3.0, 200.0, "D", -11)):
        x = depth / (C_WATER / f)
        y = np.interp(np.log(x), np.log(xs), ys)
        a.plot(x, y, mk, color="#b0442c", ms=6)
        a.annotate(f"{depth:.0f} m at {f:.0f} Hz: {y:.0f} dB", (x, y), fontsize=7,
                   xytext=(6, dy), textcoords="offset points")
    a.set_xlabel("source depth / wavelength in water")
    a.set_ylabel("loss into the air, dB")
    a.set_ylim(0, 40)
    a.legend(fontsize=7, loc="lower right")
    a.set_title("a shallow source is not 35 dB down\n(the evanescent near field crosses the surface)")

    if n == 3:
        z = np.load("results/interface_field_flat_3.0.npz")
        x, zz = z["x"], z["z"]
        air = zz > 0                          # the water halves agree to 0.5 %; show the air
        P = np.abs(z["pinn"][air]); Pr = np.abs(z["ref"][air])
        top = 20 * np.log10(P + 1e-12); bot = 20 * np.log10(Pr + 1e-12)
        half = len(x) // 2
        img = np.concatenate([top[:, :half], bot[:, half:]], 1)
        img = np.where(np.isfinite(img), img, np.nan)
        b = ax[1]
        vmax = np.nanpercentile(bot, 99.5); vmin = vmax - 40
        im = b.imshow(img, origin="lower", extent=[x[0], x[-1], 0, zz[-1]], vmin=vmin,
                      vmax=vmax, cmap="magma", aspect="auto")
        b.axvline(0, color="w", lw=0.6, ls=":")
        b.text(x[0] + 1, zz[-1] - 2.5, "PINN", color="w", fontsize=8)
        b.text(x[-1] - 6, zz[-1] - 2.5, "exact", color="w", fontsize=8)
        b.set_xlabel("x, m"); b.set_ylabel("height above the water, m")
        b.set_title("the air field of a line source 3 m down, |p| in dB:\nleft half the network, right half the exact field")
        b.grid(False)
        plt.colorbar(im, ax=b, fraction=0.04)

        c = ax[2]
        sx = np.array(pinn["surface_x"])
        c.plot(sx, 1e6 * np.array(pinn["surface_abs_p_ref"]), color="#2c6fb0", lw=1.6, label="exact")
        c.plot(sx, 1e6 * np.array(pinn["surface_abs_p_pinn"]), color="#1baf7a", lw=1.6, label="PINN")
        c.set_xlabel("x along the surface, m"); c.set_ylabel("|p| just above the surface, ×10⁻⁶")
        c.set_ylim(0, None)
        c.legend(fontsize=8)
        pd = -10 * np.log10(max(pinn["fraction"], 1e-12)); rd = -10 * np.log10(pinn["reference_2d"])
        c.set_title(f"where it goes wrong: the air pressure on the surface\n"
                    f"power into the air: exact {rd:.1f} dB, network {pd:.1f} dB")
    fig.tight_layout()
    fig.savefig("figures/06_shallow_source.png")
    plt.close(fig)
    print("ok")


if __name__ == "__main__":
    main()
