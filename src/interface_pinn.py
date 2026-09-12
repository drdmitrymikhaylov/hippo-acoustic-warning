"""Sound crossing the water surface when the source is *shallow*.

The 35.5 dB interface loss on the project page is a plane-wave number: it
holds for a source many wavelengths below the surface, where the field at
the interface is a sum of propagating plane waves and Snell's law and the
impedance ratio decide everything.  A hippopotamus does not call from many
wavelengths down.  At 50 Hz the wavelength in water is 30 m; an animal at
1-3 m is a fraction of a wavelength from the surface, and in that regime the
field at the interface is dominated by *evanescent* components, which the
plane-wave calculation does not contain.  The literature calls the result
"anomalous transparency" (Godin, 2006-2008): a shallow source loses far less
than 35 dB into the air.

This file computes that effect two ways and checks them against each other.

  reference   the exact two-medium solution for a line source under a flat
              interface, by wavenumber integration (Sommerfeld integral),
              in numpy -- no network.
  PINN        two networks, one per medium, solving the Helmholtz equation
              in each and joined by the interface conditions,

                  p_w = p_a                  (pressure is continuous)
                  (1/rho_w) dp_w/dz = (1/rho_a) dp_a/dz   (normal velocity is continuous)

              with a hard-coded free-field source term: the network
              represents only the *scattered* field, so the singularity is
              never learned.  A network is not needed for the flat interface
              -- that is what the reference is for -- but it is needed for a
              sloping bank, which has no closed form, and the flat case is
              how the network earns the right to be believed there.

Both report the same quantity: the fraction of the source's power that ends
up in the air, as a function of source depth in wavelengths, next to the
plane-wave value.
"""

import json
import sys
import time

sys.path.insert(0, 'src')

import numpy as np

from propagation import C_AIR, C_WATER, RHO_AIR, RHO_WATER, radiated_fraction_into_air

F_HZ = 50.0
K_W = 2 * np.pi * F_HZ / C_WATER
K_A = 2 * np.pi * F_HZ / C_AIR
LAM_W = C_WATER / F_HZ


# ----------------------------------------------------------------- reference
def kz(kx, k):
    """Vertical wavenumber with the radiation branch (Im >= 0)."""
    return np.sqrt(k ** 2 - kx ** 2 + 0j) * np.where(k ** 2 - kx ** 2 >= 0, 1.0, 1.0)


def reflection_transmission(kx):
    """Plane-wave coefficients for a wave incident from the water (z < 0) on
    the interface z = 0 with air above, for horizontal wavenumber kx.  Pressure
    convention: p_w = e^{i kzw z} + R e^{-i kzw z}, p_a = T e^{i kza z}."""
    kzw = np.sqrt(K_W ** 2 - kx ** 2 + 0j)
    kza = np.sqrt(K_A ** 2 - kx ** 2 + 0j)
    kzw = np.where(kzw.imag < 0, -kzw, kzw)
    kza = np.where(kza.imag < 0, -kza, kza)
    a = kzw / RHO_WATER
    b = kza / RHO_AIR
    R = (a - b) / (a + b)
    T = 2 * a / (a + b)
    return R, T, kzw, kza


def power_fraction_reference(depth_m, n=200_001, kx_max=None, geometry="2d", f_hz=None):
    """Fraction of the source's radiated power that crosses into the air, for
    a source `depth_m` below the surface.  Wavenumber integral: transmitted
    intensity through the plane z = 0+ integrated over horizontal wavenumber
    (Parseval), divided by the total power radiated (into water + air).
    geometry "2d" is a line source (what the PINN solves); "3d" is a point
    source, whose deep limit is the 35.5 dB plane-wave number of the page."""
    global K_W, K_A
    if f_hz is not None:
        kw0, ka0 = K_W, K_A
        K_W, K_A = 2 * np.pi * f_hz / C_WATER, 2 * np.pi * f_hz / C_AIR
    kx_max = kx_max or 12.0 * K_A
    kx = np.linspace(0.0, kx_max, n)
    R, T, kzw, kza = reflection_transmission(kx)
    weight = kx if geometry == "3d" else np.ones_like(kx)
    zs = -depth_m
    # incident plane-wave spectrum of a line source at depth: amplitude at z=0
    # is e^{i kzw d}; use pressure spectrum P(kx) = e^{i kzw d} / (i kzw) (2-D
    # Weyl / Sommerfeld representation, common factor dropped)
    inc = np.exp(1j * kzw * (-zs)) / (1j * kzw)
    # transmitted: T * inc, propagating in air only where kza is real
    trans = T * inc
    Ia = np.where(kza.imag == 0, np.abs(trans) ** 2 * kza.real / RHO_AIR, 0.0)
    # reflected + incident field in water: outgoing power through z -> -inf is
    # |1 + R e^{...}|^2 evaluated as downgoing spectrum: down-going amplitude
    # at large negative z is inc * (1 + R) ... careful: the down-going wave is
    # the direct wave below the source (amplitude 1/(i kzw)) plus the
    # reflected wave R e^{i kzw d} / (i kzw) e^{i kzw d}; both go down.
    down = (1.0 + R * np.exp(2j * kzw * depth_m)) / (1j * kzw)
    Iw = np.where(kzw.imag == 0, np.abs(down) ** 2 * kzw.real / RHO_WATER, 0.0)
    Pa = np.trapezoid(Ia * weight, kx)
    Pw = np.trapezoid(Iw * weight, kx)
    if f_hz is not None:
        K_W, K_A = kw0, ka0
    return float(Pa / (Pa + Pw))


def field_reference(x, z, depth_m, n=60_001, kx_max=None, eps=2e-3):
    """Complex pressure field of the line source (water below z=0, air above),
    for plotting and for scoring the PINN.  Normalised so the free-field
    source is p0 = (i/4) H0(k r) -- i.e. spectrum e^{i kzw |z - zs|}/(i kzw)
    times 1/(2 pi) integrated over kx, with the interface terms added."""
    global K_W, K_A
    kx_max = kx_max or 12.0 * K_A
    kx = np.linspace(-kx_max, kx_max, n)
    # a small loss keeps the 1/kz branch-point singularity off the grid
    kw0, ka0 = K_W, K_A
    K_W, K_A = kw0 * (1 + 1j * eps), ka0 * (1 + 1j * eps)
    R, T, kzw, kza = reflection_transmission(kx)
    K_W, K_A = kw0, ka0
    zs = -depth_m
    xf = np.asarray(x, float).ravel(); zf = np.asarray(z, float).ravel()
    out = np.empty(xf.size, complex)
    A = 1j / kzw                                # p0 = (i/4) H0 = (i/4pi) int e^{i kz|z|} e^{i kx x} / kz dkx
    dk = (kx[1] - kx[0]) / (4 * np.pi)
    for i0 in range(0, xf.size, 256):           # chunked: the full outer product does not fit in memory
        X = xf[i0:i0 + 256].reshape(-1, 1); Z = zf[i0:i0 + 256].reshape(-1, 1)
        phase = np.exp(1j * kx.reshape(1, -1) * X)
        below = (Z < 0).ravel()
        direct = np.exp(1j * kzw.reshape(1, -1) * np.abs(Z - zs))
        refl = R.reshape(1, -1) * np.exp(1j * kzw.reshape(1, -1) * (depth_m - Z))
        pw = (A * (direct + refl) * phase).sum(1) * dk
        pa = (A * T.reshape(1, -1) * np.exp(1j * kzw.reshape(1, -1) * depth_m)
              * np.exp(1j * kza.reshape(1, -1) * Z) * phase).sum(1) * dk
        out[i0:i0 + 256] = np.where(below, pw, pa)
    return out


# ---------------------------------------------------------------------- PINN
SW = 0.25                                   # scale of the water field (|p0| near the source)
SA = SW * 2 * K_W * RHO_AIR / (K_A * RHO_WATER)   # scale of the air field (plane-wave T)
# The water side of a water/air interface is a pressure-release surface to
# O(rho_a/rho_w): the field in the water is the source plus its negative
# image, and only a correction of order SA is left for a network to find.
# That correction is what the water network represents; the air network
# represents the whole (small) air field.


def free_field(x, z, zs, xs=0.0):
    """Free-field line source p0 = (i/4) H0(k_w r) and its gradient, as
    numpy complex arrays (the source is fixed; nothing here is learned)."""
    from scipy.special import hankel1
    r = np.sqrt((x - xs) ** 2 + (z - zs) ** 2) + 1e-9
    h0, h1 = hankel1(0, K_W * r), hankel1(1, K_W * r)
    p = 0.25j * h0
    g = -0.25j * K_W * h1
    return p, g * (x - xs) / r, g * (z - zs) / r


def source_and_image(x, z, zs, slope_deg=0.0):
    """Source minus its mirror image across the surface z = x tan(slope):
    the pressure-release solution, exact for an infinitely light air."""
    th = np.radians(slope_deg)
    xi, zi = zs * np.sin(2 * th), -zs * np.cos(2 * th)      # image of (0, zs)
    p, gx, gz = free_field(x, z, zs)
    pi, gxi, gzi = free_field(x, z, zi, xs=xi)
    return p - pi, gx - gxi, gz - gzi


def train_pinn(depth_m, slope_deg=0.0, steps=4000, seed=0, log=print, R_w=None, R_a=None,
               n_col=1536, lbfgs=150):
    """Two networks (water, air) for the *scattered* field, joined at the
    interface z = x tan(slope) by pressure and normal-velocity continuity.
    Each medium is a half-disc centred on the surface point above the source,
    closed by a second-order (Bayliss-Turkel) radiation condition on its arc,
    dp/dr = (i k - 1/(2 r)) p, which is what an outgoing cylindrical wave
    satisfies -- a first-order condition on a box reflected enough grazing
    energy to spoil the air field.  Returns the fraction of power crossing
    into the air."""
    import torch
    torch.set_default_dtype(torch.float64)
    torch.manual_seed(seed)
    R_w = R_w or 0.6 * LAM_W                 # 18 m at 50 Hz
    R_a = R_a or 3.0 * (C_AIR / F_HZ)        # 3 air wavelengths, 21 m
    th = np.radians(slope_deg); tan = float(np.tan(th))
    zs = -depth_m
    nvec_s = np.array([-np.sin(th), np.cos(th)])              # normal of the surface, into the air

    class Field(torch.nn.Module):
        def __init__(self, k, scale, width=96, n_feat=48, seed=0):
            super().__init__()
            g = torch.Generator().manual_seed(seed)
            # Fourier features from half the wavenumber to 2.5x it: enough for
            # the near field of a source a tenth of a wavelength down; wider
            # ranges made the second derivatives in the PDE loss noisy
            sc = torch.exp(torch.rand(1, n_feat, generator=g) * np.log(5.0) + np.log(0.5))
            self.register_buffer("B", torch.randn(2, n_feat, generator=g) * k * sc)
            self.scale = scale
            self.net = torch.nn.Sequential(torch.nn.Linear(2 * n_feat, width), torch.nn.Tanh(),
                                           torch.nn.Linear(width, width), torch.nn.Tanh(),
                                           torch.nn.Linear(width, width), torch.nn.Tanh(),
                                           torch.nn.Linear(width, 2))

        def forward(self, xz):
            f = xz @ self.B
            return self.scale * self.net(torch.cat([torch.sin(f), torch.cos(f)], 1))

    # --- the air network carries the Neumann data as a HARD constraint.
    # A soft interface loss let the network match dp/dn with a thin boundary
    # layer of tiny amplitude and violate the Helmholtz equation next to the
    # surface (residual 10-90 % of k^2 p at z = 1 m while the domain-averaged
    # loss looked small).  So the air field is written as
    #     p_a = g(t) s e^{-s/l}  +  N(t, q(s)),   q = l (sqrt(1 + (s/l)^2) - 1),
    # with s the normal distance from the surface, t the coordinate along it,
    # g(t) = (rho_a/rho_w) * dp_w/dn the analytic Neumann data (a Chebyshev
    # fit, so it is differentiable in torch), l = 1/k_a.  dp_a/dn = g at s = 0
    # by construction: q is even in s, so N's normal derivative vanishes on
    # the surface while its value there is free.
    ell = 1.0 / K_A
    e_t = np.array([np.cos(th), np.sin(th)])
    tt = np.linspace(-R_a, R_a, 4001)
    _, gx0, gz0 = source_and_image(tt * e_t[0], tt * e_t[1], zs, slope_deg)
    g_t = (RHO_AIR / RHO_WATER) * (gx0 * nvec_s[0] + gz0 * nvec_s[1])
    deg = 120
    cre = np.polynomial.chebyshev.chebfit(tt / R_a, g_t.real, deg)
    cim = np.polynomial.chebyshev.chebfit(tt / R_a, g_t.imag, deg)
    fit_err = np.max(np.abs(np.polynomial.chebyshev.chebval(tt / R_a, cre)
                            + 1j * np.polynomial.chebyshev.chebval(tt / R_a, cim) - g_t)) / np.max(np.abs(g_t))
    log(f"    Neumann data: Chebyshev fit of degree {deg}, max error {fit_err:.1e} of the peak")
    C_re, C_im = torch.tensor(cre), torch.tensor(cim)

    def chebval(u, c):
        """Clenshaw recurrence, differentiable."""
        b1 = torch.zeros_like(u); b2 = torch.zeros_like(u)
        for k in range(len(c) - 1, 0, -1):
            b1, b2 = 2 * u * b1 - b2 + c[k], b1
        return u * b1 - b2 + c[0]

    class AirField(torch.nn.Module):
        def __init__(self):
            super().__init__()
            self.N = Field(K_A, SA, seed=seed + 1)

        def forward(self, xz):
            sdist = (xz[:, 0] * nvec_s[0] + xz[:, 1] * nvec_s[1]).reshape(-1, 1)
            t = (xz[:, 0] * e_t[0] + xz[:, 1] * e_t[1]).reshape(-1, 1)
            u = torch.clamp(t / R_a, -1.0, 1.0)
            g = torch.cat([chebval(u, C_re), chebval(u, C_im)], 1)
            # the free part is a function of (t, s^2): its normal derivative
            # vanishes at the surface, so dp/dn = g there, while its value at
            # the surface is unconstrained (an earlier version multiplied the
            # network by s^2 and thereby forced p = 0 on the surface as well --
            # a Dirichlet condition nobody asked for)
            # q(s) = l (sqrt(1 + (s/l)^2) - 1): even in s, ~ s^2/2l near the
            # surface (zero normal derivative there), ~ s far from it, so the
            # Fourier features see ordinary distances
            q = ell * (torch.sqrt(1 + (sdist / ell) ** 2) - 1)
            return g * sdist * torch.exp(-sdist / ell) + self.N(torch.cat([t, q], 1))

    fw, fa = Field(K_W, SA, seed=seed), AirField()
    params = list(fw.parameters()) + list(fa.parameters())

    def c2(a):
        return torch.tensor(np.stack([a.real, a.imag], 1))

    def grads(f, xz):
        xz = xz.clone().requires_grad_(True)
        u = f(xz)
        g = [torch.autograd.grad(u[:, c].sum(), xz, create_graph=True)[0] for c in range(2)]
        gx = torch.stack([g[0][:, 0], g[1][:, 0]], 1)
        gz = torch.stack([g[0][:, 1], g[1][:, 1]], 1)
        return xz, u, gx, gz

    def laplacian(f, xz, k):
        xz, u, gx, gz = grads(f, xz)
        res = []
        for c in range(2):
            gxx = torch.autograd.grad(gx[:, c].sum(), xz, create_graph=True)[0][:, 0]
            gzz = torch.autograd.grad(gz[:, c].sum(), xz, create_graph=True)[0][:, 1]
            res.append(gxx + gzz + k ** 2 * u[:, c])
        return torch.stack(res, 1)

    def total_water(xz):
        """Total field and gradient in the water: source + image + network."""
        _, u, gx, gz = grads(fw, xz)
        p0, g0x, g0z = source_and_image(xz[:, 0].numpy(), xz[:, 1].numpy(), zs, slope_deg)
        return u + c2(p0), gx + c2(g0x), gz + c2(g0z)

    def bt1(u, gx, gz, xz, k):
        """Bayliss-Turkel radiation residual on an arc centred at the origin:
        dp/dr - (i k - 1/(2r)) p, for p = pr + i pi."""
        r = torch.sqrt(xz[:, 0] ** 2 + xz[:, 1] ** 2).reshape(-1, 1)
        nx, nz = xz[:, 0:1] / r, xz[:, 1:2] / r
        dn = gx * nx + gz * nz
        ikp = torch.stack([-k * u[:, 1], k * u[:, 0]], 1)
        return dn - ikp + u / (2 * r)

    def half_disc(n, medium, near_frac=0.5):
        """Random points in the half-disc of the medium; half of them within
        two source depths of the surface where the field is sharpest."""
        R = R_w if medium == "w" else R_a
        sgn = -1.0 if medium == "w" else 1.0
        r = R * torch.sqrt(torch.rand(n))
        phi = torch.rand(n) * np.pi
        # local frame: e_t along the surface, e_n into the air
        t = r * torch.cos(phi); nn_ = sgn * r * torch.sin(phi)
        near = torch.rand(n) < near_frac
        nn_ = torch.where(near, sgn * torch.rand(n) * 2 * depth_m, nn_)
        t = torch.where(near, (torch.rand(n) * 2 - 1) * R, t)
        x = t * np.cos(th) - nn_ * np.sin(th)
        z = t * np.sin(th) + nn_ * np.cos(th)
        return torch.stack([x, z], 1)

    def arc(n, medium):
        R = R_w if medium == "w" else R_a
        sgn = -1.0 if medium == "w" else 1.0
        phi = torch.rand(n) * np.pi
        t = R * torch.cos(phi); nn_ = sgn * R * torch.sin(phi)
        return torch.stack([t * np.cos(th) - nn_ * np.sin(th), t * np.sin(th) + nn_ * np.cos(th)], 1)

    def losses(n=n_col):
        l_pde = torch.mean(laplacian(fw, half_disc(n, "w"), K_W) ** 2) / (K_W ** 2 * SA) ** 2 \
            + torch.mean(laplacian(fa, half_disc(n, "a"), K_A) ** 2) / (K_A ** 2 * SA) ** 2
        # interface: half the points within 3 depths of the source, where the data lives
        t1 = (torch.rand(n // 2) * 2 - 1) * min(R_w, R_a)
        t2 = (torch.rand(n - n // 2) * 2 - 1) * 3 * depth_m
        t = torch.cat([t1, t2])
        xi = torch.stack([t * np.cos(th), t * np.sin(th)], 1)
        pw, gwx, gwz = total_water(xi)
        _, pa, gax, gaz = grads(fa, xi)
        dnw = gwx * nvec_s[0] + gwz * nvec_s[1]
        dna = gax * nvec_s[0] + gaz * nvec_s[1]
        l_p = torch.mean((pw - pa) ** 2) / SA ** 2
        l_v = torch.mean((dna - (RHO_AIR / RHO_WATER) * dnw) ** 2) / (K_A * SA) ** 2
        # radiation on the arcs: the source and image are exactly outgoing, so
        # in the water only the network's correction is constrained
        xw = arc(n // 2, "w"); xw_, u, gx, gz = grads(fw, xw)
        l_r = torch.mean(bt1(u, gx, gz, xw_, K_W) ** 2) / (K_W * SA) ** 2
        xa = arc(n // 2, "a"); xa_, u, gx, gz = grads(fa, xa)
        l_r = l_r + torch.mean(bt1(u, gx, gz, xa_, K_A) ** 2) / (K_A * SA) ** 2
        return l_pde, l_p, l_v, l_r

    def total_loss():
        a, b, c, d = losses()
        return a + 10.0 * b + d, (a, b, c, d)          # c (velocity) is exact by construction

    opt = torch.optim.Adam(params, lr=1e-3)
    sched = torch.optim.lr_scheduler.CosineAnnealingLR(opt, steps)
    t0 = time.time()
    for it in range(steps):
        opt.zero_grad()
        loss, (a, b, c, d) = total_loss()
        loss.backward()
        opt.step(); sched.step()
        if it % 500 == 0 or it == steps - 1:
            log(f"    step {it:5d}  pde {float(a):.2e}  p-cont {float(b):.2e}  "
                f"v-cont {float(c):.2e}  radiation {float(d):.2e}  [{time.time()-t0:.0f} s]")
    if lbfgs:
        torch.manual_seed(seed + 7)
        fixed = {}
        opt2 = torch.optim.LBFGS(params, max_iter=lbfgs, history_size=50, line_search_fn="strong_wolfe")

        def closure():
            opt2.zero_grad()
            torch.manual_seed(seed + 7)          # same collocation set every evaluation
            loss, parts = total_loss()
            loss.backward()
            fixed["parts"] = parts
            return loss
        opt2.step(closure)
        a, b, c, d = fixed["parts"]
        log(f"    L-BFGS  pde {float(a):.2e}  p-cont {float(b):.2e}  v-cont {float(c):.2e}  "
            f"radiation {float(d):.2e}  [{time.time()-t0:.0f} s]")

    # power through the arcs: I_n = Im(conj(p) dp/dn) / (2 rho omega); the
    # common 1/(2 omega) cancels in the fraction
    def flux(u, gx, gz, xz, rho):
        r = torch.sqrt(xz[:, 0] ** 2 + xz[:, 1] ** 2).reshape(-1, 1)
        dn = gx * xz[:, 0:1] / r + gz * xz[:, 1:2] / r
        return ((u[:, 0] * dn[:, 1] - u[:, 1] * dn[:, 0]) / rho).detach()

    nl = 4001
    phi = torch.linspace(0, np.pi, nl); dphi = float(phi[1] - phi[0])
    def arc_pts(R, sgn):
        t = R * torch.cos(phi); nn_ = sgn * R * torch.sin(phi)
        return torch.stack([t * np.cos(th) - nn_ * np.sin(th), t * np.sin(th) + nn_ * np.cos(th)], 1)
    xw = arc_pts(R_w, -1.0); u, gx, gz = total_water(xw)
    Pw = float(flux(u, gx, gz, xw, RHO_WATER).sum() * dphi * R_w)
    xa = arc_pts(R_a, 1.0); _, u, gx, gz = grads(fa, xa)
    Pa = float(flux(u, gx, gz, xa, RHO_AIR).sum() * dphi * R_a)
    # air power straight through the surface, for a cross-check
    t = torch.linspace(-R_a, R_a, nl); dt = float(t[1] - t[0])
    top = torch.stack([t * np.cos(th) - 1e-3 * np.sin(th), t * np.sin(th) + 1e-3 * np.cos(th)], 1)
    _, u, gx, gz = grads(fa, top)
    dn = gx * nvec_s[0] + gz * nvec_s[1]
    Pa_surf = float(((u[:, 0] * dn[:, 1] - u[:, 1] * dn[:, 0]) / RHO_AIR).detach().sum() * dt)

    R = max(R_w, R_a)
    gx_ = np.linspace(-R, R, 161); gz_ = np.linspace(-R, R, 161)
    GX, GZ = np.meshgrid(gx_, gz_)
    pts = torch.tensor(np.stack([GX.ravel(), GZ.ravel()], 1))
    with torch.no_grad():
        ua = fa(pts).numpy(); uw = fw(pts).numpy()
    p0, _, _ = source_and_image(GX.ravel(), GZ.ravel(), zs, slope_deg)
    pw_tot = uw + np.stack([p0.real, p0.imag], 1)
    below = GZ.ravel() < GX.ravel() * tan
    field = np.where(below[:, None], pw_tot, ua)
    inside = (GX.ravel() ** 2 + GZ.ravel() ** 2) <= np.where(below, R_w, R_a) ** 2
    field = np.where(inside[:, None], field, np.nan)
    return {"P_air": Pa, "P_air_through_surface": Pa_surf, "P_water": Pw, "_nets": (fw, fa, grads, laplacian, total_water),
            "fraction": Pa / (Pa + Pw), "R_w": R_w, "R_a": R_a, "depth_m": depth_m, "slope_deg": slope_deg,
            "grid_x": gx_, "grid_z": gz_, "field": field[:, 0].reshape(GZ.shape) + 1j * field[:, 1].reshape(GZ.shape)}


def main():
    plane3 = radiated_fraction_into_air()
    out = {"f_hz": F_HZ, "lambda_water_m": LAM_W, "plane_wave_fraction_3d": plane3,
           "reference_3d": {}, "reference_2d": {}, "vs_lambda_3d": {}, "pinn": {}}
    print(f"f = {F_HZ} Hz, lambda_w = {LAM_W:.1f} m; plane-wave fraction into air "
          f"{plane3:.2e} ({-10*np.log10(plane3):.1f} dB)", flush=True)
    for d in (0.25, 0.5, 1.0, 2.0, 3.0, 5.0, 10.0, 20.0, 50.0):
        out["reference_3d"][str(d)] = power_fraction_reference(d, geometry="3d")
        out["reference_2d"][str(d)] = power_fraction_reference(d, geometry="2d")
        print(f"  depth {d:5.2f} m ({d/LAM_W:.3f} lambda): into air, point source "
              f"{-10*np.log10(out['reference_3d'][str(d)]):.1f} dB, line source "
              f"{-10*np.log10(out['reference_2d'][str(d)]):.1f} dB", flush=True)
    for x in np.geomspace(0.005, 3.0, 40):
        out["vs_lambda_3d"][f"{x:.4f}"] = power_fraction_reference(x * LAM_W, geometry="3d")
    # One PINN case is run and scored against the exact field.  The slope
    # cases that motivated the network were not run: a network that does not
    # reproduce the flat surface has not earned them (see README).
    for d in (3.0,):
        print(f"PINN, flat surface, depth {d} m", flush=True)
        r = train_pinn(d, steps=6000, lbfgs=300, log=lambda s: print(s, flush=True))
        ref = out["reference_2d"][str(d)]
        out["pinn"][f"flat_{d}"] = {k: v for k, v in r.items() if k not in ("grid_x", "grid_z", "field", "_nets")}
        out["pinn"][f"flat_{d}"]["reference_2d"] = ref
        print(f"    PINN fraction {r['fraction']:.3e} ({-10*np.log10(max(r['fraction'],1e-12)):.1f} dB)"
              f" vs reference {ref:.3e} ({-10*np.log10(ref):.1f} dB);  air power through the surface "
              f"{r['P_air_through_surface']:.3e} vs through the arc {r['P_air']:.3e}", flush=True)
        GX, GZ = np.meshgrid(r["grid_x"], r["grid_z"])
        pref = field_reference(GX.ravel(), GZ.ravel(), d, n=20_001).reshape(GZ.shape)
        ok = np.isfinite(r["field"]) & (np.hypot(GX, GZ - (-d)) > 1.0)
        w = ok & (GZ < 0); a_ = ok & (GZ > 0)
        err_w = np.sqrt(np.mean(np.abs(r["field"][w] - pref[w]) ** 2) / np.mean(np.abs(pref[w]) ** 2))
        err_a = np.sqrt(np.mean(np.abs(r["field"][a_] - pref[a_]) ** 2) / np.mean(np.abs(pref[a_]) ** 2))
        out["pinn"][f"flat_{d}"]["rel_field_error_water"] = float(err_w)
        out["pinn"][f"flat_{d}"]["rel_field_error_air"] = float(err_a)
        # the air field along the surface and above the source, for the page
        xs = np.linspace(-6, 6, 25)
        with __import__("torch").no_grad():
            fa = r["_nets"][1]
            pa_net = fa(__import__("torch").tensor(np.stack([xs, np.full(xs.size, 1e-3)], 1))).numpy()
        out["pinn"][f"flat_{d}"]["surface_x"] = xs.tolist()
        out["pinn"][f"flat_{d}"]["surface_abs_p_pinn"] = np.hypot(pa_net[:, 0], pa_net[:, 1]).tolist()
        out["pinn"][f"flat_{d}"]["surface_abs_p_ref"] = np.abs(field_reference(xs, np.full(xs.size, 1e-3), d, n=20_001)).tolist()
        print(f"    relative field error: water {err_w:.3f}, air {err_a:.3f}", flush=True)
        np.savez_compressed(f"results/interface_field_flat_{d}.npz", x=r["grid_x"], z=r["grid_z"],
                            pinn=r["field"], ref=pref)
    with open("results/interface_pinn.json", "w") as fh:
        json.dump(out, fh, indent=1)


if __name__ == "__main__":
    main()
