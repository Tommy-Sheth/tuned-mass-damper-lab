"""
tmdcheck.py -- the Python model behind the numbers in this report.

Everything here is written from the equations in Section 1 of the report.  Units: SI.  Building first mode (M, k, c) with a pendulum damper
(m, kd, cd);  x = floor, y = ball (both from the ground), u = y - x.
"""
import numpy as np
from scipy import linalg, signal, optimize

G = 9.80665
M0, FS0, ZS0, MASS0 = 39_655e3, 0.150, 0.015, 660e3
MU0 = MASS0 / M0
GAMMA1 = 1.24678521          # first-mode participation factor, shear beam, zd/H = 0.87


# ----------------------------------------------------------------------------- model ----------
class Sys:
    def __init__(self, mu=MU0, f=1.0, zs=ZS0, zd=0.0, M=M0, fs=FS0, Gamma=GAMMA1):
        self.M, self.mu, self.f, self.zs, self.zd, self.Gamma = M, mu, f, zs, zd, Gamma
        self.ws = 2 * np.pi * fs
        self.k = M * self.ws ** 2
        self.c = 2 * zs * M * self.ws
        self.m = mu * M
        self.wd = f * self.ws
        self.kd = self.m * self.wd ** 2
        self.cd = 2 * zd * self.m * self.wd

    # state z = [x, y, xdot, ydot]
    def AB(self):
        M, m, k, c, kd, cd = self.M, self.m, self.k, self.c, self.kd, self.cd
        Mm = np.diag([M, m])
        K = np.array([[k + kd, -kd], [-kd, kd]])
        C = np.array([[c + cd, -cd], [-cd, cd]])
        Mi = np.linalg.inv(Mm)
        A = np.block([[np.zeros((2, 2)), np.eye(2)], [-Mi @ K, -Mi @ C]])
        return A, Mm, Mi

    def frf(self, beta):
        """X / (F0/k) for a force on the building, beta = Omega / omega_s."""
        b = np.asarray(beta, dtype=float)
        mu, f, zs, zd = self.mu, self.f, self.zs, self.zd
        d = f ** 2 - b ** 2 + 2j * zd * f * b
        return d / ((1 - b ** 2 + 2j * zs * b) * d - mu * b ** 2 * (f ** 2 + 2j * zd * f * b))

    def ratio_Y_over_X(self, beta):
        b = np.asarray(beta, dtype=float)
        f, zd = self.f, self.zd
        return (f ** 2 + 2j * zd * f * b) / (f ** 2 - b ** 2 + 2j * zd * f * b)

    def eig(self):
        A, _, _ = self.AB()
        lam = np.linalg.eigvals(A)
        lam = lam[lam.imag > 0]
        wn = np.abs(lam)
        zeta = -lam.real / wn
        o = np.argsort(wn)
        return wn[o] / self.ws, zeta[o]

    # --- white-noise variances (unit two-sided force intensity) ---------------------------------
    def lyap(self):
        A, _, Mi = self.AB()
        B = np.array([[0], [0], [Mi[0, 0]], [0]])
        P = linalg.solve_continuous_lyapunov(A, -B @ B.T)
        return P

    def h2_ratio(self):
        bare = Sys(mu=0.0, f=1.0, zs=self.zs, zd=0.0, M=self.M, fs=self.ws / 2 / np.pi)
        var_bare = 1.0 / (2 * bare.k * bare.c)
        return np.sqrt(self.lyap()[0, 0] / var_bare)

    def power_share(self):
        P = self.lyap()
        A, _, Mi = self.AB()
        pb = self.c * P[2, 2]
        u_dot_var = P[2, 2] + P[3, 3] - 2 * P[2, 3]
        pt = self.cd * u_dot_var
        return pt / (pt + pb)


def hinf(s, lo=0.5, hi=1.6, n=6001):
    b = np.linspace(lo, hi, n)
    a = np.abs(s.frf(b))
    i = int(np.argmax(a))
    r = optimize.minimize_scalar(lambda x: -abs(s.frf(x)), bounds=(b[max(i - 1, 0)], b[min(i + 1, n - 1)]),
                                 method="bounded", options=dict(xatol=1e-12))
    return -r.fun, r.x


def bare_peak(zs=ZS0):
    return 1 / (2 * zs * np.sqrt(1 - zs ** 2))


# -------------------------------------------------------------------------- closed forms -------
def den_hartog(mu):
    f = 1 / (1 + mu)
    zd_own = np.sqrt(3 * mu / (8 * (1 + mu)))            # damping measured against omega_d (this report)
    zd_text = np.sqrt(3 * mu / (8 * (1 + mu) ** 3))      # textbook: measured against omega_s
    s = np.sqrt(mu / (2 + mu))
    gP, gQ = np.sqrt((1 - s) / (1 + mu)), np.sqrt((1 + s) / (1 + mu))
    return dict(f=f, zd_own=zd_own, zd_text=zd_text, gP=gP, gQ=gQ, peak=np.sqrt(1 + 2 / mu))


def warburton(mu):
    f = np.sqrt(1 + mu / 2) / (1 + mu)
    zd = np.sqrt(mu * (1 + 0.75 * mu) / (4 * (1 + mu) * (1 + mu / 2)))
    zeq = 0.25 * np.sqrt(4 * mu * (1 + mu) / (4 + 3 * mu))
    return dict(f=f, zd=zd, zeq=zeq)


def opt_h2(mu=MU0, zs=ZS0):
    fun = lambda p: Sys(mu, p[0], zs, p[1]).h2_ratio()
    r = optimize.minimize(fun, [0.985, 0.065], method="Nelder-Mead", options=dict(xatol=1e-9, fatol=1e-12, maxiter=2000))
    return r.x[0], r.x[1], r.fun


def opt_hinf(mu=MU0, zs=ZS0):
    fun = lambda p: hinf(Sys(mu, p[0], zs, p[1]))[0]
    r = optimize.minimize(fun, [0.98, 0.08], method="Nelder-Mead", options=dict(xatol=1e-8, fatol=1e-10, maxiter=2000))
    return r.x[0], r.x[1], r.fun


# --------------------------------------------------------------------------- wind --------------
LU = 200.0 / 30.0


def S_vk(freq):
    return (1 + 70.8 * (freq * LU) ** 2) ** (-5 / 6)


def accel_var(s, S0=1.0, kind="abs"):
    """Variance of floor acceleration for force spectrum S0 * S_vk(f) (one-sided in f)."""
    fr = np.concatenate([np.logspace(-3, np.log10(0.5 * s.ws / 2 / np.pi), 4000),
                         np.linspace(0.5, 1.5, 6000) * s.ws / 2 / np.pi,
                         np.logspace(np.log10(1.5 * s.ws / 2 / np.pi), np.log10(5), 4000)])
    fr = np.unique(fr)
    beta = 2 * np.pi * fr / s.ws
    Hx = s.frf(beta) / s.k                       # displacement per newton
    Ha = (2 * np.pi * fr) ** 2 * Hx
    H = Ha if kind == "abs" else Hx
    return np.trapezoid(np.abs(H) ** 2 * S0 * S_vk(fr), fr)


def davenport(nu, T):
    a = np.sqrt(2 * np.log(nu * T))
    return a + 0.5772 / a


def wind_calibration(target_mg=7.5, T=600.0):
    bare = Sys(mu=0.0)
    gpk = davenport(FS0, T)
    sig_a = target_mg * 1e-3 * G / gpk
    S0 = sig_a ** 2 / accel_var(bare, 1.0)
    return gpk, sig_a, S0


def synth_wind(S0, n, dt, rng):
    """Gaussian force, one-sided spectrum S0*S_vk(f) in N^2/Hz, random-phase spectral representation."""
    fr = np.fft.rfftfreq(n, dt)
    df = fr[1] - fr[0]
    a = np.sqrt(2.0 * S0 * S_vk(fr) * df)
    a[0] = 0.0
    X = (n / 2.0) * a * np.exp(1j * rng.uniform(0, 2 * np.pi, len(fr)))
    return np.fft.irfft(X, n)


def simulate(s, u, dt, kind="force", x0=None):
    """Exact (first-order-hold) time response. kind='force' (N on the building) or 'ground' (m/s^2)."""
    A, Mm, Mi = s.AB()
    if kind == "force":
        B = np.array([[0], [0], [Mi[0, 0]], [0]])
    else:
        B = np.array([[0], [0], [-s.Gamma], [-1.0]])
    sysd = signal.cont2discrete((A, B, np.eye(4), np.zeros((4, 1))), dt, method="foh")
    Ad, Bd, Cd, Dd, _ = sysd
    # use dlsim (handles FOH state shift internally)
    t, y, xs = signal.dlsim((Ad, Bd, Cd, Dd, dt), u.reshape(-1, 1), x0=x0)
    return y   # columns x, y, xdot, ydot


def floor_accel(s, y, u, dt, kind="force"):
    """Absolute floor acceleration from the equation of motion."""
    x, yy, xd, yd = y.T
    if kind == "force":
        return (u - s.k * x - s.c * xd + s.kd * (yy - x) + s.cd * (yd - xd)) / s.M
    return (-s.k * x - s.c * xd + s.kd * (yy - x) + s.cd * (yd - xd)) / s.M - s.Gamma * u


# ---------------------------------------------------------------------- ground motions ---------
def pulse_motion(dt=0.02, T=60.0, Tp=5.0, gamma=2.0, nu=0.0, t0=10.0):
    t = np.arange(0, T, dt)
    fp = 1 / Tp
    tau = t - t0
    v = np.where(np.abs(tau) <= gamma / (2 * fp),
                 0.5 * (1 + np.cos(2 * np.pi * fp * tau / gamma)) * np.cos(2 * np.pi * fp * tau + nu), 0.0)
    return np.gradient(v, dt)


def long_period_motion(seed, dt=0.02, T=300.0, f_lo=0.08, f_hi=0.5, rise=30.0, strong=120.0, tau_d=40.0):
    rng = np.random.default_rng(seed)
    n = int(T / dt)
    w = rng.standard_normal(n)
    sos = signal.butter(4, [f_lo, f_hi], btype="band", fs=1 / dt, output="sos")
    a = signal.sosfiltfilt(sos, w)
    t = np.arange(n) * dt
    env = np.where(t < rise, t / rise, np.where(t < rise + strong, 1.0, np.exp(-(t - rise - strong) / tau_d)))
    return a * env


def read_at2(path):
    vals = []
    with open(path) as fh:
        lines = fh.readlines()
    for ln in lines[4:]:
        vals += [float(x) for x in ln.split()]
    return np.array(vals) * G
