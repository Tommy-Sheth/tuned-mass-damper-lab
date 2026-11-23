"""Independent Python check of the main numbers in the report.
Run from the repository root:  python verify/python_check.py
(needs numpy and scipy)"""
import sys, os
sys.path.insert(0, os.path.dirname(__file__))
import numpy as np
from tmdcheck import *

print("== tuning ==")
dh = den_hartog(MU0)
print("Design A (Den Hartog): f = %.4f  zeta_d = %.4f  peak = %.2f" % (dh['f'], dh['zd_own'], dh['peak']))
f, zd, ratio = opt_h2()
print("Design B (numerical average-sway optimum): f = %.4f  zeta_d = %.4f  average-sway ratio = %.3f" % (f, zd, ratio))
fi, zi, pk = opt_hinf()
print("lowest-peak design: f = %.4f  zeta_d = %.4f  peak = %.2f" % (fi, zi, pk))
h2 = Sys(MU0, f, ZS0, zd)
print("Design B: L = %.2f m   peak = %.2f" % (G / h2.wd**2, hinf(h2)[0]))

print("\n== wind (50 random storms) ==")
gpk, sig_a, S0 = wind_calibration(7.5)
bare = Sys(mu=1e-9, f=0.7, zd=0.1)
dt, N = 0.05, 16000
pk = {'bare': [], 'h2': []}
for seed in range(50):
    F = synth_wind(S0, N, dt, np.random.default_rng(seed))
    for name, s in (('bare', bare), ('h2', h2)):
        y = simulate(s, F, dt, 'force')
        a = floor_accel(s, y, F, dt, 'force')[int(200/dt):int(800/dt)]
        pk[name].append(np.max(abs(a)) * 1e3 / G)
b, t = np.mean(pk['bare']), np.mean(pk['h2'])
print("mean biggest acceleration: %.2f -> %.2f milli-g  (cut %.0f%%)" % (b, t, 100 * (1 - t / b)))

print("\n== tower going out of tune ==")
f0, zd0 = f, zd
ref = hinf(Sys(MU0, f0, ZS0, zd0))[0]
for d in (-0.065, -0.03, 0.0, 0.03):
    fs = FS0 * (1 + d)
    s = Sys(MU0, f0 * FS0 / fs, ZS0, zd0, fs=fs)
    p = hinf(s)[0]
    print("%+5.1f%%  peak %.2f  (%+.0f%%)  average-sway ratio %.3f" % (100 * d, p, 100 * (p / ref - 1), s.h2_ratio()))
