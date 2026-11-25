#!/usr/bin/env python3
"""Paired analysis of the DWPC 3-wire laser ensemble.

usage: python3 analyze.py [results_dir]   (default ~/dwpc_ens_all)

For every seed, dx = x_laser - x_ref is computed for each wire; the two runs
share the same thermal-noise sequence, so they coincide before the pulse.
Writes dx_summary.csv (mean, std, standard error, t-statistic vs time after the
pulse) and dx_final.csv (per-seed dx at the last output time), and a figure
dwpc_ensemble.png if matplotlib is available.
"""
import glob, os, re, sys
import numpy as np

root = os.path.expanduser(sys.argv[1] if len(sys.argv) > 1 else '~/dwpc_ens_all')
runs = sorted(glob.glob(os.path.join(root, 'run_*')))
seeds, D, t = [], [], None
t_pulse = None
for r in runs:
    fl, fr = os.path.join(r, 'laser', 'Benchmark3.data'), os.path.join(r, 'ref', 'Benchmark3.data')
    if not (os.path.isfile(fl) and os.path.isfile(fr)):
        print('skip (incomplete):', r); continue
    L, R = np.loadtxt(fl), np.loadtxt(fr)
    if L.shape != R.shape:
        print('skip (length mismatch):', r); continue
    nw = (L.shape[1] - 5) // 2                      # time, (x,phi)*nw, H, Te, Tl, Ms
    cols = [1 + 2*k for k in range(nw)]
    if t_pulse is None:
        txt = open(os.path.join(r, 'laser', 'input')).read()
        t_pulse = float(re.search(r'^TTM_t0=(\S+)', txt, re.M).group(1)) * 1e9   # ns
    pre = L[:, 0] < t_pulse - 1e-6
    mism = np.abs(L[pre][:, cols] - R[pre][:, cols]).max() if pre.any() else 0.0
    if mism > 1e-9:
        print('warning: laser/ref differ before the pulse in', r, mism)
    seeds.append(int(r.split('_')[-1])); D.append(L[:, cols] - R[:, cols]); t = L[:, 0]

if not D:
    sys.exit('no complete runs found in ' + root)
D = np.array(D)                                     # (n_seeds, n_times, n_wires)
n, _, nw = D.shape
m = t >= t_pulse
tt = t[m] - t_pulse
mu, sd = D[:, m, :].mean(0), D[:, m, :].std(0, ddof=1)
se = sd / np.sqrt(n)
tstat = np.divide(mu, se, out=np.zeros_like(mu), where=se > 0)

hdr = 'time_after_pulse_ns,' + ','.join('mean_dx%d_nm,std_dx%d_nm,se_dx%d_nm,t_dx%d' % (k, k, k, k) for k in range(nw))
np.savetxt('dx_summary.csv', np.column_stack([tt] + [np.column_stack([mu[:, k], sd[:, k], se[:, k], tstat[:, k]]) for k in range(nw)]),
           delimiter=',', header=hdr, comments='')
np.savetxt('dx_final.csv', np.column_stack([seeds, D[:, -1, :]]), delimiter=',',
           header='seed,' + ','.join('dx%d_nm' % k for k in range(nw)), comments='', fmt=['%d'] + ['%.4f'] * nw)

print('%d paired seeds, pulse at %.3f ns, %d wires' % (n, t_pulse, nw))
for k in range(nw):
    e = D[:, -1, k]
    print('wire %d at %.3f ns after pulse: mean %+.2f nm, std %.2f nm, SE %.2f nm, t = %+.2f, %d/%d positive'
          % (k, tt[-1], e.mean(), e.std(ddof=1), e.std(ddof=1)/np.sqrt(n), tstat[-1, k], (e > 0).sum(), n))
print('|t| > ~2 suggests a systematic (non-random) laser-induced shift for that wire.')

try:
    import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
except ImportError:
    sys.exit(0)
ink, ink2, c1 = '#191919', '#666666', '#257de0'
fig, ax = plt.subplots(nw, 2, figsize=(9.5, 3*nw), dpi=150, gridspec_kw=dict(width_ratios=[3, 1]), sharey='row', squeeze=False)
for k in range(nw):
    A, H = ax[k, 0], ax[k, 1]
    for a in (A, H):
        a.grid(axis='y', color='#eeeeee'); a.set_axisbelow(True); a.axhline(0, color='#bbbbbb', lw=0.8)
        for s in ('top', 'right'): a.spines[s].set_visible(False)
    for j in range(n):
        A.plot(tt, D[j, m, k], color='#9fc3ee', lw=0.6, alpha=0.6 if n > 30 else 1.0)
    A.fill_between(tt, mu[:, k]-sd[:, k], mu[:, k]+sd[:, k], color=c1, alpha=0.12, lw=0)
    A.plot(tt, mu[:, k], color=c1, lw=2.2)
    A.set_ylabel('Δx wire %d (nm)' % k); A.set_xlim(0, tt[-1])
    e = D[:, -1, k]
    H.hist(e, bins=min(30, max(8, n//4)), orientation='horizontal', color=c1, edgecolor='white')
    H.set_title('Δx at %.2f ns' % tt[-1], loc='left', fontsize=10, color=ink)
    H.text(0.98, 0.04, 'mean %+.1f nm\nσ %.1f nm\nt = %+.2f' % (e.mean(), e.std(ddof=1), tstat[-1, k]),
           transform=H.transAxes, ha='right', va='bottom', fontsize=9, color=ink)
ax[-1, 0].set_xlabel('time after pulse (ns)')
fig.suptitle('DWPC ensemble: x_laser − x_ref, %d paired seeds' % n, x=0.01, ha='left', fontsize=11, color=ink)
fig.tight_layout(); fig.savefig('dwpc_ensemble.png')
print('figure: dwpc_ensemble.png')
