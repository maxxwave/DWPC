#!/usr/bin/env python3
"""Distance scan analysis: how independent are the laser-induced DW displacements
of neighbouring wires as a function of their spacing?

usage: python3 analyze_scan.py [scan_dir]    (default ~/dwpc_scan_all)

scan_dir/gap_<rij>/run_<seed>/{laser,ref}/Benchmark3.data, from dwpc_distance_scan.sbatch.
For each gap and each seed: dx_k = x_laser,k - x_ref,k at the last output time
(same noise sequence in both runs). Across seeds, per gap:
  r_ij        Pearson correlation of dx_i and dx_j (95% CI from Fisher z)
  P(same)     probability that dx_i and dx_j have the same sign (0.5 = independent bits)
  t_k         mean(dx_k)/SE: bias of each wire (|t| > 2: not purely random)
Writes scan_summary.csv and scan_correlation.png.
"""
import glob, os, re, sys
import numpy as np

root = os.path.expanduser(sys.argv[1] if len(sys.argv) > 1 else '~/dwpc_scan_all')
rows = []
for g in glob.glob(os.path.join(root, 'gap_*')):
    tag = os.path.basename(g)[4:]
    dx = []
    for r in sorted(glob.glob(os.path.join(g, 'run_*'))):
        fl, fr = os.path.join(r, 'laser', 'Benchmark3.data'), os.path.join(r, 'ref', 'Benchmark3.data')
        if not (os.path.isfile(fl) and os.path.isfile(fr)): continue
        L, R = np.loadtxt(fl), np.loadtxt(fr)
        if L.shape != R.shape or L.ndim != 2: continue
        nw = (L.shape[1] - 5) // 2
        cols = [1 + 2*k for k in range(nw)]
        dx.append(L[-1, cols] - R[-1, cols])
    if len(dx) < 3: continue
    dx = np.array(dx); n = len(dx)
    Ly = 50e-9
    if tag == 'inf':
        d_nm = np.inf
    else:
        txt = open(glob.glob(os.path.join(g, 'run_*', 'laser', 'input'))[0]).read()
        Ly = float(re.search(r'^Ly=(\S+)', txt, re.M).group(1))
        d_nm = (float(tag) + Ly) * 1e9
    rec = dict(gap=tag, d_nm=d_nm, n=n)
    for i, j in [(0, 1), (1, 2), (0, 2)]:
        r = np.corrcoef(dx[:, i], dx[:, j])[0, 1]
        z, se = np.arctanh(np.clip(r, -0.999999, 0.999999)), 1/np.sqrt(max(n - 3, 1))
        rec['r%d%d' % (i, j)] = r
        rec['r%d%d_lo' % (i, j)], rec['r%d%d_hi' % (i, j)] = np.tanh(z - 1.96*se), np.tanh(z + 1.96*se)
        rec['psame%d%d' % (i, j)] = np.mean(np.sign(dx[:, i]) == np.sign(dx[:, j]))
    for k in range(dx.shape[1]):
        e = dx[:, k]; rec['std%d' % k] = e.std(ddof=1); rec['t%d' % k] = e.mean()/(e.std(ddof=1)/np.sqrt(n))
    rows.append(rec)

if not rows: sys.exit('no complete gap_* data in ' + root)
rows.sort(key=lambda r: r['d_nm'])
keys = list(rows[0].keys())
with open('scan_summary.csv', 'w') as f:
    f.write(','.join(keys) + '\n')
    for r in rows: f.write(','.join(str(r[k]) for k in keys) + '\n')
for r in rows:
    print('d = %7s nm  n=%3d  r01=%+.2f [%+.2f,%+.2f]  r12=%+.2f  r02=%+.2f  P(same)01=%.2f  t0=%+.1f t1=%+.1f t2=%+.1f'
          % ('inf' if np.isinf(r['d_nm']) else '%.0f' % r['d_nm'], r['n'], r['r01'], r['r01_lo'], r['r01_hi'],
             r['r12'], r['r02'], r['psame01'], r['t0'], r['t1'], r['t2']))

try:
    import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
except ImportError:
    sys.exit(0)
fin = [r for r in rows if np.isfinite(r['d_nm'])]; inf = [r for r in rows if not np.isfinite(r['d_nm'])]
ink, ink2 = '#191919', '#666666'
fig, ax = plt.subplots(figsize=(7.5, 4.5), dpi=150)
for s in ('top', 'right'): ax.spines[s].set_visible(False)
ax.grid(axis='y', color='#eeeeee'); ax.axhline(0, color='#bbbbbb', lw=0.8)
for (i, j), c, lab in [((0, 1), '#257de0', 'fire 0-1 (vecini)'), ((1, 2), '#f47230', 'fire 1-2 (vecini)'), ((0, 2), '#00a888', 'fire 0-2 (al doilea vecin)')]:
    d = np.array([r['d_nm'] for r in fin]); rr = np.array([r['r%d%d' % (i, j)] for r in fin])
    lo = np.array([r['r%d%d_lo' % (i, j)] for r in fin]); hi = np.array([r['r%d%d_hi' % (i, j)] for r in fin])
    ax.errorbar(d, rr, yerr=[rr - lo, hi - rr], color=c, marker='o', ms=6, lw=2, capsize=3, label=lab)
if inf:
    r = inf[0]; ax.axhspan(r['r01_lo'], r['r01_hi'], color='#999999', alpha=0.15, lw=0)
    ax.text(ax.get_xlim()[1], r['r01_hi'], 'fără cuplaj (95% CI)', ha='right', va='bottom', color=ink2, fontsize=9)
ax.set_xscale('log'); ax.set_xlabel('distanța centru-centru între fire (nm)')
ax.set_ylabel('corelația Δx între fire'); ax.legend(frameon=False, fontsize=9)
fig.tight_layout(); fig.savefig('scan_correlation.png'); print('figure: scan_correlation.png')
