"""Figures for the DWPC manual. Every curve is read from DWPC output files or computed
from them; reference curves are independent solutions (SciPy) of the same equations."""
import numpy as np, json, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from matplotlib.ticker import FixedLocator, FixedFormatter, NullLocator

# ---- house style: sizes match the LaTeX page (text width 16 cm = 6.30 in, 11 pt body)
FULL, HALF = 6.3, 3.6
C = ['#257de0', '#f47230', '#00a884', '#c98a00']          # blue, orange, green, ochre
ink, ink2, grid = '#191919', '#5f6368', '#ececec'
plt.rcParams.update({
    'font.family': 'serif', 'font.serif': ['cmr10'], 'axes.formatter.use_mathtext': True, 'axes.unicode_minus': True,
    'mathtext.fontset': 'cm', 'font.size': 9, 'axes.labelsize': 9, 'xtick.labelsize': 8,
    'ytick.labelsize': 8, 'legend.fontsize': 8, 'axes.edgecolor': '#9a9a9a', 'axes.linewidth': 0.7,
    'axes.labelcolor': ink, 'xtick.color': ink2, 'ytick.color': ink2, 'xtick.major.width': 0.7,
    'ytick.major.width': 0.7, 'lines.linewidth': 1.4, 'legend.frameon': False,
    'legend.handlelength': 1.8, 'pdf.fonttype': 42, 'savefig.bbox': 'tight', 'savefig.pad_inches': 0.03})
OUT = '/tmp/manual/figures/'; R = {}

def style(a):
    a.grid(color=grid, lw=0.6); a.set_axisbelow(True)
    for s in ['top', 'right']: a.spines[s].set_visible(False)
def panel(a, s, x=-0.13, y=1.02):
    a.text(x, y, s, transform=a.transAxes, fontsize=10, fontname='cmb10', va='bottom', ha='left', color=ink)
def tag(a, s):
    a.text(0.985, 0.95, s, transform=a.transAxes, ha='right', va='top', fontsize=8, color=ink,
           bbox=dict(boxstyle='round,pad=0.25', fc='white', ec='none', alpha=0.9))
def load(f):
    return np.array([[float(v) for v in l.split()] for l in open(f) if not l.startswith('#') and l.strip()])
trap = getattr(np, 'trapezoid', None) or np.trapz
def save(fig, name):
    fig.savefig(OUT + name + '.pdf'); fig.savefig('/tmp/rc/prev_' + name + '.png', dpi=110); plt.close(fig)

# ---- A: driven DW at 0 K and 100 K
fig, ax = plt.subplots(2, 1, figsize=(FULL, 3.5), sharex=True)
for a, (d, lab, p) in zip(ax, [('b1_0K', '$T$ = 0 K', '(a)'), ('b1_100K', '$T$ = 100 K', '(b)')]):
    r = load('/tmp/bm/' + d + '/output'); style(a)
    for s in (1, -1): a.axhline(s * 62.78, color=ink2, lw=0.7, ls=(0, (3, 3)))
    a.plot(r[:, 0], r[:, 1], color=C[0], lw=1.1)
    a.set_ylabel('$x$ (nm)'); a.set_ylim(-150, 150); panel(a, p, x=-0.075); a.set_title(lab, loc='right', fontsize=8.5, color=ink, pad=3)
ax[1].set_xlabel('time (ns)'); ax[1].set_xlim(0, 20)
fig.tight_layout(h_pad=0.6); save(fig, 'driven_dw')

# ---- B: time-step convergence
dts = ['8e-12', '4e-12', '2e-12', '1e-12', '5e-13', '2.5e-13']; d = np.array([float(x) for x in dts])
xr = load('/tmp/bm/b2b/ref/output')[-1, 1]; conv = {}
fig, a = plt.subplots(figsize=(HALF, 2.9)); style(a)
for i, (integ, order) in enumerate([('HEUN', 2), ('RK4', 4), ('RK45', 5)]):
    e = np.array([abs(load('/tmp/bm/b2b/%s_%s/output' % (integ, dt))[-1, 1] - xr) for dt in dts]); conv[integ] = e.tolist()
    a.loglog(d * 1e12, e, 'o-', color=C[i], ms=4, label='%s (order %d)' % (integ, order))
    a.loglog(d * 1e12, e[0] * (d / d[0]) ** order, color=C[i], lw=0.7, ls=(0, (4, 3)))
a.axhline(1e-9, color=ink2, lw=0.6); a.text(0.27, 7.5e-10, 'round-off floor', fontsize=7.5, color=ink2, va='top')
a.set_ylim(3e-10, 3e-1); a.set_xlabel(r'time step $\Delta t$ (ps)'); a.set_ylabel(r'$|x-x_{\rm ref}|$ at 4 ns (nm)')
a.xaxis.set_major_locator(FixedLocator([0.25, 0.5, 1, 2, 4, 8])); a.xaxis.set_major_formatter(FixedFormatter(['0.25', '0.5', '1', '2', '4', '8'])); a.xaxis.set_minor_locator(NullLocator())
a.legend(loc='upper left'); fig.tight_layout(); save(fig, 'convergence'); R['conv'] = conv; R['dts'] = dts

# ---- C: equipartition in one well
a2, a4 = -1.2866e-6, 1.632e8; kB = 1.380649e-23; k = -4 * a2; x0 = np.sqrt(-a2 / (2 * a4))
E = lambda x: a2 * x * x + a4 * x ** 4
def exact(Tk):
    b = 1 / (kB * Tk); u = np.linspace(0, 4 * x0 * 1e9, 400001); w = np.exp(-b * (E(u * 1e-9) - E(x0)))
    Z = trap(w, u); m1 = trap(u * w, u) / Z; m2 = trap(u * u * w, u) / Z
    return np.sqrt(m2 - m1 ** 2)
eq = {}
fig, a = plt.subplots(figsize=(HALF, 2.9)); style(a)
Ts = np.linspace(2, 22, 60)
a.plot(Ts, [exact(t) for t in Ts], color=ink, lw=1.2, label='Boltzmann, full $E(x)$')
a.plot(Ts, np.sqrt(kB * Ts / k) * 1e9, color=ink2, lw=1, ls=(0, (4, 3)), label="harmonic $\\sqrt{k_BT/E''(x_0)}$")
for i, integ in enumerate(['RK4', 'HEUN']):
    pts = []
    for Tk in [5, 10, 20]:
        r = load('/tmp/bm/b3/%s_%dK/output' % (integ, Tk)); xs = r[r[:, 0] > 20, 1]; pts.append((Tk, xs.std()))
        eq['%s_%d' % (integ, Tk)] = [xs.std(), exact(Tk), np.sqrt(kB * Tk / k) * 1e9]
    p = np.array(pts); a.plot(p[:, 0] + (i - 0.5) * 0.5, p[:, 1], 'o', color=C[i], ms=5, mec='white', mew=0.6, label='DWPC, ' + integ)
a.set_xlabel('temperature $T$ (K)'); a.set_ylabel(r'$\sigma_x$ in one well (nm)'); a.set_xlim(0, 23)
a.legend(loc='upper left'); fig.tight_layout(); save(fig, 'equipartition'); R['eq'] = eq

# ---- D: two-temperature model vs SciPy
g, Cl, G, P0, t0, s, T0 = 6e3, 3.9e6, 8e17, 1.76356e22, 0.5e-12, 2.548e-14, 300
Tc, be, al = 635, 0.341, 2.322
meq = lambda T: (1 - (min(T, Tc) / Tc) ** al) ** be
fig, ax = plt.subplots(2, 1, figsize=(FULL, 4.1), sharex=True, gridspec_kw=dict(height_ratios=[1.15, 1]))
for a in ax: style(a)
ttm = {}
for j, tm in enumerate([200e-15, 300e-15]):
    r = load('/tmp/bm/b4/tau_%se-15/output' % ('200' if j == 0 else '300'))
    t = r[:, 0] * 1e-9; Te, Tl, Ms = r[:, -3], r[:, -2], r[:, -1]
    def f(tt, y):
        P = P0 * np.exp(-(tt - t0) ** 2 / (2 * s * s)); return [(-G * (y[0] - y[1]) + P) / (g * y[0]), G * (y[0] - y[1]) / Cl, (meq(y[0]) - y[2]) / tm]
    sol = solve_ivp(f, (0, t[-1]), [T0, T0, meq(T0)], t_eval=t, rtol=1e-11, atol=1e-10, max_step=2e-15, method='DOP853')
    ttm['tau_%d' % round(tm * 1e15)] = [float(abs(Te - sol.y[0]).max()), float(abs(Tl - sol.y[1]).max()), float(abs(Ms / 509407.45 - sol.y[2]).max()), float(Te.max()), float(Ms.min() / (509407.45 * meq(300)))]
    tp = (t - t0) * 1e12
    sel = np.r_[np.arange(0, len(tp), 15)[tp[::15] < 1.0], np.where(tp >= 1.0)[0][::90]]
    if j == 0:
        ax[0].plot(tp, Te, color=C[1], label='$T_e$ (DWPC)'); ax[0].plot(tp, Tl, color=C[2], label='$T_l$ (DWPC)')
        ax[0].plot(tp[sel], sol.y[0][sel], 'o', mfc='none', color=ink, ms=3.2, mew=0.6, label='SciPy DOP853')
        ax[0].plot(tp[sel], sol.y[1][sel], 'o', mfc='none', color=ink, ms=3.2, mew=0.6)
        ax[0].axhline(Tc, color=ink2, lw=0.7, ls=(0, (3, 3))); ax[0].text(4.97, Tc + 6, '$T_c$ = 635 K', ha='right', va='bottom', fontsize=7.5, color=ink2)
    ax[1].plot(tp, Ms / 1e3, color=C[0], ls='-' if j == 0 else (0, (5, 2)), label=r'$\tau_m$ = %d fs (DWPC)' % round(tm * 1e15))
    ax[1].plot(tp[sel], sol.y[2][sel] * 509.40745, 'o', mfc='none', color=ink, ms=3.2, mew=0.6, label='SciPy DOP853' if j == 0 else None)
ax[0].set_ylabel('temperature (K)'); ax[0].legend(loc='lower left', bbox_to_anchor=(0.04, 1.0), ncol=3, borderaxespad=0.1); panel(ax[0], '(a)', x=-0.075)
ax[1].set_ylabel('$M_s$ (kA/m)'); ax[1].set_xlabel('time after the pulse maximum (ps)'); ax[1].legend(loc='lower right'); ax[1].set_xlim(-0.5, 5); panel(ax[1], '(b)', x=-0.075)
fig.tight_layout(h_pad=0.8); save(fig, 'ttm_ni'); R['ttm'] = ttm
r = load('/tmp/bm/b4/tau_200e-15/output'); Te, Tl = r[-1, -3], r[-1, -2]
R['energy'] = [g / 2 * (Te ** 2 - T0 ** 2) + Cl * (Tl - T0), P0 * s * np.sqrt(2 * np.pi)]

# ---- E: repulsion between DWs in two wires
fig, a = plt.subplots(figsize=(HALF, 2.9)); style(a); cp = {}
for i, gap in enumerate(['50e-9', '100e-9', '150e-9', '250e-9']):
    r = load('/tmp/bm/b6/gap_%s/Benchmark3.data' % gap); sep = r[:, 1] - r[:, 3]; cp[gap] = float(sep[-1])
    a.semilogy(r[:, 0], sep, color=C[i], label='%d nm' % round(float(gap) * 1e9))
a.set_xlabel('time (ns)'); a.set_ylabel('DW separation $x_0-x_1$ (nm)'); a.legend(title='wire spacing', title_fontsize=8, loc='upper left'); a.set_xlim(0, 2)
fig.tight_layout(); save(fig, 'coupling'); R['coupling'] = cp
json.dump(R, open('/tmp/manual/results.json', 'w'), indent=1, default=float)

# ---- F: sine/square classification (DWPC Signal.out at 0 K and 100 K)
Nv = 10
def rc_load(D):
    dd = np.loadtxt(D + 'sine_square_init'); x, y = dd[:, 0], dd[:, 1]; n = len(x) // 2
    S = np.fromfile(D + 'Signal.out', sep=' ').reshape(n, Nv); V = np.fromfile(D + 'Valid_Signal.out', sep=' ').reshape(n, Nv)
    W = np.linalg.solve(S.T @ S + 1e-6 * np.eye(Nv), S.T @ y[:n])
    acc = lambda p, t: int(((p > 0.5).astype(int) == t.astype(int)).sum())
    return x, y[:n], S, S @ W, {'train': acc(S @ W, y[:n]), 'valid': acc(V @ W, y[n:]), 'n': n}
x, yt, S0, p0, r0 = rc_load('/tmp/rc/sine0/'); _, _, S1, p1, r1 = rc_load('/tmp/rc/sine/')
json.dump({'0K': r0, '100K': r1}, open('/tmp/rc/sine_res.json', 'w'))
kk = np.arange(160)
fig, ax = plt.subplots(4, 1, figsize=(FULL, 5.6), sharex=True, gridspec_kw=dict(height_ratios=[0.8, 1, 1, 1]))
a = ax[0]; style(a)
for i in kk:
    if yt[i] > 0.5: a.axvspan(i - 0.5, i + 0.5, color=C[1], alpha=0.13, lw=0)
a.step(kk, x[:160], where='mid', color=ink, lw=1.0); a.set_ylabel('input $u_k$'); a.set_ylim(-1.3, 1.3); panel(a, '(a)', x=-0.075)
a.text(55.5, 1.1, 'square', fontsize=7.5, color=C[1], ha='center', va='bottom')
a = ax[1]; style(a)
for j in range(Nv): a.plot(kk, S0[:160, j], color=C[0], lw=0.7, alpha=0.6)
a.set_ylabel('state $s_{ki}$ (nm)'); panel(a, '(b)', x=-0.075); a.set_title('10 virtual nodes, $T$ = 0 K', loc='right', fontsize=8.5, color=ink, pad=3)
for a, p, T, pl in [(ax[2], p0, '0 K', '(c)'), (ax[3], p1, '100 K', '(d)')]:
    style(a); a.step(kk, yt[:160], where='mid', color=ink2, lw=0.9, label='target')
    a.plot(kk, p[:160], color=C[1], lw=1.1, label=r'readout $\mathbf{W}\cdot\mathbf{s}_k$')
    a.axhline(0.5, color=ink2, lw=0.6, ls=(0, (3, 3))); a.set_ylabel('output'); a.set_ylim(-0.45, 1.3)
    panel(a, pl, x=-0.075); a.set_title('readout, $T$ = ' + T, loc='right', fontsize=8.5, color=ink, pad=3)
ax[2].legend(loc='upper right', ncol=2); ax[3].set_xlabel('sample index $k$'); ax[3].set_xlim(-1, 160)
fig.tight_layout(h_pad=0.5); save(fig, 'rc_sine')

# ---- G: spoken digits (DWPC Processed_spoken_signal.txt at 0 K)
Rs = json.load(open('/tmp/rc/spoken_dt1e-12_0K.json')); CM = np.array(Rs['cm'])
raw = np.loadtxt('/tmp/rc/spoken/Signal_in.data'); U = raw[:, :1024]; lab = raw[:, 1024].astype(int)
X = np.loadtxt('/tmp/rc/spoken0/Processed_spoken_signal.txt')
fig = plt.figure(figsize=(FULL, 3.0)); gs = fig.add_gridspec(2, 2, width_ratios=[1.35, 1], wspace=0.32, hspace=0.18)
i0 = np.where(lab == 3)[0][0]; kq = np.arange(96)
a = fig.add_subplot(gs[0, 0]); style(a); a.plot(kq, U[i0, :96], color=ink2, lw=0.9, marker='o', ms=2)
a.set_ylabel('input $u_{ik}$'); a.tick_params(labelbottom=False); panel(a, '(a)', x=-0.16)
for c in range(1, 6): a.axvline(16 * c - 0.5, color=grid, lw=1.2, zorder=0)
b = fig.add_subplot(gs[1, 0], sharex=a); style(b); b.plot(kq, X[i0, :96], color=C[0], lw=0.9, marker='o', ms=2)
for c in range(1, 6): b.axvline(16 * c - 0.5, color=grid, lw=1.2, zorder=0)
b.set_ylabel('$X_{ik}$ (nm)'); b.set_xlabel('feature index $k$ (16 per channel)'); b.set_xlim(-1, 96)
c_ = fig.add_subplot(gs[:, 1]); c_.imshow(CM, cmap='Blues', vmin=0, vmax=100)
c_.set_xticks(range(10)); c_.set_yticks(range(10)); c_.set_xlabel('predicted digit'); c_.set_ylabel('true digit')
for s_ in c_.spines.values(): s_.set_visible(False)
c_.tick_params(length=0)
for i in range(10):
    for j in range(10):
        if CM[i, j]: c_.text(j, i, CM[i, j], ha='center', va='center', fontsize=6.5, color='white' if CM[i, j] > 55 else ink)
panel(c_, '(b)', x=-0.2)
save(fig, 'rc_spoken')
print('ok', R['conv']['RK4'][:2], r0, r1)
