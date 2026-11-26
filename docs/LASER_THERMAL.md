# DWPC laser / thermal extension (branch `laser-thermal`)

This branch extends DWPC with temperature-dependent magnetic parameters, a
two-temperature model (TTM) for ultrafast laser heating, fixes to the
multi-wire (coupled DW) code, and scripts for running paired-noise ensembles on
a Slurm cluster. All new options are off by default: with the shipped `input`
the single-DW output columns are bit-identical to `master` (the time column is
written with more digits and three temperature columns are appended).

## 1. Temperature-dependent material parameters

Ms(T) = Ms0 · (1 − (T/Tc)^a)^β , A(T) = A0 · (Ms/Ms0)^n , E_pin(x,T) = E_pin(x) · (Ms(T)/Ms(T0))^p

| Key | Meaning | Default |
| --- | --- | --- |
| `use_Ms_T` | enable Ms(T); `Ms` is then Ms0 (0 K) | 0 |
| `Tc` | Curie temperature (K) | 0 |
| `Ms_exponent` | β | 0 |
| `Ms_rescaling_exponent` | a (Evans et al., PRB 91, 144425 (2015), Eq. 5); 1 = classical form | 1 |
| `Aex_exponent` | n; `Aex` is then A0 (0 K); 0 = constant A | 0 |
| `pinning_exponent` | p; pinning polynomial taken at the base temperature; 0 = no scaling | 0 |

Everything that depends on Ms (muMs, H_demag, prefac1–4, A, pinning scale) is
recomputed by `stor::update_magnetic_parameters()` whenever Ms changes. Before
this branch, H_demag and muMs were evaluated once from the hard-coded defaults,
so changing `Ms`, `Ly` or `Lz` in the input did not reach them.

## 2. Time-dependent temperature

Temperature, Ms and everything derived from it are updated at every evaluation
of the right-hand side (`gradient`, `Jacobian`, `euler`), including the
intermediate Runge-Kutta stages. The thermal noise uses the instantaneous
(spin) temperature.

**Gaussian temperature pulse** (`use_T_pulse=1`):
T(t) = Temperature + `T_pulse_peak` · exp(−(t − `T_pulse_t0`)² / 2`T_pulse_sigma`²)

**Two-temperature model** (`use_TTM=1`, `src/ttm.cpp`):

    γe·Te dTe/dt = −G (Te − Tl) + P(t),        P(t) = P0 · exp(−(t − t0)² / 2σ²)
    Cl    dTl/dt =  G (Te − Tl) − Cl (Tl − T0)/τ_cool
    dm/dt        = (m_eq(Ts) − m)/τ_m          (optional, τ_m > 0)

| Key | Meaning |
| --- | --- |
| `TTM_gamma_e`, `TTM_C_l`, `TTM_G` | J m⁻³ K⁻², J m⁻³ K⁻¹, W m⁻³ K⁻¹ |
| `TTM_P0` or `TTM_fluence` | peak absorbed power (W m⁻³) or absorbed fluence (J m⁻², uniform over Lz) |
| `TTM_t0`, `TTM_sigma` | pulse centre and standard deviation (s); FWHM = 2.355 σ |
| `TTM_tau_cool` | lattice cooling to the substrate (s), 0 = off |
| `TTM_tau_m` | magnetisation relaxation time (s), 0 = Ms follows m_eq(Ts) instantly |
| `TTM_spin_temperature` | `electron` (default) or `phonon`: drives Ms(T) and the noise |
| `TTM_dt_max` | largest internal TTM step (s), default 1e-14 |

The TTM is integrated with its own adaptive RK4 sub-steps (limited by the
relaxation times and by σ/20 around the pulse) and a short history of states,
so every DW integrator stage gets T at its exact time, including RK45 stages
that step back in time and programs that restart time at 0. The run stops with
an explanation if Ms falls below 0.1 % of Ms0 (the 1D DW model is not defined
for a demagnetised wire). `use_TTM` and `use_T_pulse` are mutually exclusive.

## 3. Output

All time-series outputs (Benchmark1/2/3, spin_current, Lyapunov, multi-DW
`Current_signal.txt`) get three trailing columns `T_e(K)`, `T_l(K)`, `Ms(A/m)`.
The time column is written with 10 decimals in ns (was 4, i.e. 0.1 ps), so
fs output intervals are resolved.

## 4. Multi-wire (coupled DW) code

- `Benchmark3` writes x_i and phi_i of every wire, the applied field and the
  temperature columns (RK4 only).
- `multi_dw` RK4: stages 2–3 are evaluated at t + dt/2 and stage 4 at t + dt
  (previously all at t). Results of earlier multi-DW runs change slightly.
- `DW_coupling`: for inner wires, the interaction with the previous wire used
  the distance to the next wire; each neighbour now uses its own distance. The
  loop no longer reads `X_DW[i+1]` past the last wire.
- `multi_dw::setup` allocates `u_dw` and `j_dens_dw` (they were empty in
  `Benchmark3`); `x_coord[1]` is set only when `Nwires > 1`.
- New keys: `rij` (edge-to-edge gap; centre distance = rij + Ly),
  `DW_x0_<i>`, `DW_phi0_<i>` (initial positions/angles), `rng_seed`
  (fixed noise seed, 0 = clock).

## 5. Cluster scripts (`cluster/`)

| File | Purpose |
| --- | --- |
| `make_input.sh` | input generator for the Ni laser study (all physics parameters in one place; `RIJ=` env var, `inf` = no coupling) |
| `dwpc_ensemble.sbatch` | Slurm array: one seed per task, run with and without the pulse (same noise) |
| `dwpc_distance_scan.sbatch` | Slurm array over (gap, seed) |
| `deploy.sh`, `collect.sh` | copy/build DWPC on all nodes and gather results over `ssh -p 2222` |
| `analyze.py` | paired Δx = x_laser − x_ref statistics (mean, σ, SE, t) and figure |
| `analyze_scan.py` | inter-wire correlation of Δx vs spacing (95 % CI), same-sign probability, bias |

## 6. Verification done

- Single-DW regression: with all new options off, output identical to `master`.
- TTM vs an independent SciPy DOP853 solution: max |ΔTe| ≈ 2·10⁻⁶ K for the RK4, HEUN, RK45 and EULER DW integrators; τ_m mode within 2·10⁻⁴ Ms0.
- TTM energy balance: γe/2 (Te² − T0²) + Cl (Tl − T0) equals the absorbed energy density.
- Ms(T), Ms at the pulse maximum and the T ≥ Tc / Ms → 0 guards checked against hand calculations.
- Coupling: two wires without pinning or field at 0 K repel symmetrically; exactly aligned DWs stay at rest (unstable equilibrium), a 0.01 nm offset grows.

## 7. Known limitations / to do

- The dipolar moment `my` in `DW_coupling` is constant; it should scale with Ms(T).
- Thermal noise uses `std::default_random_engine` (minstd_rand0 in libstdc++); switch to `std::mt19937_64` before any randomness (NIST) analysis.
- Coupling is nearest-neighbour only and uses point monopoles (q = 2 Ms S): valid qualitatively, quantitatively only for centre distances of several wire widths.
- τ_m is phenomenological (to be calibrated against atomistic m(t)); pinning and A(T) scalings are assumptions to be replaced by VAMPIRE/MuMax3 values.
- The makefile still builds with `-g -pg` (profiling).
