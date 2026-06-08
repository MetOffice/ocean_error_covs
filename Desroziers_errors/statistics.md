# Statistics: Sample Collection and Estimators

This document describes how observations are pooled into samples, how the
accumulated statistics are stored in `binned.nc`, and how the covariance
estimators (including bias variants) are derived from those statistics.

---

## 1. Sample collection

### 1.1 Per-timestep accumulation

The tool iterates over all input files in chronological order. Each file
corresponds to one time step (e.g. one day or one month, as set by `data_freq`
in the observation `.ini` file). For every time step:

1. All required variables are loaded from the NetCDF file.
2. For **each predictor bin** $s = 0, 1, \ldots, n_\text{classes}-1$:
   - The `mask_func` (or default predictor mask) is applied to select
     observations belonging to bin $s$.
   - If the number of selected observations exceeds `max_obs_per_day`, a
     random subset of that size is drawn without replacement.
3. For each ordered pair $(i, j)$ where observation $i$ comes from the first
   observation type and $j$ from the second (same type for auto-covariances):
   - The great-circle distance between $i$ and $j$ is computed using the
     haversine formula and the result is assigned to the nearest distance bin
     $k$ (see §1.2).
   - The innovation vectors at $i$ and $j$ are accumulated into the running
     sums described in §1.3.

This process repeats for every time step and the running sums grow monotonically.
All arrays in `binned.nc` store the **total accumulated sums** over the entire
date range, not per-timestep values.

### 1.2 Spatial distance bins

The correlation grid is a 1-D array of evenly spaced distances:

$$
\mathbf{d} = [0,\; \Delta d, 2\Delta d, \ldots, d_\text{max} - \Delta d]
\quad \text{[km]}
$$

set by `max_distance` and `distance_interval` in the `[Grid]` section of the
observation `.ini` file. Each pair $(i, j)$ is assigned to the bin whose
centre is closest to the haversine distance between the two observations.
Pairs at zero separation (i.e. the same location) fall in bin $k=0$ and
contribute to the variance estimate.

### 1.3 Accumulated sums

For predictor bin $s$ and distance bin $k$, with sample count
$c_{s,k}$ = `total_obs_bin[s, k]`, the following sums are accumulated:

| Symbol | `binned.nc` variable | Definition |
|--------|----------------------|------------|
| $c_{s,k}$ | `total_obs_bin` | Number of observation pairs $(i,j)$ in bin $(s,k)$ |
| $S^i_{ob}$ | `d_ob_i` | $\displaystyle\sum d_{ob}(i) = \sum [y - \mathbf{H}(\mathbf{x}_b)]_i$ |
| $S^j_{ob}$ | `d_ob_j` | $\displaystyle\sum d_{ob}(j) = \sum [y - \mathbf{H}(\mathbf{x}_b)]_j$ |
| $S_{oa}$ | `d_oa` | $\displaystyle\sum d_{oa}(j) = \sum [y - \mathbf{H}(\mathbf{x}_a)]_j$ |
| $S_{ab}$ | `d_ab` | $\displaystyle\sum d_{ab}(j) = \sum [\mathbf{H}(\mathbf{x}_a) - \mathbf{H}(\mathbf{x}_b)]_j$ |
| $E_R$ | `R` | $\displaystyle\sum d_{ob}(i)\cdot d_{oa}(j)$ |
| $E_{HBH}$ | `HBH` | $\displaystyle\sum d_{ob}(i)\cdot d_{ab}(j)$ |
| $E_{R+HBH}$ | `R_HBH` | $\displaystyle\sum d_{ob}(i)\cdot d_{ob}(j)$ |

The subscript $i$ denotes the "starting" observation of the pair and $j$ the
"ending" observation. For auto-covariances both come from the same observation
type; for cross-covariances they come from different types.

Only the sums corresponding to the `estimate` options selected in `config.ini`
are accumulated and saved.

---

## 2. Covariance estimators

Covariance estimates are derived from the accumulated sums once all time steps
have been processed. Bins with $c_{s,k} \le c_\text{min}$ (default 100, set by
`ss_lim` in `CorrEstimator.do_compute`) are masked to `NaN`.

### 2.1 Bias-corrected (default) estimate

The standard bias-corrected (centred) sample cross-covariance estimate removes
the contribution of non-zero mean innovations:

$$
\hat{C}[v]_{s,k} = \frac{1}{c - 1}
\left(E_v - \frac{S^i_{ob}S^j_w}{c}\right)
$$

where $v \in \{R, HBH, R+HBH\}$ and $S^j_w$ is the corresponding
accumulation sum at the $j$ point ($S_{oa}$, $S_{ab}$, or $S^j_{ob}$
respectively). This is the primary estimate written to `cov_R`, `cov_HBH`, and
`cov_R_HBH`.

The subtracted term $S^i_{ob}S^j_w / c$ is the sample estimate of
$c \cdot \bar{d}_{ob} \cdot \bar{d}_w$, i.e. the bias in the un-centred
cross-product that arises when the mean innovations are non-zero.

### 2.2 No bias removed

$$
\hat{C}_\text{no bias}[v]_{s,k} = \frac{E_v}{c - 1}
$$

This is the un-centred cross-product divided by $(c-1)$. It equals the
bias-corrected estimate **plus** the mean-product term. It is provided so
users can assess the magnitude of the innovation bias contribution. Only
computed for `HBH` and `R+HBH`.

### 2.3 With innovation bias

The mean innovation for predictor bin $s$, averaged over all distance bins, is:

$$
\bar{d}_{ob}[s] = \frac{1}{n_\text{grid}}
\sum_k \frac{S^i_{ob}[s,k]}{c[s,k]}
$$

This mean is then added back as an isotropic bias covariance:

$$
\hat{C}_\text{with innov bias}[v]_{s,k} =
\hat{C}[v]_{s,k} + \frac{c}{c-1}\bar{d}_{ob}[s]^2
$$

This represents the covariance that would be estimated if the innovation bias
were treated as a signal rather than removed. Only computed for `HBH` and
`R+HBH`.

### 2.4 Bias term

The bias at zero separation ($k = 0$) for each innovation vector is:

$$
\text{bias}_v[s] = \frac{S^j_v[s, 0]}{c[s, 0]}
$$

This is the **mean** of the respective innovation/residual variable across
all observations that fall in predictor bin $s$ and distance bin $k=0$.
A non-zero value indicates a systematic offset in the assimilation system:
- `bias_d_oa` ≠ 0 → systematic observation-minus-analysis residual (observation bias)
- `bias_d_ab` ≠ 0 → systematic analysis increment (background bias)
- `bias_d_ob_j` ≠ 0 → systematic innovation (combined bias)

---

## 3. Standard deviation and correlation

### 3.1 Unsmoothed correlation

The zero-separation covariance ($k=0$) estimates the error **variance**. The
standard deviation and the raw (unsmoothed) normalised spatial correlation for
distance bin $k$ are:

$$
\sigma_v[s] = \sqrt{\hat{C}[v]_{s,0}}
$$

$$
\rho_v[s,k] = \min\left(\frac{\hat{C}[v]_{s,k}}{\sigma_v[s]^2}, 1\right)
$$

The $\min(\cdot, 1)$ clamp prevents numerical artefacts from producing
super-unity correlations. By construction $\rho_v[s, 0] = 1$.

### 3.2 Moving-average smoothing

The raw correlation $\rho_v$ is noisy at large separation distances where the
sample count $c_{s,k}$ is smaller. A centred moving-average filter is therefore
applied to all bins $k \ge 1$ before the result is saved to `cov.nc`.

For a smoothing half-width $h = \lfloor w/2 \rfloor$ (default $w = 6$,
$h = 3$), the smoothed correlation at bin $k$ is:

$$
\tilde{\rho}_v[s,k] =
\frac{1}{k_1 - k_0}
\sum_{k'=k_0}^{k_1 - 1} \rho_v[s,k'],
\qquad k \ge 1
$$

where the window bounds are clipped to the valid index range:

$$
k_0 = \max(0, k - h), \qquad k_1 = \min(n_\text{grid}, k + h)
$$

The zero-separation bin is always set to exactly 1:

$$
\tilde{\rho}_v[s, 0] = 1
$$

**Window width at the boundaries.** In the interior of the grid the window
contains $2h$ bins ($w = 6$ bins by default). Near the edges the window is
truncated:

| Distance bin $k$ | Window bins used | Effective width |
|-----------------|------------------|-----------------|
| 1 | $[0, 1, 2, 3]$ | 4 |
| 2 | $[0, 1, 2, 3, 4]$ | 5 |
| $3 \le k \le n-4$ | $[k-3, \ldots, k+2]$ | 6 (full) |
| $n-3$ | $[n-6, \ldots, n-2]$ | 5 |
| $n-2$ | $[n-5, \ldots, n-2]$ | 4 |
| $n-1$ | $[n-4, \ldots, n-2]$ | 3 |

Because the window uses integer half-width $h = w/2$, it contains $h$ bins to
the left of $k$ and $h-1$ bins to the right (i.e. it is slightly left-heavy
for even $w$). The smoothing window width $w$ can be changed by modifying the
`smooth_window` parameter in `CorrEstimator._convert_cov_to_corr`.

The smoothed values $\tilde{\rho}_v$ are what is written to the `corr_*`
variables in `cov.nc`.
