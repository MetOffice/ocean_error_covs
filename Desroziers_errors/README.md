# Desroziers Diagnostics for Ocean Data Assimilation

The functions here are used to estimate error variances and spatial correlations of observation
and background error in ocean data assimilation systems using
[Desroziers' statistics](https://doi.org/10.1256/qj.05.108). The implementation is based on

Fowler, A.M., Skákala, J. & Ford, D.(2023) Validating and improving the uncertainty assumptions for the assimilation of ocean-colour-derived chlorophyll a into a marine biogeochemistry model of the Northwest European Shelf Seas. Quarterly Journal of the Royal Meteorological Society, 149(750), 300–324. Available from: https://doi.org/10.1002/qj.4408

## Overview

Error covariance matrices **R** (observation errors) and **HBH**ᵀ (background
errors in observation space) are key inputs to any variational or ensemble data
assimilation system, yet they are rarely well known. This tool applies the
Desrozier's consistency diagnostics to estimate these matrices directly from the
innovation vectors produced by an assimilation system.

The diagnostics exploit three difference vectors available from the assimilation
cycle:

| Symbol | Meaning |
|--------|---------|
| **d_ob** = y − H(x_b) | observation minus background (innovation) |
| **d_oa** = y − H(x_a) | observation minus analysis (residual) |
| **d_ab** = H(x_a) − H(x_b) | analysis increment in observation space |

Under the assumption of linear statistical estimation theory:

$$
\mathbf{R} \approx \mathbb{E}[\mathbf{d_{oa}} \mathbf{d_{ob}}^\top], \qquad
\mathbf{HBH}^\top \approx \mathbb{E}[\mathbf{d_{ab}}, \mathbf{d_{ob}}^\top], \qquad
\mathbf{HBH^\top+R} \approx \mathbb{E}[\mathbf{d_{ob}} \mathbf{d_{ob}}^\top]
$$

Expectations are formed by pooling observations over time **and** over pairs of
observations that fall within the same spatial distance bin (up to a configurable
maximum distance). This yields spatially resolved error correlation functions as
well as error variances. See [statistics.md](statistics.md) for the full
derivation of the estimators.

### Predictor-dependent covariances

Error statistics can vary with oceanographic regime. This tool supports
stratifying the estimates by any discrete **predictor** variable, such as:

- **Optical Water Types (OWT)** - for ocean colour observations
- **Ocean bathymetry** - for depth-dependent error structures
- **Any integer-indexed categorical variable** present in the input files

The predictor can also be time or regions/space. For continuous predictors the
user bins them into integer categories before passing them to this tool. See
examples below.

The predictor can be read directly from a variable in the input data. The variable
has the same dimension as the observation vector, and each element of the predictor
variable should be the index of the predictor for the given observation.
Alternatively, a custom [`mask_func`](#custom-predictor-mask) can derive the bin
assignment from any combination of data variables at runtime.

### Multi-observation-type support

Multiple observation types can be processed simultaneously. All pairwise
cross-covariances between types are computed automatically, enabling the
construction of a full, blocked error covariance matrix.

---

## Input Data

The code reads **NetCDF files** containing observation-space vectors for each assimilation cycle. One set of input files is required per observation type; the file locations and variable names are specified in the per-observation-type `.ini` configuration file (see `configs/` for examples).

### Required variables

Each input file must contain **either** the horizontal coordinate pair **or** a single vertical / 1D coordinate:

| Variable | Description | Mode |
|----------|-------------|------|
| `lon` | Observation longitude [degrees] | horizontal (default) |
| `lat` | Observation latitude [degrees] | horizontal (default) |
| `vert` | Vertical / along-track coordinate (e.g. depth [m], pressure [hPa]) | vertical (`is_horizontal = False`) |

See [1D / vertical-coordinate mode](#1d--vertical-coordinate-mode) for setup details.

And at least **two** of the three Desroziers difference vectors — the third is derived automatically as the sum or difference of the other two:

| Variable | Meaning |
|----------|---------|
| `d_ob` = y − H(x_b) | Observation minus background (innovation) |
| `d_oa` = y − H(x_a) | Observation minus analysis (residual) |
| `d_ab` = H(x_a) − H(x_b) | Analysis increment in observation space |

The variable names in the NetCDF files do not need to match these defaults; they are mapped via the `d_ob_name`, `d_oa_name`, and `d_ab_name` keys in the observation `.ini` file. These variables
must be given in one-dimensional arrays with the size of the number of observations.

### Optional variables

A **predictor variable** can be included to stratify the statistics into discrete bins (e.g. Optical Water Type, bathymetry category). Its name in the NetCDF file is set with `pred_name` in the observation `.ini` file. Each element of the predictor array must be an integer bin index for the corresponding observation.

Any additional variables needed to build a custom predictor mask (see [Custom predictor mask](#custom-predictor-mask)) can be listed under `variables` in the `.ini` file and will be passed to `mask_func` at runtime.

### File naming and time coverage

Files are located using the `input_dir`, `filename_format`, `data_start`, `data_end`, and `data_freq` settings in the observation `.ini` file:

| Key | Description |
|-----|-------------|
| `input_dir` | Directory containing the input NetCDF files |
| `filename_format` | `strftime`-style format string used to generate file names from dates (e.g. `data_%Y%m%dT%H%MZ.nc`) |
| `data_freq` | Time step between files using [NumPy datetime units](https://numpy.org/doc/stable/reference/arrays.datetime.html#datetime-units) (e.g. `1D` for daily, `1M` for monthly). Leave blank to read all `.nc` files in `input_dir` |
| `data_start` | Start of the date range in ISO 8601 format (e.g. `2015-02-01T00:00:00.00`) |
| `data_end` | End of the date range (inclusive) in ISO 8601 format |

For example, the configuration below reads monthly files from February to March 2015:

```ini
[Input]
input_dir        = /path/to/data
filename_format  = obs_%%Y%%m.nc
data_freq        = 1M
data_start       = 2015-02-01T00:00:00.00
data_end         = 2015-03-01T00:00:00.00
variables        = d_ob, d_ab, d_oa, lon, lat
d_ob_name        = d_ob
d_ab_name        = d_ab
d_oa_name        = d_oa
lon_name         = lon
lat_name         = lat
```

> **Note:** In `configparser` `.ini` files the `%` character must be escaped as `%%` inside `filename_format`.

---

## 1D / vertical-coordinate mode

By default the tool works in **horizontal** mode: observation pairs are binned by
great-circle distance and the output describes spatially resolved error
correlation functions. To process **one-dimensional data** (e.g. vertical profiles
indexed by depth or pressure, or any along-track data with a single scalar
coordinate) set `is_horizontal = False` in the top-level `config.ini`:

```ini
[Uncertainty]
estimate        = R, HBH, R+HBH
is_horizontal   = False   # default is True
```

In vertical mode:

- Each observation-type `.ini` file must supply `vert_name` instead of `lon_name` / `lat_name`.
- Observation pairs are binned by the **absolute difference of their `vert` values** rather than great-circle distance.
- `max_distance` and `distance_interval` in `[Grid]` are interpreted in the **same units as `vert_name`** (e.g. metres if `vert_name` points to a depth variable in metres).
- The `bins` coordinate in the output files is likewise in those units.

### Example observation config (vertical mode)

```ini
[Input]
input_dir        = /path/to/profiles
filename_format  = profiles_%%Y%%m%%d.nc
data_freq        = 1D
data_start       = 2015-02-01T00:00:00.00
data_end         = 2015-03-01T00:00:00.00
variables        = d_ob, d_ab, d_oa, depth
d_ob_name        = d_ob
d_ab_name        = d_ab
d_oa_name        = d_oa
vert_name        = depth

[Grid]
max_distance      = 500    # metres
distance_interval = 10     # metres
```

---

## Usage

The tool is controlled by a top-level `config.ini` file and one `.ini` file per
observation type. See [configs/README.md](configs/README.md) for a full
reference before getting started.

### Python API

```python
from desroziers import calc_desroziers

calc_desroziers('my_config.ini')
```

### Custom predictor mask

For complex predictor binning logic (e.g. filtering by data quality flags before
selecting the dominant OWT), pass a `mask_func` callable:

```python
import numpy as np
from desroziers import calc_desroziers

def mask_func_owt(name, data, i_bin):
    total_owt_limits = [0.5, 1.5]
    mask = data['OWT_dom'].copy()
    mask[data['OWT_tot'] < total_owt_limits[0]] = -1
    mask[data['OWT_tot'] > total_owt_limits[1]] = -1
    return mask == i_bin

calc_desroziers('my_config.ini', mask_func=mask_func_owt)
```

The `mask_func` signature must be:

| Parameter | Type | Description |
|-----------|------|-------------|
| `name` | `str` | Observation type name (from `[General] name` in its `.ini` file) |
| `data` | `dict[str, np.ndarray]` | All variables loaded for the current time step |
| `i_bin` | `int` | Zero-based predictor bin index |

It must return a boolean `np.ndarray` of length `len(data['lon'])`.
The variables in data are specified in observation type configuration files.
It always contains at least ``'lon'`` and ``'lat'`` or ``'vert'``, and
``'d_ob_j'``. Depending on the required estimates it can also contain ``'d_oa'``,
and ``'d_ab'``.


#### Example: combined hemisphere and seasonal predictor

To have predictor by **hemisphere** (north / south) and **half-year** (Jan–Jun /
Jul–Dec), giving 4 bins in total, set `n_classes = 4` in the `[Predictor]`
section and include `month` in the `variables` list of the obs `.ini` file so
that the month index is available in `data`.

```ini
# in the observation .ini file
[Input]
variables = d_ob, d_ab, d_oa, lon, lat, month

[Predictor]
n_classes = 4
```

The bin layout is:

| `i_bin` | Hemisphere | Half-year |
|---------|------------|-----------|
| 0 | North (lat ≥ 0) | Jan – Jun |
| 1 | North (lat ≥ 0) | Jul – Dec |
| 2 | South (lat < 0) | Jan – Jun |
| 3 | South (lat < 0) | Jul – Dec |

```python
import numpy as np
from desroziers import calc_desroziers

def mask_func_region_time(name, data, i_bin):
    # 0 = Northern hemisphere, 1 = Southern hemisphere
    hemisphere = (data['lat'] < 0).astype(int)
    # 0 = first half of year (months 1-6), 1 = second half (months 7-12)
    half_year = (data['month'] > 6).astype(int)
    # encode into a single bin index: 0, 1, 2, 3
    combined = hemisphere * 2 + half_year
    return combined == i_bin

calc_desroziers('my_config.ini', mask_func=mask_func_region_time)
```

This pattern generalises to any number of discrete categorical predictors: encode
them into a single integer index `combined = predictor_a * n_b + predictor_b`
and set `n_classes = n_a * n_b`.

---

## Output

Two NetCDF files are written per observation-type pair. For a pair `obs_a` and
`obs_b` they are named `<binned_filename>_obs_a-obs_b.nc` and
`<cov_filename>_obs_a-obs_b.nc` (both names are set in `config.ini`). When only
one observation type is configured, both `obs_a` and `obs_b` are the same name.

### Dimensions (both files)

| Dimension | Description |
|-----------|-------------|
| `n_grid` | Number of spatial distance bins/grid |
| `n_predictors` | Number of predictor bins |

### `binned.nc` — accumulated statistics

Intermediate file storing the raw sums accumulated over all time steps.
Re-running with new data can be supported by loading and adding to these arrays.

| Variable | Shape | Description |
|----------|-------|-------------|
| `bins` | `(n_grid,)` | Right edge of each distance bin [km] |
| `total_obs_bin` | `(n_predictors, n_grid)` | Number of observation pairs accumulated in each predictor × distance bin |
| `d_ob_i` | `(n_predictors, n_grid)` | Sum of $y - \mathbf{H}(\mathbf{x}_b)$ at the *first* point of each observation pair |
| `d_ob_j` | `(n_predictors, n_grid)` | Sum of $y - \mathbf{H}(\mathbf{x}_b)$ at the *second* point of each observation pair |
| `d_oa` | `(n_predictors, n_grid)` | Sum of $y - \mathbf{H}(\mathbf{x}_a)$ at the second point |
| `d_ab` | `(n_predictors, n_grid)` | Sum of $\mathbf{H}(\mathbf{x}_a) - \mathbf{H}(\mathbf{x}_b)$ at the second point |
| `R` | `(n_predictors, n_grid)` | Sum of $(y-\mathbf{H}(\mathbf{x}_b)) \cdot (y-\mathbf{H}(\mathbf{x}_a))^\top$ — numerator for the **R** estimator |
| `HBH` | `(n_predictors, n_grid)` | Sum of $(y-\mathbf{H}(\mathbf{x}_b)) \cdot (\mathbf{H}(\mathbf{x}_a)-\mathbf{H}(\mathbf{x}_b))^\top$ — numerator for the **HBH** estimator |
| `R_HBH` | `(n_predictors, n_grid)` | Sum of $(y-\mathbf{H}(\mathbf{x}_b)) \cdot (y-\mathbf{H}(\mathbf{x}_b))^\top$ — numerator for the **R+HBH** estimator |

Only the variables corresponding to the `estimate` options selected in
`config.ini` are written.

### `cov.nc` — covariance estimates

Final file produced from the accumulated sums in `binned.nc`. Bins with fewer
than `ss_lim` (default 100) observation pairs are set to `NaN`. See
[docs/statistics.md](docs/statistics.md) for the full derivation of each
estimator.

**Coordinates**

| Variable | Shape | Description |
|----------|-------|-------------|
| `bins` | `(n_grid,)` | Right edge of each distance bin [km] |

**Bias-corrected estimates** (mean innovations subtracted — default)

| Variable | Shape | Description |
|----------|-------|-------------|
| `cov_R` | `(n_predictors, n_grid)` | observation error covariance as a function of separation distance |
| `cov_HBH` | `(n_predictors, n_grid)` | background error covariance in observation space |
| `cov_R_HBH` | `(n_predictors, n_grid)` | sum of observation and background error covariance  |
| `std_dev_R` | `(n_predictors,)` | observation error standard deviation (square root of zero-separation covariance) |
| `std_dev_HBH` | `(n_predictors,)` | background error standard deviation |
| `std_dev_R_HBH` | `(n_predictors,)` | sum of observation and background error standard deviation |
| `corr_R` | `(n_predictors, n_grid)` | observation error spatial correlation (smoothed, normalised by `std_dev_R`) |
| `corr_HBH` | `(n_predictors, n_grid)` | background error spatial correlation |
| `corr_R_HBH` | `(n_predictors, n_grid)` | sum of observation and background error spatial correlation |

**No bias removed** (un-centred cross-product; `HBH` and `R+HBH` only)

| Variable | Shape | Description |
|----------|-------|-------------|
| `cov_HBH_no_bias_removed` | `(n_predictors, n_grid)` | Expectation of the product of `d_ob`*`d_ab` |
| `cov_R_HBH_no_bias_removed` | `(n_predictors, n_grid)` | Expectation of the product of `d_ob`*`d_ob` |
| `std_dev_HBH_no_bias_removed` | `(n_predictors,)` | Standard deviation derived from `cov_HBH_no_bias_removed` |
| `std_dev_R_HBH_no_bias_removed` | `(n_predictors,)` | Standard deviation derived from `cov_R_HBH_no_bias_removed` |
| `corr_HBH_no_bias_removed` | `(n_predictors, n_grid)` | Correlation derived from `cov_HBH_no_bias_removed` |
| `corr_R_HBH_no_bias_removed` | `(n_predictors, n_grid)` | Correlation derived from `cov_R_HBH_no_bias_removed` |

**With innovation bias** (bias-corrected estimate plus squared mean innovation; `HBH` and `R+HBH` only)

| Variable | Shape | Description |
|----------|-------|-------------|
| `cov_HBH_with_innov_bias` | `(n_predictors, n_grid)` | `cov_HBH` with mean innovation as bias |
| `cov_R_HBH_with_innov_bias` | `(n_predictors, n_grid)` | `cov_R_HBH` with mean innovation as bias |
| `std_dev_HBH_with_innov_bias` | `(n_predictors,)` | Standard deviation derived from `cov_HBH_with_innov_bias` |
| `std_dev_R_HBH_with_innov_bias` | `(n_predictors,)` | Standard deviation derived from `cov_R_HBH_with_innov_bias` |
| `corr_HBH_with_innov_bias` | `(n_predictors, n_grid)` | Correlation derived from `cov_HBH_with_innov_bias` |
| `corr_R_HBH_with_innov_bias` | `(n_predictors, n_grid)` | Correlation derived from `cov_R_HBH_with_innov_bias` |

**Bias terms** (mean of the innovation/residual vector at zero separation)

| Variable | Shape | Description |
|----------|-------|-------------|
| `bias_d_oa` | `(n_predictors,)` | Mean $d_{oa} = y - \mathbf{H}(\mathbf{x}_a)$ per predictor bin; non-zero indicates observation bias |
| `bias_d_ab` | `(n_predictors,)` | Mean $d_{ab} = \mathbf{H}(\mathbf{x}_a) - \mathbf{H}(\mathbf{x}_b)$ per predictor bin; non-zero indicates background bias |
| `bias_d_ob_j` | `(n_predictors,)` | Mean $d_{ob} = y - \mathbf{H}(\mathbf{x}_b)$ per predictor bin; non-zero indicates combined bias |

Only variables corresponding to the `estimate` options selected in `config.ini`
are written. The `+` in `R+HBH` is replaced by `_` in all NetCDF variable names
(e.g. `cov_R_HBH`, `corr_R_HBH_with_innov_bias`).

---

## Contributions

Feel free to ask questions, or ask for features in issues and make pull request. You can also contact us by email:
Yumeng Chen (yumeng.chen@reading.ac.uk) or Alison Fowler (a.m.fowler@reading.ac.uk)
