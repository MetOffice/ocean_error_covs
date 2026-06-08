# Configuration Reference

The tool is driven by two levels of INI configuration files parsed with Python's
standard `configparser` library.

---

## Top-level `config.ini`

The top-level file controls which observation types to process, which
covariance matrices to estimate, and where to write the output.

```ini
[ObsTypes]
# Comma-separated list of observation config file names (stems or full paths).
# Each entry may be given with or without the .ini extension.
# Resolved first relative to the directory of config.ini, then relative to the
# current working directory.
obs_config_files = chlo, carbon

[Uncertainty]
# Comma-separated selection of covariance matrices to estimate.
# R       – observation error covariance
# HBH     – background error covariance in observation space
# R+HBH   – sum of both
estimate = R, HBH, R+HBH

# Set to False to switch from horizontal (geographic) distance binning to
# vertical / 1-D mode, where pairs are binned by the absolute difference of
# a single scalar coordinate (e.g. depth or pressure).
# Default: True
# When False, every observation-type .ini must supply vert_name instead of
# lon_name and lat_name, and the [Grid] distances are in the same units as
# that coordinate.
is_horizontal = True

[Output]
# Directory where output NetCDF files are written.
output_dir = output/all
# Stem of the intermediate binned-statistics file (no .nc extension).
binned_filename = binned
# Stem of the final covariance-estimates file (no .nc extension).
cov_filename = cov
```

---

## Per-observation-type `<obs>.ini`

Each observation type listed in `obs_config_files` has its own INI file with
four sections.

### `[General]`

```ini
[General]
# Short identifier used internally and in output variable names.
name = chlo
# Human-readable label written to output file attributes.
long_name = Phytoplankton Chlorophyll-a Concentration
```

### `[Input]`

```ini
[Input]
# Directory containing the input NetCDF innovation files.
input_dir = /path/to/innovation/files

# strftime-style filename pattern (% must be doubled because configparser
# treats % as an escape character).
# See https://docs.python.org/3/library/datetime.html#strftime-and-strptime-format-codes
filename_format = owt_chlo_%%Y%%m.nc

# Temporal frequency of the input files, using NumPy datetime units:
# e.g. '6h' (6-hourly), '1D' (daily), '1M' (monthly), '1Y' (yearly).
# Leave empty to load every NetCDF file found in input_dir.
data_freq = 1M

# Date range to process (NumPy datetime format).
data_start = 2015-02-01T00:00:00.00
data_end   = 2016-12-01T00:00:00.00

# Comma-separated list of variables to load from each file.
# Must include at least the innovation/residual variables and lon, lat.
# Include any additional fields needed by the predictor mask function
# (e.g. OWT_dom, OWT_tot for optical water types).
variables = d_ob, d_ab, d_oa, lon, lat

# Mapping from required internal names to variable names in the NetCDF files.
d_ob_name = d_ob   # observation minus background
d_ab_name = d_ab   # analysis increment in observation space
d_oa_name = d_oa   # observation minus analysis
lon_name  = lon    # required in horizontal mode (is_horizontal = True)
lat_name  = lat    # required in horizontal mode (is_horizontal = True)
# vert_name = depth  # required in vertical mode (is_horizontal = False); replaces lon_name/lat_name

# Variable holding the integer predictor bin index for each observation.
# Omit if no predictor is used (all observations treated as one bin),
# or if `mask_func` is used.
pred_name = mask
```

### `[Predictor]`

```ini
[Predictor]
# Name of the predictor (used in output variable names).
name = OWT

# Total number of predictor bins (bin indices run from 0 to n_classes-1).
n_classes = 14

# Maximum number of observations to load per time step.
# A random subset is drawn when the file exceeds this limit.
max_obs_per_day = 10000

# minimum number of samples to compute covariance estimates
n_samples = 100

# window size for smoothing correlations, better if even
smooth_window = 6
```

### `[Grid]`

Defines the discrete distance bins used to compute correlations.
In **horizontal mode** (`is_horizontal = True`) distances are great-circle
distances in **kilometres**. In **vertical mode** (`is_horizontal = False`)
distances are absolute differences of the `vert` coordinate, interpreted in
**whatever units `vert_name` uses** (e.g. metres for depth).

If the distance grid is `[0, 200, 400, ...]`, then the covariance at lag 200 is
estimated from all observation pairs whose separation falls in the interval
(0, 200], and so on.

```ini
[Grid]
# Upper bound of the distance grid (km in horizontal mode; vert units in vertical mode).
max_distance = 10000

# Width of each distance bin (same units as max_distance).
distance_interval = 200
```

---

## Vertical / 1-D mode

Setting `is_horizontal = False` in `[Uncertainty]` switches the tool from
horizontal geographic binning to **1-D distance binning** along a single scalar
coordinate (e.g. depth, pressure, or any along-track index).

### Changes to `[Input]`

Replace `lon_name` and `lat_name` with `vert_name`. The named variable must be
a 1-D array with one value per observation.

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
vert_name        = depth   # replaces lon_name / lat_name
```

### Changes to `[Grid]`

`max_distance` and `distance_interval` are now in the same units as `vert_name`.

```ini
[Grid]
max_distance      = 500   # metres (if depth is in metres)
distance_interval = 10
```

### Output

The `bins` coordinate in both output files (`binned.nc` and `cov.nc`) is in the
same units as `vert_name` rather than kilometres.

---

## Provided example configs

| File | Predictor | Notes |
|------|-----------|-------|
| `chlo_all.ini` / `carbon_all.ini` | No predictor (global pooling) | Pooled over all bins |
| `chlo_lat_month.ini` / `carbon_lat_month.ini` | Latitude band × month | Seasonal/latitudinal stratification |
