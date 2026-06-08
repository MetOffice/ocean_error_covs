"""Grid class for handling spatial grid information.

Author: Y Chen, University of Reading, 2025
"""
import configparser

import numpy as np

class Grid:
    """
    Class to hold grid information."""
    def __init__(self, config: configparser.SectionProxy, is_horizontal: bool) -> None:
        self.lon:dict[str, np.ndarray] = {}
        self.lat:dict[str, np.ndarray] = {}
        self.vert:dict[str, np.ndarray] = {}
        self.is_horizontal = is_horizontal
        # The distance, d, defines the resolution of spatial correlations.
        # For example, if bin = [0, 50, 100], the resulting covariance
        # will calculate the spatial covariance for distances 0-50km and 50-100km.
        # Observations with less than 50 km will be used for variance, and
        # observations with distance between 50 and 100 km will be used for
        # covariance of 50 km.
        max_distance = config.getfloat('max_distance', 350)
        distance_interval = config.getfloat('distance_interval', 5)
        self.d   = np.arange(0, max_distance, distance_interval)
        # number of grid points
        self.n_gp = len(self.d)

    def update_coord(self, data: dict[str, dict[str, np.ndarray]]) -> None:
        """Get coordinate arrays.

        Returns
        -------
        lon : dict[str, np.ndarray]
            Dictionary of longitude arrays for each observation type.
        lat : dict[str, np.ndarray]
            Dictionary of latitude arrays for each observation type.
        vert : dict[str, np.ndarray]
            Dictionary of vertical coordinate arrays for each observation type.
        """
        if self.is_horizontal:
            self.update_latlon(data)
        else:
            self.update_vert(data)

    def update_vert(self, data: dict[str, dict[str, np.ndarray]]) -> None:
        """Get vertical coordinate array.

        Returns
        -------
        vert : np.ndarray
            1D array of vertical coordinates.
        """
        for name in data:
            self.vert[name] = data[name]['vert']


    def update_latlon(self, data: dict[str, dict[str, np.ndarray]]) -> None:
        """Get latitude and longitude arrays.

        Returns
        -------
        lon : np.ndarray
            1D array of longitudes.
        lat : np.ndarray
            1D array of latitudes.
        """
        for name in data:
            self.lon[name] = data[name]['lon']
            self.lat[name] = data[name]['lat']


    def calc_distance(self, samples_a:dict, samples_b: dict) -> np.ndarray:
        """Calculate distances between observations in samples_a and samples_b.

        If samples_b is None, calculate distances between all pairs in samples_a.

        Returns
        -------
        distances : np.ndarray
            2D array of distances in km between observation pairs.
        """
        if self.is_horizontal:
            return self.calc_gcd(samples_a, samples_b)
        else:
            return self.calc_vert_distance(samples_a, samples_b)

    def calc_vert_distance(self, samples_a:dict, samples_b: dict) -> np.ndarray:
        """Calculate distances in km using haversine formula, and
        compute the closest distance bin for each pair.

        Uses chunking to reduce memory footprint for large sample sizes.

        Parameters
        ----------
        samples_a : dict
            Dictionary of observation indices for the first set of observations
        samples_b : dict
            Dictionary of observation indices for the second set of observations

        chunk_size : int, optional
            Number of rows to process per chunk (default: 500).
            Reduce for very large datasets to save memory; increase for speed.

        Returns
        -------
        bin_idx : np.ndarray
            The index of distance bin for each pair of observations
        """
        r = 6371.0  # Earth radius in km
        name_a, = samples_a.keys()
        name_b, = samples_b.keys()
        idx_a = samples_a[name_a]
        idx_b = samples_b[name_b]
        x_a = self.vert[name_a][idx_a]
        x_b = self.vert[name_b][idx_b]
        dx = x_a[:, None] - x_b[None, :]
        # Find nearest bin efficiently using searchsorted + check both neighbors
        bin_chunk = self._find_nearest_bins_fast(dx.ravel())

        return bin_chunk

    def calc_gcd(self, samples_a:dict, samples_b: dict) -> np.ndarray:
        """Calculate great circle distances in km using haversine formula, and
        compute the closest distance bin for each pair.

        Uses chunking to reduce memory footprint for large sample sizes.

        Parameters
        ----------
        samples_a : dict
            Dictionary of observation indices for the first set of observations
        samples_b : dict
            Dictionary of observation indices for the second set of observations

        chunk_size : int, optional
            Number of rows to process per chunk (default: 500).
            Reduce for very large datasets to save memory; increase for speed.

        Returns
        -------
        bin_idx : np.ndarray
            The index of distance bin for each pair of observations
        """
        r = 6371.0  # Earth radius in km
        name_a, = samples_a.keys()
        name_b, = samples_b.keys()
        idx_a = samples_a[name_a]
        idx_b = samples_b[name_b]
        lat_a = np.radians(self.lat[name_a][idx_a])
        lon_a = np.radians(self.lon[name_a][idx_a])
        lat_b = np.radians(self.lat[name_b][idx_b])
        lon_b = np.radians(self.lon[name_b][idx_b])

        # Pre-compute cosines for efficiency
        cos_lat_a = np.cos(lat_a)
        cos_lat_b = np.cos(lat_b)

        dlat = lat_a[:, None] - lat_b[None, :]
        dlon = lon_a[:, None] - lon_b[None, :]

        sin_dlat = np.sin(0.5 * dlat)
        sin_dlon = np.sin(0.5 * dlon)

        a = sin_dlat * sin_dlat + \
            cos_lat_a[:, None] * cos_lat_b[None, :] * sin_dlon * sin_dlon
        a = np.clip(a, 0.0, 1.0)
        dist = r*2.0*np.arctan2(np.sqrt(a), np.sqrt(1.0 - a))
        # Find nearest bin efficiently using searchsorted + check both neighbors
        bin_chunk = self._find_nearest_bins_fast(dist.ravel())

        return bin_chunk

    def _find_nearest_bins_fast(self, distances: np.ndarray) -> np.ndarray:
        """Find nearest bin for each distance using binary search.

        More efficient than creating a full (n_distances, n_bins) array.

        Parameters
        ----------
        distances : np.ndarray
            Flattened array of pairwise distances

        Returns
        -------
        bin_idx : np.ndarray
            Index of nearest bin for each distance
        """
        # Use searchsorted to find insertion points
        idx_right = np.searchsorted(self.d, distances)
        idx_left = np.clip(idx_right - 1, 0, len(self.d) - 1)
        idx_right = np.clip(idx_right, 0, len(self.d) - 1)

        # Choose left or right bin based on which is closer
        dist_left = np.abs(self.d[idx_left] - distances)
        dist_right = np.abs(self.d[idx_right] - distances)
        bin_idx = np.where(dist_left <= dist_right, idx_left, idx_right)

        return bin_idx

if __name__ == "__main__":
    pass  # For module testing purposes
