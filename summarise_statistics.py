"""
Analyise the statistics for a fitting file
"""

import numpy as np
import numpy.ma as ma
import netCDF4 as nc4


class StatsObject:
    """
    Object to hold the statistics for a single variable
    """

    def __init__(self, name, data):
        self.name = name

        data_demasked = data.compressed()
        self.mean = np.mean(data_demasked)
        self.median = np.median(data_demasked)
        self.quantile = np.quantile(data_demasked, [0.25, 0.75])

        self.mean_iqr = np.mean(
            data_demasked[(data_demasked > self.quantile[0]) &
                          (data_demasked < self.quantile[1])])


class FitAnalysis:
    """
    A class for statistical analysis of a fitting file
    """
    ###############################
    # Initalise
    ###############################

    def __init__(self, file):
        """
        Initalise by opening file
        INPUTS:
            file: netcdf fitting file to open
        """

        self._fit_file = nc4.Dataset(file)
        self._num_params = self._get_num_params()

        self._stats_calcuated = False

    def _get_num_params(self):
        """
        Calculate the number of parameters in the fit
        """

        num = 0
        for key in self._fit_file.variables.keys():
            num += 1 if key[0:9] == "Magnitude" else 0
        return num

    ###################################################
    # Destroy
    ###################################################

    def __del__(self):
        self._fit_file.close()

    ###################################################
    # Output
    ###################################################

    def print_stats(self):
        """
        Print stats to screen
        """

        if not self._stats_calcuated:
            self._calc_stats()

        for n in range(self._num_params):
            print(f"\n\n\nStatistics for parameter {n}")
            print(f"\nMagnitude{n+1}")
            print(f"   Mean:     {self._mag_stats[n].mean}")
            print(f"   Median:   {self._mag_stats[n].median}")
            print(f"   IQR:      {self._mag_stats[n].quantile[0]}, "
                  + f"{self._mag_stats[n].quantile[1]}")
            print(f"   IQR mean: {self._mag_stats[n].mean_iqr}")
            print(f"\nScale{n+1}")
            print(f"   Mean:     {self._scale_stats[n].mean}")
            print(f"   Median:   {self._scale_stats[n].median}")
            print(f"   IQR:      {self._scale_stats[n].quantile[0]}, "
                  + f"{self._scale_stats[n].quantile[1]}")
            print(f"   IQR mean: {self._scale_stats[n].mean_iqr}")

    ###################################################
    # Calaculate stats
    ###################################################

    def _calc_stats(self):
        """Call routine to create stats objects"""

        self._mag_stats = []
        self._scale_stats = []
        for n in range(self._num_params):

            mag_name = f"Magnitude{n+1}"
            self._mag_stats += [StatsObject(
                mag_name, self._fit_file.variables[mag_name][:])]

            scale_name = f"Scale{n+1}"
            self._scale_stats += [StatsObject(
                scale_name, self._fit_file.variables[scale_name][:])]

        self._stats_calcuated = True


if __name__ == "__main__":

    fit_file = ("/data/users/frwe/ERROR_COVS_REPO/ocean_error_covs/immerse/"
                + "MSG_SEVIRI_fitting_stats.nc")
    stats = FitAnalysis(fit_file)

    stats.print_stats()

