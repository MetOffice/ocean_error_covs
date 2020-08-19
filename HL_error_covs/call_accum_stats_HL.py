import glob
from time import time
import numpy as np
import master

#input_dir = "/data/users/dcarneir/data/ncfiles/test4"
#file_prefix = "sst_BiasCorrfb_oo_fdbk"
#obs_type = "SST"
#source_types = [34]
#ilat = 40
#elat = 65
#dlat = 1.0
#ilon = -20
#elon = 13
#dlon = 1.0
#bin_list = np.array([10, 15, 20, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500, 600, 700, 800, 900, 1000])
#profile = False
#depths = []

input_dir = "/data/users/dcarneir/data/ncfiles/test4"
file_prefix = "profb"
obs_type = "POTM"
source_types = []
ilat = -90
elat = 90
dlat = 5.0
ilon = -180
elon = 180
dlon = 5.0
bin_list = np.array([10, 15, 20, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500, 600, 700, 800, 900, 1000])
depths = np.array([50, 150, 300, 500, 750, 1000, 1500, 2000])
outfilename="/data/users/dcarneir/data/ncfiles/output_error_covs/test_accum_stats_03.nc"
nproc=4
files = input_dir + '/' + '2017*' +  file_prefix + '*'
list_of_files = sorted(glob.glob(files))

start = time()
master.HL_cov_accum_stats(list_of_files, obs_type=obs_type, source_types=source_types,
                 grid_def=[[ilat,elat,dlat],[ilon,elon,dlon]], bins=bin_list,
                 depth_boundaries=depths, outfilename=outfilename, nproc=nproc)
end = time()
print("Time to calc covs in minutes: {}".format((end - start)/60.0))

