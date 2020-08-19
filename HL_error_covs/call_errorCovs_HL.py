import glob
import master

input_dir = "/data/users/dcarneir/data/ncfiles/output_error_covs"
files = input_dir + '/' + 'test_accum_stats*.nc'
list_of_files = sorted(glob.glob(files))
outfilename = '/data/users/dcarneir/data/ncfiles/output_error_covs/test_errorcovs.nc'
master.HL_error_covs(list_of_files, outfilename=outfilename)

