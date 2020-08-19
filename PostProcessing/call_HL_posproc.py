import numpy as np
import netCDF4 as nc4
from scipy.interpolate import griddata
import time
import master

FILE='/data/users/dcarneir/data/ncfiles/output_error_covs/test3/errorcovs/sst_id34_errorcovs.nc'
OUTFIG='/data/users/dcarneir/data/ncfiles/output_error_covs/figures'
MIN_ROSS = 21.0

#big_lscale
dat=nc4.Dataset(FILE)
lons=dat.variables["longitude"][:]
lats=dat.variables["latitude"][:]
dat.close()
nlons=len(lons)
nlats=len(lats)
big_lscale = np.ones((nlats,nlons)) * 4. *  111.  #4 degrees in latitude
for pos,lat in zip(np.arange(0,nlats),lats):
    big_lscale[pos,:] += 4. * 111. * np.cos(np.deg2rad(lat))
big_lscale /= 2.

#small lscale
rossby=nc4.Dataset("/data/users/dcarneir/data/ncfiles/rossby_radii/rossby_radii.nc")
ross_rad = rossby.variables["xscale"][:,:].flatten() /1000. #convert to km
nav_lon = rossby.variables["lon"][:,:].flatten()
nav_lat = rossby.variables["lat"][:,:].flatten()

outlon,outlat = np.meshgrid(lons,lats)
ross_rad[ross_rad<MIN_ROSS] = MIN_ROSS   # needed to remove coasts before interp

ross_lscale = griddata((nav_lon,nav_lat), ross_rad, (outlon,outlat), method="nearest")

plot_locs = [(128,48),(126,132),(160,81)]

print('STARTING FITTING STEP')
start = time.time()
master.HL_fitting_function(FILE, "/data/users/dcarneir/data/ncfiles/output_error_covs/test1.nc", 
                          "MultiGauss", lenscale=(big_lscale,ross_lscale), plot=plot_locs, 
                           outfig="/data/users/dcarneir/data/ncfiles/output_error_covs/figures1",  
                           nproc=4, num_funcs=2, min_num_obs=20, max_iter=1000)
end = time.time()
print('TOTAL TIME IN MINUTES: ', (end-start)/60.0)
