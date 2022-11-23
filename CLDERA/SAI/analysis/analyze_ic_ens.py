import pdb
import glob
from sai_analysis import sai_analyzer

mean_clim = '/global/cscratch1/sd/jhollo/E3SM/E3SMv2_cases/hsw_cases/E3SM_ne16_L72_FIDEAL_10year_spinup/run/E3SM_ne16_L72_FIDEAL_10year_spinup.eam.h0.0001-01-01-00000.nc'
slice_range = (5*360, 10*360)

members = glob.glob('/global/cscratch1/sd/jhollo/E3SM/E3SMv2_cases/sai_cases/ic_ens/*ens0[2-6]/run/*eam.h2*')
#members.append(glob.glob('/global/cscratch1/sd/jhollo/E3SM/E3SMv2_cases/sai_cases/ic_ens/*ens09/run/*eam/h0*')[0])

dat_savedest='/global/cscratch1/sd/jhollo/E3SM/E3SMv2_cases/sai_cases/tmp_analysis'

A = sai_analyzer(data_savedest=dat_savedest)
A.set_mean_climate(mean_clim, slice_range=slice_range)
A.set_ensemble_members(members)
A.compute_ensemble_mean()
A.compute_anomaly()

pdb.set_trace()
