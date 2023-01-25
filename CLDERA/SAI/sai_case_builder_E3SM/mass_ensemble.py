import pdb
import sys
import pathlib
import numpy as np

sys.path.append('/global/homes/j/jhollo/repos/cime_ensemble_automator')
from namelist_lattice import namelist_lattice

# root case is ens05 (which was named ens10 before renaming...) 
# (has symmetric IC in zonal-mean zonal wind)
root_case = '/global/homes/j/jhollo/repos/climate_analysis/CLDERA/SAI/sai_case_builder_E3SM/cases/ic_ens/HSW_SAI_ne16pg2_L72_900day_180delay__ens10'
top_clone_dir = '/global/homes/j/jhollo/repos/climate_analysis/CLDERA/SAI/sai_case_builder_E3SM/cases/mass_ens'
top_output_dir='/global/cscratch1/sd/jhollo/E3SM/E3SMv2_cases/sai_cases/mass_ens'
cime_dir='/global/homes/j/jhollo/E3SM/CLDERA-E3SM_SAI/cime/scripts'
clone_prefix='HSW_SAI_ne16pg2_L72_1200day_180delay_ens05'

# clones for 0.25, 0.5, 1.5, 2x SO2 loading
lattice = namelist_lattice('eam')

# for first round of cloning (to ensure max submits for debug=3)
#lattice.expand('cldera_sai_mso2', values=np.array([0.25, 0.5]) * 17.0)
#sfx = ['MSO2_0.25X', 'MSO2_0.5X']
# for second round of cloning (to ensure max submits for debug=3)
lattice.expand('cldera_sai_mso2', values=np.array([1.5, 2]) * 17.0)
sfx = ['MSO2_1.5X', 'MSO2_2.0X']

# clones will each have resubmits=3 for 1200 days total run time
lattice.create_clones(root_case, top_clone_dir, top_output_dir, cime_dir, clone_prefix, 
                      resubmits=3, overwrite=True, clone_sfx=sfx)
lattice.submit_clone_runs(dry=False)

