import pdb
import sys
import pathlib
import numpy as np

sys.path.append('/global/homes/j/jhollo/repos/cime_ensemble_automator')
from namelist_lattice import namelist_lattice

# root case is the counterfactual from the 011423 ensemble 
root_case = '/global/homes/j/jhollo/repos/climate_analysis/CLDERA/SAI/sai_case_builder_E3SM/cases/pertlim_ens/HSW_SAI_ne16pg2_L72_300day_pertT'
top_clone_dir = '/global/homes/j/jhollo/repos/climate_analysis/CLDERA/SAI/sai_case_builder_E3SM/cases/pertlim_ens'
top_output_dir='/global/cscratch1/sd/jhollo/E3SM/E3SMv2_cases/sai_cases/pertlim_ens'
cime_dir='/global/homes/j/jhollo/E3SM/CLDERA-E3SM_SAI/cime/scripts'
clone_prefix='HSW_ne16pg2_L72_300day'

# clones for 0.25, 0.5, 1.5, 2x SO2 loading
lattice = namelist_lattice('eam')

# for first round of cloning (to ensure max submits for debug=3)
#lattice.expand('pertlim', values=np.array([1e-4+2e-9, 1e-4+1e-9, 1e-4-1e-9, 1e-4-2e-9]))
#sfx = ['pert01', 'pert02', 'pert03', 'pert04']
# for second round of cloning (to ensure max submits for debug=3)
#lattice.expand('pertlim', values=np.array([1e-4+2.5e-9, 1e-4+1.5e-9, 1e-4-1.5e-9, 1e-4-2.5e-9]))
#sfx = ['pert05', 'pert06', 'pert07', 'pert08']
# for third round of cloning (to ensure max submits for debug=3)
lattice.expand('pertlim', values=np.array([1e-4+1.7e-9, 1e-4-1.7e-9]))
sfx = ['pert09', 'pert10']

# clones will each have resubmits=3 for 1200 days total run time
lattice.create_clones(root_case, top_clone_dir, top_output_dir, cime_dir, clone_prefix, 
                      resubmits=0, overwrite=True, clone_sfx=sfx)
lattice.submit_clone_runs(dry=False)
