import pdb
import sys
import glob
import pathlib

sys.path.append('/global/homes/j/jhollo/repos/cime_ensemble_automator')
from ensembler import ensembler

# for testing used in reporting results of variability meeting 03/09/23
# run starts at day 90 post-pertubation, with no injection delay; 
# guess for sufficient level of variability
#ic_dir = "/global/cscratch1/sd/jhollo/E3SM/E3SMv2_cases/sai_cases/pertlim_ens/pertlim_ics_day90"
#root_case = '/global/homes/j/jhollo/repos/climate_analysis/CLDERA/SAI/sai_case_builder_E3SM/cases/pertlim_ic_ens/HSW_SAI_ne16pg2_L72_300day_ptens'
#clone_prefix='HSW_SAI_ne16pg2_L72_300day'

# for production runs of tighter variability
# run starts at day 1 post-pertubation, 90-day injection delay;
# achieves same atm state at time of injection as testing configuration above, with lead time included
ic_dir = "/global/cscratch1/sd/jhollo/E3SM/E3SMv2_cases/sai_cases/pertlim_ens/pertlim_ics_day1"
root_case = '/global/homes/j/jhollo/repos/climate_analysis/CLDERA/SAI/sai_case_builder_E3SM/cases/pertlim_ic_ens/HSW_SAI_ne16pg2_L72_1200day_ptens'
clone_prefix='HSW_SAI_ne16pg2_L72_1200day_90delay'

top_clone_dir = '/global/homes/j/jhollo/repos/climate_analysis/CLDERA/SAI/sai_case_builder_E3SM/cases/pertlim_ic_ens'
top_output_dir='/global/cscratch1/sd/jhollo/E3SM/E3SMv2_cases/sai_cases/pertlim_ic_ens'
cime_dir='/global/homes/j/jhollo/E3SM/CLDERA-E3SM_SAI/cime/scripts'

ens = ensembler()
ens.add_members(ic_dir)
# uncomment if any other namelist settings need to be modified along with IC specification
#ens.lattice.expand('STOP_N', values=[180], xmlchange=True)

ens.create_members(root_case, top_clone_dir, top_output_dir, cime_dir, 
                   clone_prefix, overwrite=True, resubmits=0)
ens.submit_members(dry=False)
