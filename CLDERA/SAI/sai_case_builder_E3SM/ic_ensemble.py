import pdb
import sys
import pathlib

sys.path.append('/global/homes/j/jhollo/repos/cime_ensemble_automator')
from ensembler import ensembler

ic_dir = "/global/cscratch1/sd/jhollo/E3SM/E3SMv2_cases/hsw_cases/"\
         "E3SM_ne16_L72_FIDEAL_3year_ensICs/run"
# root case after tuning
root_case = '/global/homes/j/jhollo/repos/climate_analysis/CLDERA/SAI/sai_case_builder_E3SM/cases/HSW_SAI_ne16pg2_L72_900day_180delay'
top_clone_dir = '/global/homes/j/jhollo/repos/climate_analysis/CLDERA/SAI/sai_case_builder_E3SM/cases/ic_ens'
top_output_dir='/global/cscratch1/sd/jhollo/E3SM/E3SMv2_cases/sai_cases/ic_ens'
cime_dir='/global/homes/j/jhollo/E3SM/CLDERA-E3SM_SAI/cime/scripts'
clone_prefix='HSW_SAI_ne16pg2_L72_900day_180delay'

ens = ensembler()
ens.add_members(ic_dir, globstr="*newCoordNames")

ens.create_members(root_case, top_clone_dir, top_output_dir, cime_dir, 
                   clone_prefix, overwrite=True, resubmits=3)
ens.submit_members(dry=False)

