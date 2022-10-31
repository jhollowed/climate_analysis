import pdb
import sys
import pathlib

sys.path.append('/global/homes/j/jhollo/repos/CIME_namelist_automator')
from ensembler import ensembler

ic_dir = "/global/cscratch1/sd/jhollo/E3SM/E3SMv2_cases/hsw_cases/"\
         "E3SM_ne16_L72_FIDEAL_3year_ensICs/run"
# root case has power law heating with gamma=10
root_case = '/global/homes/j/jhollo/repos/climate_analysis/CLDERA/SAI/sai_case_builder_E3SM/current/cases/gamma_clones/HSW_SAI_ne16pg2_L72_allActive__cldera_sai_gammaq_10__cldera_sai_gammatau_10'
top_clone_dir = '/global/homes/j/jhollo/repos/climate_analysis/CLDERA/SAI/sai_case_builder_E3SM/current/cases/ic_ens'
top_output_dir='/global/cscratch1/sd/jhollo/E3SM/E3SMv2_cases/sai_cases/ic_ens'
cime_dir='/global/homes/j/jhollo/E3SM/CLDERA-E3SM_SAI/cime/scripts'
clone_prefix='HSW_SAI_ne16pg2_L72_allActive_gamma10'

ens = ensembler()
ens.add_members(ic_dir, globstr="*newCoordNames")

ens.create_members(root_case, top_clone_dir, top_output_dir, cime_dir, 
                   clone_prefix, overwrite=True, resubmits=4)
ens.submit_members(dry=False)

