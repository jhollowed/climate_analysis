import pdb
import sys
import pathlib

sys.path.append('/global/homes/j/jhollo/repos/cime_ensemble_automator')
from ensembler import ensembler

# model dir
cime_dir='/global/homes/j/jhollo/E3SM/CLDERA-E3SM_PV/cime/scripts'

#res = 'ne30'
res = 'ne16'

if(res == 'ne16'):
    # location of ICs for initializing ensemble members
    ic_dir = "/pscratch/sd/j/jhollo/E3SM/E3SMv2_cases/hsw_cases/E3SM_ne16_L72_FIDEAL_3year_ensICs/run"
    # root case to clone
    root_case = '/global/u2/j/jhollo/repos/climate_analysis/CLDERA/PV/cases/HSW_PVPT_ne16pg2_L72_30day'
    # clone dirs
    top_clone_dir = '/global/u2/j/jhollo/repos/climate_analysis/CLDERA/PV/cases/highvar_PVPTtest_ens'
    top_output_dir='/pscratch/sd/j/jhollo/E3SM/E3SMv2_cases/pv_cases/highvar_PVPTtest_ens'
elif(res == 'ne30'):
    ic_dir = "/pscratch/sd/j/jhollo/E3SM/E3SMv2_cases/hsw_cases/E3SM_ne30_L72_FIDEAL_3year_ensICs/run"
    root_case = '/global/u2/j/jhollo/repos/climate_analysis/CLDERA/PV/cases/HSW_PVPT_ne30pg2_L72_30day'
    top_clone_dir = '/global/u2/j/jhollo/repos/climate_analysis/CLDERA/PV/cases/highvar_PVPTtest_ne30_ens'
    top_output_dir='/pscratch/sd/j/jhollo/E3SM/E3SMv2_cases/pv_cases/highvar_PVPTtest_ne30_ens'
else:
    raise RuntimeError('unsupported res: {}'.format(res))

# --------------------------- Injection members --------------------------

# ---- create ensembler, add members
ens = ensembler()
ens.add_members(ic_dir, globstr="*newCoordNames")
# ---- create ensemble member cases, submit jobs
ens.create_members(root_case, top_clone_dir, top_output_dir, cime_dir, 
                   clone_prefix='injection', overwrite=True, resubmits=0)
ens.submit_members(dry=False)

# --------------------------- Counterfactual members --------------------------

# ---- create ensembler, add members
ens = ensembler()
ens.add_members(ic_dir, globstr="*newCoordNames")

# ---- add dimension to lattice to produce counterfactuals
ens.lattice.expand('cldera_sai_mso2,cldera_sai_mash', values=['0,0'], 
                    group=True, group_labels='Minject')

# ---- create ensemble member cases, submit jobs
ens.create_members(root_case, top_clone_dir, top_output_dir, cime_dir, 
                   clone_prefix='counterfactual', overwrite=True, resubmits=0)
ens.submit_members(dry=False)
