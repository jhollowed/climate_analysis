import pdb
import sys
import pathlib

sys.path.append('/global/homes/j/jhollo/repos/CIME_namelist_automator')
from namelist_lattice import namelist_lattice

lattice = namelist_lattice('eam', nofill=True)
lattice.expand('cldera_sai_MSO2', values=[17*2, 17*5, 17*10])

root_case = '/global/homes/j/jhollo/repos/climate_analysis/CLDERA/SAI/sai_case_builder_E3SM/current/cases/HSW_SAI_ne16pg2_L72_heatingEns_e90st80'
top_clone_dir = '/global/homes/j/jhollo/repos/climate_analysis/CLDERA/SAI/sai_case_builder_E3SM/current/cases/HSW_SAI_ne16pg2_L72_heatingEns_e90st80_CLONES'
top_output_dir='/global/cscratch1/sd/jhollo/E3SM/E3SMv2_cases/sai_cases/HSW_SAI_ne16pg2_L72_heatingEns_e90st80_CLONES'
cime_dir='/global/homes/j/jhollo/E3SM/CLDERA-E3SM_SAIBranch/cime/scripts'
clone_prefix=None

lattice.create_clones(root_case, top_clone_dir, top_output_dir, cime_dir=cime_dir, overwrite=True)
lattice.submit_clone_runs()

