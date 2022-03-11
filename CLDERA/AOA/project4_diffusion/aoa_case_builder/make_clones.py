import pdb
import sys
import pathlib
import numpy as np
from namelist_lattice import namelist_lattice

# =============================================================================
# =============================================================================

top_clone_dir = '/glade/u/home/jhollowed/CAM/cases_589/project4/clones'
top_output_dir = '/glade/scratch/jhollowed/CAM/cases_589/project4'
cime = '/glade/u/home/jhollowed/CAM/CAM_dev/cime/scripts'

if(0):
    # =============================== SE_NU ==========================================

    # define case and cime locations
    root_case = '/glade/u/home/jhollowed/CAM/cases_589/project4/cases/SE_ne16L72_mod3'
    stdout = '{}/automator_nu.log'.format(root_case)

    # construct lattice object
    lattice = namelist_lattice('cam', nofill=False)

    # add desired namelist variations
    # default, x10, x50
    lattice.expand('se_nu', values = [0.4e16, 0.4e17, 2e17])
    # default x 0.1
    lattice.expand('se_nu_div', values = [0.2e16])
    lattice.expand('RESUBMIT', values=[13], xmlchange=True)

    # create clones
    lattice.create_clones(root_case, top_clone_dir, top_output_dir, cime_dir=cime, stdout=stdout)

    # submit runs
    lattice.submit_clone_runs()

    # ================================= SE_NU_DIV ====================================

    root_case = '/glade/u/home/jhollowed/CAM/cases_589/project4/cases/SE_ne16L72_mod3'
    stdout = '{}/automator_nu_div.log'.format(root_case)
    lattice = namelist_lattice('cam', nofill=False)
    # default, x10, x50
    lattice.expand('se_nu_div', values=[2.0e16, 2.0e17, 1.0e18])
    # default x 0.1
    lattice.expand('se_nu', values = [0.4e15])
    lattice.expand('RESUBMIT', values=[13], xmlchange=True)
    lattice.create_clones(root_case, top_clone_dir, top_output_dir, cime_dir=cime, stdout=stdout)
    lattice.submit_clone_runs()


    # =================================== FV3 RF ENSEMBLE =============================

    root_case = '/glade/u/home/jhollowed/CAM/cases_589/project4/cases/FV3_C48L72_mod4'
    stdout = '{}/automator_nu_div.log'.format(root_case)
    lattice = namelist_lattice('cam', nofill=False)
    lattice.expand('pertlim', values = (1e-2 + (np.random.rand(4) - 0.5)/1000))
    lattice.expand('RESUBMIT', values=[13], xmlchange=True)
    lattice.create_clones(root_case, top_clone_dir, top_output_dir, cime_dir=cime, stdout=stdout,
                          clone_sfx=['ens1', 'ens2', 'ens3', 'ens4'])
    lattice.submit_clone_runs()

# =================================== FV3 CJ ENSEMBLE =============================

root_case = '/glade/u/home/jhollowed/CAM/cases_589/project4/cases/FV3_C48L72_mod3'
stdout = '{}/automator_nu_div.log'.format(root_case)
lattice = namelist_lattice('cam', nofill=False)
lattice.expand('pertlim', values = (1e-2 + (np.random.rand(4) - 0.5)/1000))
lattice.expand('RESUBMIT', values=[13], xmlchange=True)
lattice.create_clones(root_case, top_clone_dir, top_output_dir, cime_dir=cime, stdout=stdout, 
                      clone_sfx=['ens1', 'ens2', 'ens3', 'ens4'])
lattice.submit_clone_runs()


# =================================== FV3 CJ ENSEMBLE w/ less diffusive vort =============================

root_case = '/glade/u/home/jhollowed/CAM/cases_589/project4/cases/FV3_C48L72_mod3'
stdout = '{}/automator_nu_div.log'.format(root_case)
lattice = namelist_lattice('cam', nofill=False)
lattice.expand('pertlim', values = (1e-2 + (np.random.rand(4) - 0.5)/1000))
lattice.expand('fv3_hord_vt', values = ['8'])
lattice.expand('RESUBMIT', values=[13], xmlchange=True)
lattice.create_clones(root_case, top_clone_dir, top_output_dir, cime_dir=cime, stdout=stdout, 
                      clone_sfx=['ens1_vt8', 'ens2_vt8', 'ens3_vt8', 'ens4_vt8'], overwrite=True)
lattice.submit_clone_runs()

