#!/bin/bash -fe

# Joe Hollowed
# 20230720
# This run script generates a single member of a CLDERA "limited variability" ensemble, which
# starts on June1st 1991, just before the eruption of Mt. Pinatubo. This run script was copied
# and modified from a script provided by Benj Wagman (hereafter BMW).
# Inline comments added by me begin with 'JH' 
# Original header comments from BW below
#
# -------------------------------------------------------------------------------
#
# Benjamin M. Wagman
# 20230309
# Limited variability ensemble member.
# All ensemble members are Hybrid start on 1991-06-01 from initial condition 1988-06-01 in 
# v2.LR.WCYCL20TR.0211.trc.pmcpu.ens1 + pertlim. 
# Including diagnostic PV. 
# Enable hybrid start with mismatched refdate and startdate by doing on all mpas restarts: 
# ncrename -v xtime,xtime.orig {restart file}
#
# -------------------------------------------------------------------------------


main() {

# For debugging, uncomment line below
#set -x

# --- Configuration flags ----

# ----- Machine and project
readonly MACHINE=pm-cpu
readonly PROJECT="m4014" # CLDERA NERSC account.

# ----- Simulation
readonly COMPSET="WCYCL20TR"
readonly RESOLUTION="ne30pg2_EC30to60E2r2"
readonly CASE_NAME="v2.LR.WCYCL20TR.pmcpu.limvar.ens3"

# ----- Run options
# JH: E3SM model start options:
#   -- "initial"  - use specificed initial conditions?
#   -- "continue" - continue from existing restart files
#   From the CESM docs at 
#   https://www2.cesm.ucar.edu/models/cesm1.1/cesm/doc/modelnl/env_run.html#run_start :
#   -- "branch"   - In a branch run, all components are initialized using a consistent
#                  set of restart files from a previous run (determined by the
#                  RUN_REFCASE and RUN_REFDATE variables in env_run.xml).  The case name
#                  is generally changed for a branch run, although it does not have to
#                  be. In a branch run, setting RUN_STARTDATE is ignored because the
#                  model components obtain the start date from their restart datasets.
#                  Therefore, the start date cannot be changed for a branch run. Branch runs are
#                  typically used when sensitivity or parameter studies are required, or
#                  when settings for history file output streams need to be modified
#                  while still maintaining bit-for-bit reproducibility. Under this
#                  scenario, the new case is able to produce an exact bit-for-bit restart
#                  in the same manner as a continuation run IF no source code or
#                  component namelist inputs are modified. All models use restart files
#                  to perform this type of run.  RUN_REFCASE and RUN_REFDATE are required
#                  for branch runs. To set up a branch run, locate the restart tar file or restart
#                  directory for RUN_REFCASE and RUN_REFDATE from a previous run, then
#                  place those files in the RUNDIR directory.
#   -- "hybrid"   - A hybrid run indicates that the model is initialized more like a
#                  startup, but uses initialization datasets FROM A PREVIOUS case.  This
#                  is somewhat analogous to a branch run with relaxed restart
#                  constraints.  A hybrid run allows users to bring together combinations
#                  of initial/restart files from a previous case (specified by
#                  RUN_REFCASE) at a given model output date (specified by
#                  RUN_REFDATE). Unlike a branch run, the starting date of a hybrid run
#                  (specified by RUN_STARTDATE) can be modified relative to the reference
#                  case. In a hybrid run, the model does not continue in a bit-for-bit
#                  fashion with respect to the reference case. The resulting climate,
#                  however, should be continuous provided that no model source code or
#                  namelists are changed in the hybrid run.  In a hybrid initialization,
#                  the ocean model does not start until the second ocean coupling
#                  (normally the second day), and the coupler does a cold start without
#                  a restart file.
readonly MODEL_START_TYPE="hybrid"  # 'initial', 'continue', 'branch', 'hybrid'
readonly START_DATE="1991-06-01" 

# ----- Additional options for 'branch' and 'hybrid'
readonly GET_REFCASE=TRUE
readonly RUN_REFDIR="/pscratch/sd/w/wagmanbe/E3SM_simulations/CLDERA/v2.LR.WCYCL20TR.0211.trc.pmcpu.ens1.limvar/archive/rest/1988-06-01-00000_mpas_time_renamed"
readonly RUN_REFCASE="v2.LR.WCYCL20TR.0211.trc.pmcpu.ens1.limvar"
readonly RUN_REFDATE="1988-06-01"

# ----- Set paths
# point to my branch which implements PV, PT tracers (not in CLDERA master as of 7/20/23)
readonly CODE_ROOT="/global/homes/j/jhollo/E3SM/CLDERA-E3SM_PV"
readonly CASE_ROOT="/pscratch/sd/j/jhollo/E3SM/E3SMv2_cases/historical_cases/${CASE_NAME}"

# ----- Sub-directories
readonly CASE_BUILD_DIR=${CASE_ROOT}/build
readonly CASE_ARCHIVE_DIR=${CASE_ROOT}/archive

# ----- Define type of run
#  short tests: 'XS_2x5_ndays', 'XS_1x10_ndays', 'S_1x10_ndays', 
#               'M_1x10_ndays', 'M2_1x10_ndays', 'M80_1x10_ndays', 'L_1x10_ndays'
#  or 'production' for full simulation
#readonly run='S_1x10_ndays'  
readonly run='production' 
if [ "${run}" != "production" ]; then

  # Short test simulations
  tmp=($(echo $run | tr "_" " "))
  layout=${tmp[0]}
  units=${tmp[2]}
  resubmit=$(( ${tmp[1]%%x*} -1 ))
  length=${tmp[1]##*x}

  readonly CASE_SCRIPTS_DIR=${CASE_ROOT}/tests/${run}/case_scripts
  readonly CASE_RUN_DIR=${CASE_ROOT}/tests/${run}/run
  readonly PELAYOUT=${layout}
  readonly WALLTIME="0:30:00"
  readonly STOP_OPTION=${units}
  readonly STOP_N=${length}
  readonly REST_OPTION=${STOP_OPTION}
  readonly REST_N=${STOP_N}
  readonly RESUBMIT=${resubmit}
  readonly DO_SHORT_TERM_ARCHIVING=false

else

  # ----- Production simulation
  # JH: These values modified for 4-month runs with PV,PT,E90,ST80 tracers enabled
  readonly CASE_SCRIPTS_DIR=${CASE_ROOT}/case_scripts
  readonly CASE_RUN_DIR=${CASE_ROOT}/run
  readonly PELAYOUT="custom-58"  # 23 SYPD according to NDK. BMW getting ~7.5 SYPD  
  readonly WALLTIME="1:00:00"    # 1-hour limit for regular queue on pm-cpu
                                 # JH: BMW ran for 11 hours for ~6.5 years, we only need 4 months
  readonly STOP_OPTION="nmonths"
  readonly STOP_N="4"            # 1991-06-01 to 1991-10-01
  readonly REST_OPTION="nmonths"
  readonly REST_N="1"            # How often to write a restart file
  readonly RESUBMIT="0"          # Submissions after initial one
  readonly DO_SHORT_TERM_ARCHIVING=false
fi

# ---- Coupler history
# sets coupler snapshot history file frequency 
readonly HIST_OPTION="nyears"
readonly HIST_N="5"

# ----- Leave empty (unless you understand what it does)
readonly OLD_EXECUTABLE=""

# ----- enable packages
# CLDERA_STRAT_VOLC: turn cldera stratospheric aerosol treatment on (1) or off (0). Off by default.
# CLDERA_PASSIVE_TRCRS: turn cldera e90, st80 tracers on (1) or off (0). Off by default.
# CLDERA_DYNAMIC_TRCRS: turn cldera pv, pt tracers on (1) or off (0). Off by default.
export CLDERA_STRAT_VOLC=1
export CLDERA_PASSIVE_TRCRS=1 # Adds E90j, ST80_25j, AOA
export CLDERA_DYNAMIC_TRCRS=1 # Adds PV_TRCR, PT_TRCR
echo 'cldera_strat_volc =' $CLDERA_STRAT_VOLC
echo 'cldera_passive_trcs =' $CLDERA_PASSIVE_TRCRS
echo 'cldera_dynamic_trcs =' $CLDERA_DYNAMIC_TRCRS

# ----- Toggle flags for what to do
do_create_newcase=true
do_case_setup=true
do_case_build=true
do_case_submit=true

# ----- call functions
#
# Make directories created by this script world-readable
umask 022
# Create case
create_newcase
# Custom PE layout
custom_pelayout
# Setup
case_setup
# Build
case_build
# Configure runtime options
runtime_options
# Copy script into case_script directory for provenance
copy_script
# Submit
case_submit
# All done
echo $'\n----- All done -----\n'

}

# =======================
# Custom user_nl settings
# =======================

user_nl() {


cat << EOF >> user_nl_eam

!!! ensemble generation utilities. For each enesmble, change seed_custom (i.e., seed_custom=1, seed_custom=2, ...) !!!
 pertlim = 1D-14
 new_random = .true.
 seed_clock = .false.
 seed_custom = 3

!! JH: Here we will only be outputting PV, PT, E90, ST80 quantities. This run should match the 
!!     limvar ensemble members which are already on disk, and so all other variables are already
!!     contianed there...
!!     T100 will be included in the monthly outputs for later verification of the statement above

!! JH: suppress h3, h4, h5; not needed
!!                 h0, h1,  h2
 nhtfrq =          0, -24, -24
 mfilt  =          1,  30,  30
 avgflag_pertape = 'A','A', 'I'

fincl1 = 'T100','E90j','ST80_25j','PV','PV_TRCR','PT','PT_TRCR'

fincl2 = 'E90j', 'ST80_25j','PV','PV_TRCR','PT','PT_TRCR'

fincl3 = 'PV','PV_TRCR','PT','PT_TRCR'


 ext_frc_specifier              = 'SO2         -> /global/cfs/cdirs/m4014/data/E3SM/model_input/emissions/volc/merge_cmip6_and_volc_so2/merged_cmip6_and_VolcanEESMv3.11_1850-2014_180x360_c20221019.nc',
         'SOAG        -> /global/cfs/cdirs/e3sm/inputdata/atm/cam/chem/trop_mozart_aero/emis/DECK_ne30/cmip6_mam4_soag_elev_1850-2014_c180205.nc',
         'bc_a4       -> /global/cfs/cdirs/e3sm/inputdata/atm/cam/chem/trop_mozart_aero/emis/DECK_ne30/cmip6_mam4_bc_a4_elev_1850-2014_c180205.nc',
         'num_a1      -> /global/cfs/cdirs/e3sm/inputdata/atm/cam/chem/trop_mozart_aero/emis/DECK_ne30/cmip6_mam4_num_a1_elev_1850-2014_c180205.nc',
         'num_a2      -> /global/cfs/cdirs/e3sm/inputdata/atm/cam/chem/trop_mozart_aero/emis/DECK_ne30/cmip6_mam4_num_a2_elev_1850-2014_c180205.nc',
         'num_a4      -> /global/cfs/cdirs/e3sm/inputdata/atm/cam/chem/trop_mozart_aero/emis/DECK_ne30/cmip6_mam4_num_a4_elev_1850-2014_c180205.nc',
         'pom_a4      -> /global/cfs/cdirs/e3sm/inputdata/atm/cam/chem/trop_mozart_aero/emis/DECK_ne30/cmip6_mam4_pom_a4_elev_1850-2014_c180205.nc',
         'so4_a1      -> /global/cfs/cdirs/e3sm/inputdata/atm/cam/chem/trop_mozart_aero/emis/DECK_ne30/cmip6_mam4_so4_a1_elev_1850-2014_c180205.nc',
         'so4_a2      -> /global/cfs/cdirs/e3sm/inputdata/atm/cam/chem/trop_mozart_aero/emis/DECK_ne30/cmip6_mam4_so4_a2_elev_1850-2014_c180205.nc'


 prescribed_volcaero_filetype ='VOLC_MIXING_RATIO'
 prescribed_volcaero_datapath ='/global/cfs/cdirs/m4014/data/E3SM/model_input/emissions/volc/zero-presc-volc-mmr'
 prescribed_volcaero_file='BMW_volcanic_1850-3009_all_zero.nc'


 use_hetfrz_classnuc  = .true.  !Default is true anyway
 hist_hetfrz_classnuc = .true.  !Default is false unless specified. 

  seasalt_emis_scale = 0.36
 dust_emis_fact = 3.255D0
  
 rad_climate            = 'A:Q:H2O', 'N:O2:O2', 'N:CO2:CO2', 'A:O3:O3', 'N:N2O:N2O', 'N:CH4:CH4', 'N:CFC11:CFC11', 'N:CFC12:CFC12', 
                           'M:mam4_mode1:/global/cfs/cdirs/m4014/data/E3SM/model_input/emissions/modal/modal_optics_mode1_hb.nc', 
                           'M:mam4_mode2:/global/cfs/cdirs/e3sm/inputdata/atm/cam/physprops/mam4_mode2_rrtmg_c130628.nc', 
                           'M:mam4_mode3:/global/cfs/cdirs/m4014/data/E3SM/model_input/emissions/modal/modal_optics_mode3_hb.nc', 
                           'M:mam4_mode4:/global/cfs/cdirs/e3sm/inputdata/atm/cam/physprops/mam4_mode4_rrtmg_c130628.nc'

rad_diag_1            = 'A:Q:H2O', 'N:O2:O2', 'N:CO2:CO2', 'A:O3:O3', 'N:N2O:N2O', 'N:CH4:CH4', 'N:CFC11:CFC11', 'N:CFC12:CFC12'

EOF

cat << EOF >> user_nl_elm

 ! Override
 check_finidat_fsurdat_consistency = .false.

 hist_nhtfrq =           0, -24, -6
 hist_mfilt =            1, 365, 120
 hist_avgflag_pertape = 'A','A', 'A'

 hist_fincl2 = 'EFLX_LH_TOT', 'FGR', 'FIRA', 'FIRE', 'FPSN', 'FSA', 'FSD24', 'FSH', 'FSI24', 'FSNO', 'FSR', 'H2OSNO', 'H2OSOI', 'QDRAI', 'QDRAI_XS', 'QH2OSFC', 'QIRRIG_REAL', 'QRUNOFF', 'QSNOMELT', 'QSOIL', 'QTOPSOIL', 'QVEGE', 'QVEGT', 'RAIN', 'SNOW', 'SNOW_DEPTH', 'SOILICE', 'SOILLIQ', 'SOILWATER_10CM', 'TG', 'TLAI', 'TSAI', 'TSOI', 'TV', 'TWS'

hist_fincl3 = 'EFLX_LH_TOT', 'FPSN', 'FSH', 'QH2OSFC', 'QSOIL', 'QVEGE', 'QVEGT', 'SOILWATER_10CM'

 
EOF

cat << EOF >> user_nl_mosart

EOF


}

# =====================================
# Customize MPAS stream files if needed
# =====================================

patch_mpas_streams() {

echo

}

# =====================================================
# Custom PE layout: custom-N where N is number of nodes
# =====================================================

custom_pelayout() {

if [[ ${PELAYOUT} == custom-58 ]];
then
    echo $'\n CUSTOMIZE 58 NODE CONFIGURATION for Perlmutter from NDK:'  # BMW 20221117


    pushd ${CASE_SCRIPTS_DIR}

   ./xmlchange NTHRDS=1
   ./xmlchange MAX_MPITASKS_PER_NODE=128

   ./xmlchange NTASKS_ATM=5504
   ./xmlchange NTASKS_CPL=688
   ./xmlchange NTASKS_OCN=1920
   ./xmlchange NTASKS_WAV=1
   ./xmlchange NTASKS_GLC=1
   ./xmlchange NTASKS_ICE=5248
   ./xmlchange NTASKS_ROF=256
   ./xmlchange NTASKS_LND=5248
   ./xmlchange NTASKS_ESP=1
   ./xmlchange NTASKS_IAC=1

   ./xmlchange NTHRDS_ATM=1
   ./xmlchange NTHRDS_CPL=1
   ./xmlchange NTHRDS_OCN=1
   ./xmlchange NTHRDS_WAV=1
   ./xmlchange NTHRDS_GLC=1
   ./xmlchange NTHRDS_ICE=1
   ./xmlchange NTHRDS_ROF=1
   ./xmlchange NTHRDS_LND=1
   ./xmlchange NTHRDS_ESP=1
   ./xmlchange NTHRDS_IAC=1

   ./xmlchange ROOTPE_ATM=0
   ./xmlchange ROOTPE_CPL=0
   ./xmlchange ROOTPE_OCN=5504
   ./xmlchange ROOTPE_WAV=0
   ./xmlchange ROOTPE_GLC=0
   ./xmlchange ROOTPE_ICE=0
   ./xmlchange ROOTPE_ROF=5248
   ./xmlchange ROOTPE_LND=0
   ./xmlchange ROOTPE_ESP=0
   ./xmlchange ROOTPE_IAC=0

   ./xmlchange PSTRID_CPL=8

    popd

fi

}

######################################################
### Most users won't need to change anything below ###
######################################################

create_newcase() {

    if [ "${do_create_newcase,,}" != "true" ]; then
        echo $'\n----- Skipping create_newcase -----\n'
        return
    fi

    echo $'\n----- Starting create_newcase -----\n'

    if [[ ${PELAYOUT} == custom-* ]];
    then
        layout="M" # temporary placeholder for create_newcase
    else
        layout=${PELAYOUT}

    fi
    ${CODE_ROOT}/cime/scripts/create_newcase \
        --case ${CASE_NAME} \
        --output-root ${CASE_ROOT} \
        --script-root ${CASE_SCRIPTS_DIR} \
        --handle-preexisting-dirs u \
        --compset ${COMPSET} \
        --res ${RESOLUTION} \
        --machine ${MACHINE} \
        --project ${PROJECT} \
        --walltime ${WALLTIME} \
        --pecount ${layout}

    if [ $? != 0 ]; then
      echo $'\nNote: if create_newcase failed because sub-directory already exists:'
      echo $'  * delete old case_script sub-directory'
      echo $'  * or set do_newcase=false\n'
      exit 35
    fi

}

#-----------------------------------------------------
case_setup() {

    if [ "${do_case_setup,,}" != "true" ]; then
        echo $'\n----- Skipping case_setup -----\n'
        return
    fi

    echo $'\n----- Starting case_setup -----\n'
    pushd ${CASE_SCRIPTS_DIR}

    # Setup some CIME directories
    ./xmlchange EXEROOT=${CASE_BUILD_DIR}
    ./xmlchange RUNDIR=${CASE_RUN_DIR}

    # Short term archiving
    ./xmlchange DOUT_S=${DO_SHORT_TERM_ARCHIVING^^}
    ./xmlchange DOUT_S_ROOT=${CASE_ARCHIVE_DIR}

    # Build with COSP, except for a data atmosphere (datm)
    if [ `./xmlquery --value COMP_ATM` == "datm"  ]; then 
      echo $'\nThe specified configuration uses a data atmosphere, so cannot activate COSP simulator\n'
    else
      echo $'\nConfiguring E3SM to use the COSP simulator\n'
      ./xmlchange --id CAM_CONFIG_OPTS --append --val='-cosp'
    fi

    # Extracts input_data_dir in case it is needed for user edits to the namelist later
    local input_data_dir=`./xmlquery DIN_LOC_ROOT --value`

    # wagmanbe
    # If cldera_chem_tracers, use the custom chem_mech.in file. 
    if [[ $CLDERA_CHEM_TRACERS -eq 1 ]]; then
	./xmlchange --id CAM_CONFIG_OPTS --append --val=' -cldera_chem_tracers'
	./xmlchange --id CAM_CONFIG_OPTS --append --val=' -usr_mech_infile $USR_MECH_INFILE'
    fi 

    # wagmanbe
    # If cldera_passive_tracers, pass the flag. 
    if [[ $CLDERA_PASSIVE_TRCRS -eq 1 ]]; then
	./xmlchange --id CAM_CONFIG_OPTS --append --val=' -cldera_passive_trcs'
    fi 
    # JH:
    # If cldera_dynamic_tracers, pass the flag. 
    if [[ $CLDERA_DYNAMIC_TRCRS -eq 1 ]]; then
	./xmlchange --id CAM_CONFIG_OPTS --append --val=' -cldera_dynamic_trcs'
    fi 

    #hybrown
    # Add CLDERA prog flag
    ./xmlchange --id CAM_CONFIG_OPTS --append --val=' -cldera_strat_volc $CLDERA_STRAT_VOLC'

    # Custom user_nl
    user_nl

    # Finally, run CIME case.setup
    ./case.setup --reset

    popd
}

#-----------------------------------------------------
case_build() {

    pushd ${CASE_SCRIPTS_DIR}

    # do_case_build = false
    if [ "${do_case_build,,}" != "true" ]; then

        echo $'\n----- case_build -----\n'

        if [ "${OLD_EXECUTABLE}" == "" ]; then
            # Ues previously built executable, make sure it exists
            if [ -x ${CASE_BUILD_DIR}/e3sm.exe ]; then
                echo 'Skipping build because $do_case_build = '${do_case_build}
            else
                echo 'ERROR: $do_case_build = '${do_case_build}' but no executable exists for this case.'
                exit 297
            fi
        else
            # If absolute pathname exists and is executable, reuse pre-exiting executable
            if [ -x ${OLD_EXECUTABLE} ]; then
                echo 'Using $OLD_EXECUTABLE = '${OLD_EXECUTABLE}
                cp -fp ${OLD_EXECUTABLE} ${CASE_BUILD_DIR}/
            else
                echo 'ERROR: $OLD_EXECUTABLE = '$OLD_EXECUTABLE' does not exist or is not an executable file.'
                exit 297
            fi
        fi
        echo 'WARNING: Setting BUILD_COMPLETE = TRUE.  This is a little risky, but trusting the user.'
        ./xmlchange BUILD_COMPLETE=TRUE

    # do_case_build = true
    else

        echo $'\n----- Starting case_build -----\n'

        # Turn on debug compilation option if requested
        if [ "${DEBUG_COMPILE^^}" == "TRUE" ]; then
            ./xmlchange DEBUG=${DEBUG_COMPILE^^}
        fi

        # Run CIME case.build
        ./case.build

        # Some user_nl settings won't be updated to *_in files under the run directory
        # Call preview_namelists to make sure *_in and user_nl files are consistent.
        ./preview_namelists

    fi

    popd
}

#-----------------------------------------------------
runtime_options() {

    echo $'\n----- Starting runtime_options -----\n'
    pushd ${CASE_SCRIPTS_DIR}

    # Set simulation start date
    ./xmlchange RUN_STARTDATE=${START_DATE}

    # Segment length
    ./xmlchange STOP_OPTION=${STOP_OPTION,,},STOP_N=${STOP_N}

    # Restart frequency
    ./xmlchange REST_OPTION=${REST_OPTION,,},REST_N=${REST_N}

    # Coupler history
    ./xmlchange HIST_OPTION=${HIST_OPTION,,},HIST_N=${HIST_N}

    # Coupler budgets (always on)
    ./xmlchange BUDGETS=TRUE

    # Set resubmissions
    if (( RESUBMIT > 0 )); then
        ./xmlchange RESUBMIT=${RESUBMIT}
    fi

    # Run type
    # Start from default of user-specified initial conditions
    if [ "${MODEL_START_TYPE,,}" == "initial" ]; then
        ./xmlchange RUN_TYPE="startup"
        ./xmlchange CONTINUE_RUN="FALSE"

    # Continue existing run
    elif [ "${MODEL_START_TYPE,,}" == "continue" ]; then
        ./xmlchange CONTINUE_RUN="TRUE"

    elif [ "${MODEL_START_TYPE,,}" == "branch" ] || [ "${MODEL_START_TYPE,,}" == "hybrid" ]; then
        ./xmlchange RUN_TYPE=${MODEL_START_TYPE,,}
        ./xmlchange GET_REFCASE=${GET_REFCASE}
	./xmlchange RUN_REFDIR=${RUN_REFDIR}
        ./xmlchange RUN_REFCASE=${RUN_REFCASE}
        ./xmlchange RUN_REFDATE=${RUN_REFDATE}
        echo 'Warning: $MODEL_START_TYPE = '${MODEL_START_TYPE} 
	echo '$RUN_REFDIR = '${RUN_REFDIR}
	echo '$RUN_REFCASE = '${RUN_REFCASE}
	echo '$RUN_REFDATE = '${START_DATE}
 
    else
        echo 'ERROR: $MODEL_START_TYPE = '${MODEL_START_TYPE}' is unrecognized. Exiting.'
        exit 380
    fi

    # Patch mpas streams files
    patch_mpas_streams

    popd
}

#-----------------------------------------------------
case_submit() {

    if [ "${do_case_submit,,}" != "true" ]; then
        echo $'\n----- Skipping case_submit -----\n'
        return
    fi

    echo $'\n----- Starting case_submit -----\n'
    pushd ${CASE_SCRIPTS_DIR}
    
    # Run CIME case.submit
    ./case.submit

    popd
}

#-----------------------------------------------------
copy_script() {

    echo $'\n----- Saving run script for provenance -----\n'

    local script_provenance_dir=${CASE_SCRIPTS_DIR}/run_script_provenance
    mkdir -p ${script_provenance_dir}
    local this_script_name=`basename $0`
    local script_provenance_name=${this_script_name}.`date +%Y%m%d-%H%M%S`
    cp -vp ${this_script_name} ${script_provenance_dir}/${script_provenance_name}

}

#-----------------------------------------------------
# Silent versions of popd and pushd
pushd() {
    command pushd "$@" > /dev/null
}
popd() {
    command popd "$@" > /dev/null
}

# Now, actually run the script
#-----------------------------------------------------
main
