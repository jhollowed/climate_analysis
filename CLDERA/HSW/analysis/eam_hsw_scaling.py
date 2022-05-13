from matplotlib import pyplot as plt
import subprocess
import numpy as np
import glob
import pdb

# -------------------------------------------------------

def gather_run_times(topDir):
    '''
    Return total and eam run times for all case directories located in topDir.
    It is assumed that each of these cases has a timing/ directory

    Parameters
    ----------
    topDir : str
        location of top directory containing E3SM CIME cases
    
    Return
    ------
    times : (N, 3) array
        (nproc, total run time , eam run time) for each case in topDir, where 
        the run times are in seconds
    '''

    cases = glob.glob('{}/*'.format(topDir))
    timings = []

    for case in cases:

        # get pecount
        caseName = case.split('{}/'.format(topDir))[-1]
        caseName.split()
        pe = np.array(subprocess.check_output(['cd {0}; {0}/xmlquery NTASKS'.format(case)], 
                                               shell=True).decode('utf-8').strip().split())
        pe = int(pe[[('ATM' in p) for p in pe]][0].replace("'","").replace(",","").split(':')[-1])
        
        # get timing
        try:
            timingFile = glob.glob('{}/timing/e3sm_timing.*'.format(case))[0]
        except IndexError:
            continue

        with open(timingFile) as tf:
            tfl = tf.readlines()
            tfl = np.array([line.rstrip() for line in tfl])
        totl_mask = [('TOT Run Time' in l) for l in tfl]
        atml_mask = [('ATM Run Time' in l) for l in tfl]
        totl = tfl[totl_mask][0]
        atml = tfl[atml_mask][0]

        # get wall seconds per model day, scale to wall minutes per model year
        tot_spmd = float(totl.split()[5])
        atm_spmd = float(atml.split()[5])
        tot_mpmy = (tot_spmd * 365) / 60
        atm_mpmy = (atm_spmd * 365) / 60

        timings.append((pe, tot_mpmy, atm_mpmy))

    timings = np.array(timings)[np.argsort([t[0] for t in timings])]
    return timings

# -------------------------------------------------------

def plot_scaling(times, title):
    '''
    Plot model scaling

    Parameters
    ----------
    times : (N, 3) array
        output of gather_run_times() above
    title : str
        title of model to label in plot title
    '''

    pe = [time[0] for time in times]
    tot = [time[1] for time in times]
    atm = [time[2] for time in times]
    perfect = np.hstack([tot[0], (tot[0] / (2*np.arange(len(pe))))[1:]])
    
    fig = plt.figure()
    ax = fig.add_subplot(111)

    ax.plot(pe, tot, '-^', color='b', label='total run time')
    ax.plot(pe, atm, '-^', color='r', label='EAM run time')
    ax.plot(pe, perfect, '-k', lw=0.75, label='perfect scaling')
    ylim = ax.get_ylim()
    ax.plot([1536, 1536], ylim, '--k', label='number of \ncubedsphere grid cells')
    ax.set_ylim(ylim)
    ax.set_xscale('log')
    ax.set_yticks(np.arange(2, 20, 2))
    ax.set_xticks(pe)
    ax.set_xticklabels(['%d' % p for p in pe])
    ax.grid()
    ax.minorticks_off()
    ax.set_xlabel('number of processes')
    ax.set_ylabel('wall minutes per simulated year')
    ax.legend()
    ax.set_title('E3SM scaling for {}'.format(title))
    plt.show()

    


        

if __name__ == '__main__':
    topDir = '/global/homes/j/jhollo/repos/climate_analysis/CLDERA/HSW/hsw_case_builder/cases'
    times = gather_run_times(topDir)
    plot_scaling(times, 'HSW, ne16np4 L72')
    pdb.set_trace()
        

    
