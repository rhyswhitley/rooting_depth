#!/usr/bin/env python3

import pandas as pd
import os, re
import subprocess as sub
#from joblib import Parallel, delayed
from multiprocessing import cpu_count

def run_SPAbin(dirpath):
    # Echo to user
    print("\t> executing simulation")

    # take the head of the file path to tell the subproces which directory
    # to execute in
#    bin_cwd = os.path.split(bin_fpath)[0]
#
#    # execture the SPA Fortran binary
#    proc = sub.Popen(bin_fpath, cwd=bin_cwd, bufsize=1, shell=False,
#                        close_fds=False)

    bin_fpath = dirpath + "/SPA"
    if not os.path.exists(bin_fpath):
        print("\t\tSimlink missing! Creating one for you now. =)")
        os.symlink(BINLOC, bin_fpath)

#    proc = sub.Popen(bin_fpath, cwd=dirpath, bufsize=1, shell=False,
#                        close_fds=False)
#
#    # Slight pause to allow python to catchup with the terminal
#    proc.wait()

def create_inf_file(input_inf, soil, odir, slen, spinup=1, inc_spin=0):
    """
    Creates the temporary configuration file for the experiment. This is
    overwritten when the following simulation is executed.
    """

    mod_odir = os.path.basename(odir)

    file_list = input_inf + [soil, mod_odir, slen-1, spinup, inc_spin]

    index_l = ['meto_file', 'phen_file', 'vege_file', 'soil_file', \
             'oput_dir', 'N.of.days', 'SpinUp', 'IncSpin']

    spa_cfg = pd.DataFrame({'value': file_list}, index=index_l)

    return spa_cfg

def create_dir(dirpath, soilfile):
    """
    Creates the output directory for the experiment
    """
    depth_str = os.path.splitext(soilfile)[0].split("_")[-1]
    output_dir = "{0}/outputs_{1}".format(dirpath, depth_str)
    return output_dir

def get_simlen(fname):
    """
    Fast count of the number of lines in an ascii formatted file
    """
    with open(fname) as f:
        for (i, _) in enumerate(f):
            pass
    return i + 1

def main():

    # Get the number of available cores for multi-proc
    num_cores = cpu_count()

    spa_dirs = [os.path.join(dp, d) for (dp, dn, _) in os.walk(DIRPATH) \
                    for d in dn if d not in ['inputs']]

    # get the phenology paths
    phen_paths = [os.path.join(dp, f) \
                    for (dp, _, fn) in os.walk(DIRPATH) if fn \
                        for f in fn if re.search(r'phenology', f)]

    # find the length of the files -> equates to the length of the simulation
    site_dlens = [get_simlen(ppath) for ppath in phen_paths




    test = [dn for (_, dn, _) in os.walk(DIRPATH, topdown=True) for d in dn if d == 'inputs']
    print(test)
    return




    # collect all input file names for each site
    input_files = [[f for f in fn if re.search(r'^((?!soils|DS_Store).)*.csv$', f)] \
                    #for d in dn if d in ['inputs'] \
                    for (_, dn, fn) in os.walk(DIRPATH, topdown=True) \
                    for d in dn]

    # now take from the list only the soil profiles
    soilprof_files = [[f for f in fn if re.search(r'soils', f)] \
                    for (_, dn, fn) in os.walk(DIRPATH) if not dn]

    # create output directories for each soil profile at each site
    soil_outdirs = [[create_dir(sdir, soil) for soil in soil_list] \
                    for (sdir, soil_list) in zip(spa_dirs, soilprof_files)]

    # create all the SPA configuration files for the batch run
    cfg_files = [[create_inf_file(inp, soil, sdir, simlen) \
                    for (soil, sdir) in zip(soil_list, dir_list)] \
                    for (inp, soil_list, dir_list, simlen) \
                        in zip(input_files, soilprof_files, soil_outdirs, site_dlens)]

    print(input_files)
    return 1

    # run all the SPA rooting depth experiments
    for (cur_dir, cfg_ex, soil_out) in zip(spa_dirs, cfg_files, soil_outdirs):
        # message to echo to user
        print("Running root experiments for site: {0}".format(os.path.basename(cur_dir)))
        # run each experiment based off the configuration file
        for (cfg, out_dir) in zip(cfg_ex, soil_out):

            # create the file path to write to for the configuration file
            fpath = cur_dir + '/file_names.txt'
            # if it already exists remove it; prevents appending to the file ad infinitum
            if os.path.exists(fpath):
                os.remove(fpath)
                cfg.to_csv(fpath, header=None, index=True, sep=',', mode='w')
            else:
                cfg.to_csv(fpath, header=None, index=True, sep=',', mode='w')

            # create the output directories if they don't already exit
            if not os.path.exists(out_dir):
                os.makedirs(out_dir)

            # now execute the binary
            run_SPAbin(cur_dir)


    # Execute each binary collected from the above path search
#    for (i, bin_fpath) in enumerate(binary_locs):
#        run_SPAbin(bin_fpath, i+1)

    # Execute each binary collected from the above path search
#    Parallel(n_jobs=num_cores)(delayed(run_SPAbin)(bin_path, i+1) \
#        for (i, bin_path) in enumerate(binary_locs))

    return 1

if __name__ == "__main__":

    DIRPATH = os.path.expanduser("~/Savanna/Models/SPA1/outputs/rootexp/rootdepth/")

    BINLOC = os.path.expanduser("~/Repositories/EcoModels/spa/src/SPA_i686")

    main()

