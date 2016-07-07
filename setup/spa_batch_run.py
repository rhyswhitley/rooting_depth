#!/usr/bin/env python2.7

import os, re
import subprocess as sub
from joblib import Parallel, delayed
from multiprocessing import cpu_count

def run_SPAbin(bin_fpath, nex):
    # Echo to user
    print("Running SPA simulation for site {0}".format(nex))

    # take the head of the file path to tell the subproces which directory
    # to execute in
    bin_cwd = os.path.split(bin_fpath)[0]

    # execture the SPA Fortran binary
    proc = sub.Popen(bin_fpath, cwd=bin_cwd, bufsize=1, shell=False,
                        close_fds=False)

    # Slight pause to allow python to catchup with the terminal
    proc.wait()

def main():

    # Get the number of available cores for multi-proc
    num_cores = cpu_count()

    # Get paths to each SPA binary under the set directory path
    binary_locs = [os.path.join(dp, f) for (dp, _, fn) \
                    in os.walk(DIRPATH) for f in fn \
                    if re.match('SPA', f)]

    # Execute each binary collected from the above path search
#    for (i, bin_fpath) in enumerate(binary_locs):
#        run_SPAbin(bin_fpath, i+1)

    # Execute each binary collected from the above path search
    Parallel(n_jobs=num_cores)(delayed(run_SPAbin)(bin_path, i+1) \
        for (i, bin_path) in enumerate(binary_locs))

    return 1

if __name__ == "__main__":

    DIRPATH = os.path.expanduser("~/Savanna/Models/SPA1/outputs/baseline/")

    main()

