#!/usr/bin/env python3

import pickle
import os, errno, re
import numpy as np
import pandas as pd

from collections import OrderedDict
from natsort import natsorted

def create_phen_file(lai_ts):
    total = lai_ts.resample('D').mean()
    lai_veg = treegrass_frac(total, 30)
    site_phen = divide_leaves(lai_veg)

    site_phen['rootbiomass'] = [3840.]*len(total)
    site_phen['nitrogen'] = [lai*2 for lai in total]
    site_phen['TimeStep'] = total.index.dayofyear

    site_phen_rot = site_phen.iloc[:, [13, 11] + list(range(0, 11)) + [12]]
    site_phen_out = pd.concat([site_phen_rot, site_phen.iloc[:, 1:11]], axis=1)
    site_phen_out.columns = ['DOY', 'rootbiom', 'lai'] \
        + ['lfr_{0}'.format(i) for i in range(1, 11)] \
        + ['nit'] \
        + ['nfr_{0}'.format(i) for i in range(1, 11)]
    return site_phen_out

def divide_leaves(lai_part):
    # was 0.7 (trees) and 1.3 (grass) obviously wrong
    trees_frac = [lai_part['tree']/lai_part['total']/5. for i in range(5)]
    grass_frac = [lai_part['grass']/lai_part['total']/5. for i in range(5)]
    leaf_alloc = pd.concat([lai_part['total']] + trees_frac + grass_frac, axis=1)
    return leaf_alloc

def treegrass_frac(ndvi, day_rs):
    """
    Process based on Donohue et al. (2009) to separate out tree and grass cover,
    using moving windows (adapted here for daily time-step)
    """
    # first calculate the 7-month moving minimum window across the time-series
    # changed period to 3 to kill grass in dry season
    fp1 = moving_something(np.min, ndvi, period=7, day_rs=day_rs)
    fp2 = moving_something(lambda x: sum(x)/(9*day_rs), fp1, period=9, day_rs=day_rs)
    fr1 = ndvi - fp2

    ftree = [p2 - np.abs(r1) if r1 < 0 else p2 for p2, r1 in zip(fp2, fr1)]
    fgrass = ndvi - ftree

    return pd.DataFrame({'total':ndvi, 'tree':ftree, 'grass':fgrass})

def moving_something(_fun, tseries, period, day_rs=16, is_days=True):
    """
    Applies a function to a moving window of the time-series:
    ft_ = function([ f(t-N), f(t). f(t+N)])
    """
    # if the time-series is at a day-time step, update the window to a step-size of 16 days
    if is_days:
        p0 = period*day_rs
    else:
        p0 = period

    # find upper and lower bounds of the moving window
    half = p0//2
    tlen = len(tseries)

    twin = [0]*tlen
    for im in range(tlen):
        # find the something for the window that satisfy the edge conditions
        if im < half:
            # fold back onto the end of the time-series
            twin[im] = _fun(np.hstack([tseries[tlen-(half-im):tlen],\
                                        tseries[0:im+half]]))
        elif im > tlen-half:
            # fold back into the beginning of the time-series
            twin[im] = _fun(np.hstack([tseries[im-half:tlen],\
                                        tseries[0:half-(tlen-im)]]))
        else:
            twin[im] = _fun(tseries[im-half:im+half])

    return twin

def import_tower_data(file_name):
    """
    Imports dingo meteorology file and puts an index on the time column
    """
    # read in file
    dataset = pd.read_csv(file_name, parse_dates=True, index_col=['DT'])

    # add PAR data
    dataset["PAR"] = dataset["Fsd_Con"]*2.3

    # add timestep column for light interception geometric calculations
    dataset["TimeStep"] = [d.dayofyear + (d.hour + d.minute/60.)/24 \
                              for d in dataset.index]

    # return a dictionary of both drivers and
    driver_columns = ["TimeStep", "Ta_Con", "Cc", "Ws_CSAT_Con", "Fsd_Con", \
                      "VPD_Con", "PAR", "Precip_Con", "Lai_1km_new_smooth"]

    valid_columns = ["GPP_Con", "Fe_Con", "Fre_Con", "Fc_ustar", "Sws_Con", \
                     "Fe_Con_QCFlag", "Fc_Con_QCFlag"]

    # backfill NaN and remove remaining rows with NaN
    data_nonull = dataset.ix[:, driver_columns + valid_columns] \
        .interpolate().fillna(method='bfill') \
        [pd.notnull(dataset["Lai_1km_new_smooth"])]

    # get driver and target dataframes and backfill NaN values
    drivers = data_nonull.ix[:, driver_columns].apply(remove_nonphysical)
    targets = data_nonull.ix[:, valid_columns]

    return {'drivers': drivers, 'targets': targets}

def remove_nonphysical(dstream, win=30):
    """
    Non-physical values are removed by comparing values with the max of a 95% CI
    monthly moving window.

    Negative values are set to 0
    """
    # rolling mean
    #mean = pd.rolling_mean(dstream, window=win*48).fillna(method="bfill")
    mean = dstream.rolling(window=win*48).mean().fillna(method="bfill")
    # rolling standard deviation
    #std = pd.rolling_std(dstream, window=win*48).fillna(method="bfill")
    std = dstream.rolling(window=win*48).std().fillna(method="bfill")
    # determined rolling ci
    ci99 = [m + 2.5*s for (m, s) in zip(mean, std)]
    # max CI 99
    top_val = np.max(ci99)
    # clean values
    #dstream_clean = [np.min([top_val, ds[i]]) for (i, ds) in enumerate(dstream)]
    dstream_clean = np.minimum(top_val, np.maximum(dstream, 0))
    # return cleaned data stream
    return dstream_clean

def ensure_dir(path):
    # Create folders for storage if they done exist
    try:
        if not os.path.exists(path):
            os.makedirs(path)
    except OSError as exc: # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else: raise

def get_site_name(filepath):
    base = os.path.basename(filepath)
    split = base.split("_")
    return split[-2]

def main():

    print("Loading datasets from tower files\n")

    # find the filepaths for the tower data for the SITES we are interested
    file_paths = [os.path.join(fp, f) for (fp, _, fn) in os.walk(DIRPATH) \
        for f in fn if re.search("|".join(SITES), f)]

    # load the data, and split into SPA driver and target information
    data_sets = OrderedDict(natsorted({get_site_name(fn): import_tower_data(fn) \
                 for fn in file_paths}.items()))

    print("Creating each site's phenology profile\n")

    # universal phenology file (LAI inf. is always last column)
    site_phen = OrderedDict(natsorted({site: create_phen_file(ds['drivers'].ix[:, -1]) \
                 for (site, ds) in data_sets.items()}.items()))

    print("Writing files to experiment folders\n")

    # write the corresponding input file information for each site
    for ((slab1, dataset), (slab2, phen)) \
        in zip(data_sets.items(), site_phen.items()):
            if slab1 == slab2:
                print("Writing input files for site: {0}".format(slab1))
                input_path = "{0}/{1}/inputs".format(SAVEPATH, slab1)
                # if the path does not exist the create it
                ensure_dir(input_path)
                # save meteorology file
                dataset['drivers'].ix[:, :-1] \
                    .to_csv("{0}/{1}_met_drivers.csv".format(input_path, slab1), \
                            sep=",", index=False) #, line_terminator=LT)
                # save phenology file
                phen.to_csv("{0}/{1}_phenology.csv".format(input_path, slab2), \
                            sep=",", index=False) #, line_terminator=LT)
            else:
                print("{0} != {1}".format(slab1, slab2))
                print("ERROR: dictionaries not matched correctly!")
                break

    # save the data_sets dictionary object for later use in emulation and testing
    pickle.dump(data_sets, open("{0}natt_datasets.pkl".format(SAVEPATH), 'wb'), \
                protocol=2)

    return 1

if __name__ == "__main__":

    DIRPATH = os.path.expanduser("~/Google Drive/Dingo_v11/")
    SAVEPATH = os.path.expanduser("~/Savanna/Models/SPA1/outputs/rootexp/rootdepth/")

    SITES = ["HowardSprings", "AdelaideRiver", "DalyUncleared", "DryRiver", \
             "SturtPlains"]

    # line terminator for CSV
    LT = '\r\n'

    main()


