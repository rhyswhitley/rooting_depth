#!/usr/bin/env python3

from scipy.optimize import minimize

import os
import pandas as pd
import numpy as np

def root_biomass(rz_a, rz_b, kd, rdsurf):
    """
    rz_a = root layer at depth A
    rz_b = root layer at depth B
    kd = decay coefficient
    rdsurf = root density at the surface
    """
    return rdsurf*(-1/kd*np.exp(-kd*rz_a) + 1/kd*np.exp(-kd*rz_b))

def root_biom_vx(kd_, rd_surf_, x_depth):
    return [root_biomass(x_depth[i], x_depth[i-1], kd_, rd_surf_) \
                          for i in range(len(x_depth)) if i > 0]

def cost_fun(kd_, rd_surf_, r_totbiom_, x_depth):

    root_bk = root_biom_vx(kd_, rd_surf_, x_depth)

    return abs(r_totbiom_ - sum(root_bk))

def get_root_dist(rb_total_, rd_surf_, x_depth_):

    # now find the coefficient
    res = minimize(cost_fun, 5, args=(rd_surf_, rb_total_, x_depth_), \
                   method='nelder-mead', options={'xtol': 1e-8, 'disp': False})
    # optimal parameter
    kd_0 = res.x[0]

    # now we can determine the root density per soil layer
    root_biom_x = root_biom_vx(kd_0, rd_surf_, x_depth_)
    # return to user
    return root_biom_x

def create_profile(rb_total, rd_surf, slayers, soil_text):
    """
    Creates a set of soil profile input file for the first set of rooting
    depth experiments that adjusts the rooting depth at each NATT site.
    """

    # default number of soil layers
    nlayer = 20

    # create 20 soil layers
    soil_layer = ["layer_{0}".format(i+1) for i in range(nlayer)] + ['core']

    # set the thicknesses of the 20 soil layers, divide accordingly to match with field measurements
    soil_thick = [0.1 if z < 4 else 0.2 if (z >= 4) & (z < 12) \
                  else 0.1 if (z == len(soil_layer) - 1) else 1.0 \
                  for (z, _) in enumerate(soil_layer)]

    # calculate soil depth (include the surface)
    soil_depth = np.cumsum([0] + soil_thick)

    # create a dataframe to work with
    soil_df = pd.DataFrame({'thick': soil_thick, 'depth': soil_depth[:-1], \
                            'layer': soil_layer})

    # determine the root biomass for each soil layer
    root_biom = get_root_dist(rb_total, rd_surf, soil_depth)

    # set the relative fraction of roots per soil layer
    soil_df['root_frac'] = [rb/rb_total if (z <= slayers) else 0 \
                            for (z, rb) in enumerate(root_biom)]

    # add soil texture information
    soil_df['sand%'] = [soil_text['Bs'] if (np.float32(soil_text['ABz']) <= z) \
                        else soil_text['As'] for z in soil_depth[:-1]]
    soil_df['clay%'] = [soil_text['Bc'] if (np.float32(soil_text['ABz']) <= z) \
                        else soil_text['Ac'] for z in soil_depth[:-1]]

    # fraction above porosity at which drainage occurs
    soil_df['draincheck'] = 0.5

    # not really relevant:: the soil texture information above will dominate any hydraulic movement
    soil_df['organic_frac'] = 0.08
    soil_df['mineral_frac'] = 0.92

    # add extra information for state variables
    soil_df['init_water_frac'] = 0.1
    soil_df['init_soil_temp'] = 288.15
    soil_df['init_ice_prop'] = 0.0

    # remainder parameter information
    soil_df['rootrad'] = 1e-4
    soil_df['rootl'] = slayers + 1
    soil_df['biomass'] = rb_total
    soil_df['snowweight'] = 0
    soil_df['snowheight'] = 0
    # parameter information doesn't need to be repeated
    soil_df.ix[1:, range(11, 17)] = None

    # set the index (not really necessary)
    soil_df.set_index(['layer'], inplace=True)

    return soil_df

def reorder_transpose(soil_df, site):

    # The ordering the columns need to be in to
    col_order = [1, 6, 7, 8, 9, 10, 2, 12, 13, 3, 4, 11, 5, 14, 15]

    # transpose to match the format need to be read in by the model
    new_df = soil_df.ix[:, col_order].T

    # get depth on max root index
    depth = np.round(soil_df.depth[soil_df.rootl[0].astype(int)], 2)

    # create filepath string
    fpath = "{0}/{1}/inputs/{1}_soils_depth{2}m.csv".format(SAVEPATH, site, depth)

    # save to file
    new_df.to_csv(fpath, sep=",") #, line_terminator=LT)

    return None

def main():

    # soil texture information for each site
    HS_soiltex = {'As': 60, 'Ac': 15, 'Bs': 30, 'Bc': 40, 'ABz': 0.30}
    AR_soiltex = {'As': 50, 'Ac': 20, 'Bs': 20, 'Bc': 55, 'ABz': 0.30}
    DA_soiltex = {'As': 40, 'Ac': 20, 'Bs': 30, 'Bc': 35, 'ABz': 0.20}
    DR_soiltex = {'As': 25, 'Ac': 50, 'Bs': 20, 'Bc': 55, 'ABz': 0.15}
    SP_soiltex = {'As': 40, 'Ac': 20, 'Bs': 30, 'Bc': 35, 'ABz': 0.20}
    # store in list
    soiltex_natt = [HS_soiltex, AR_soiltex, DA_soiltex, DR_soiltex, SP_soiltex]

    # two for passes required: 1) Site -> 2) root layer depth
    set_slayers = [6] + list(range(11, 17))

    # create a set of soil profiles for each site for variable rooting depths
    profile_list = {site: [create_profile(1980, 10e3, slay, stex) \
                        for slay in set_slayers] \
                        for (site, stex) in zip(SITES, soiltex_natt)}

    # transpose the dataframes into a format that SPA will recognise and write to file
    _ = [[reorder_transpose(df, site) for df in slay_list] \
                    for (site, slay_list) in profile_list.items()]

if __name__ == "__main__":

    SAVEPATH = os.path.expanduser("~/Savanna/Models/SPA1/outputs/rootexp/rootdepth")

    SITES = ["HowardSprings", "AdelaideRiver", "DalyUncleared", "DryRiver", \
             "SturtPlains"]

    # line terminator for CSV (doesn't seem to be necessary)
    LT = '\r\n'

    main()
