import numpy as np
import natsort
import pandas as pd
import sys, os, re
import time
import netCDF4 as nc

from scipy import integrate

from etl.spa_ncdf import spa_netCDF4

__author__ = 'Rhys Whitley'
__email__ = 'rhys.whitley@gmail.com'
__created__ = "2015-10-15"
__modified__ = time.strftime("%c")
__version__ = '1.0'
__status__ = 'prototype'

class SPAoutput_ETL(object):
    """
    A set of functions that perform pre-processing routines on the raw SPA v1
    output CSV files. These files are imported into pandas dataframes at a
    30min or 1hr time-step. Each output file is saved as a dataframe and then
    pickled for later use.

        dirpath = directory where spa outputs are located
        savepath = directory where pickle and CSVs will be saved

    The file processed as follows:
    1. hourly.csv
    2. l1-10.csv
    3. soilwater.csv
    4. upfrac.csv
    """
    def __init__(self):
        super(SPAoutput_ETL, self).__init__()
        self.start_date = "2001-01-01"
        self.time_freq = "30min"
        self.soil_depth = None
        self.glay = 5 # index in canopy layers where grasses begin

    def reduce_canopy(self, data_list):
        """
        Generalise the 10 canopy layer files into a dictionary of tree and
        grass dataframes
        """
        # split layer files into trees and grasses and
        # concatenate these together
        trees_raw = pd.concat(data_list[:self.glay], axis=0)
        grass_raw = pd.concat(data_list[self.glay:], axis=0)

        # columns to be summed
        not_col = lambda x, y: x.columns - x.columns[y]
        sums_col = [1, 2, 6, 8, 9, 13]
        mean_col = not_col(trees_raw, sums_col)

        # lambda for quick convert and group
        convgrp = lambda x: x.convert_objects(convert_numeric=True)\
            .groupby(level=0)

        # do separate groupbys for sums and means so we can exploit
        # the faster cython operability

        # trees
        trees_sums = convgrp(trees_raw.ix[:, sums_col]).sum()
        trees_mean = convgrp(trees_raw.ix[:, mean_col]).mean()
        trees_grp = pd.concat([trees_sums, trees_mean], axis=1)

        # grasses
        grass_sums = convgrp(grass_raw.ix[:, sums_col]).sum()
        grass_mean = convgrp(grass_raw.ix[:, mean_col]).mean()
        grass_grp = pd.concat([grass_sums, grass_mean], axis=1)

        return {'trees' : trees_grp, 'grass' : grass_grp}

    def flatten_dictlist(self, dlist):
        """Flattens the tree/grass dictionary into a single dataframe"""
        # concatenate the two dataframes
        tg_df = pd.concat(dlist, axis=1)
        # set new column names
        tg_df.columns = ["{0}_{1}".format(h1, h2) for (h1, h2) in tg_df.columns]
        # return to user
        return tg_df

    def import_spa_canopy(self, file_name, clayer):
        """
        Imports SPA canopy layer outputs into a pandas dataframe
        * adds indexes on time and canopy layer
        """
        # read in file
        data = pd.read_csv(file_name, sep=r',\s+', engine='python', na_values=["Infinity"])
        # remove header whitespace
        data.columns = [dc.replace(' ', '') for dc in data.columns]
        # add layer label (need to reduce on this)
        data['clayer'] = clayer
        # re-do the datetime column as it is out by 30 minutes
        data['DATE'] = pd.date_range(start=self.start_date, periods=len(data), \
                                    freq=self.time_freq)
        # set dataframe index on this column
        data.drop('Time', axis=1, inplace=True)
        # set new indices
        data_dated = data.set_index(['DATE', 'clayer'])
        # return to user
        return data_dated.bfill()

    def import_spa_output(self, file_name):
        """
        Imports hourly SPA land-surface flux outputs with the index on the time
        column
        """
        # read in file
        data = pd.read_csv(file_name, sep=r',', engine='python', na_values=["Infinity"])
        # remove header whitespace
        data.columns = [dc.replace(' ', '') for dc in data.columns]
        # re-do the datetime column as it is out by 30 minutes
        data['DATE'] = pd.date_range(start=self.start_date, periods=len(data), \
                                    freq=self.time_freq)
        # set dataframe index on this column
        data_dated = data.set_index(['DATE'])
        # return to user
        return data_dated.bfill()

    def import_soil_profile(self, file_name):
        # import the soil profile input file
        soil_prof = pd.read_csv(file_name)
        # take the first row which has the soil layer thicknesses
        soil_layer = soil_prof.ix[0, 1:-1].values.astype('float')
        # calculate the cumulative thickness, equivalent to soil depth at each layer
        cumul_depth = np.cumsum(soil_layer)
        # return to user
        return cumul_depth
    def load_gasex_raw(self, filepaths):
        """
        Special wrapper function to perform a mass import of the SPA canopy
        layer outputs from a simulation folders. The canopy layer files are
        individually imported and then aggregated into tree and grass
        components using the reduce_canopy and flatten functions. These are
        then added with the land-surface outputs to simulation ncdf files.
        """
        # Import the canopy layer output files from each model experiment
        canopy_data = [self.import_spa_canopy(fname, i) \
                        for (i, fname) in enumerate(filepaths)]

        # Reduce and flatten the canopy layers into tree and grass outputs
        treegrass_data = self.reduce_canopy(canopy_data)
        treegrass_flat = self.flatten_dictlist(treegrass_data)

        # Return to user
        return treegrass_flat

    def load_hourly_raw(self, filepaths):
        """
        Wrapper function to perform a mass import of SPA outputs from a
        simulation files using Pandas. Returns a dictionary that can be used
        later to process outputs into ncdf file.
        """
        # import the land-surface output files from each model experiment
        land_surf = {fname.split(".")[0].split("/")[-1]: \
                    self.import_spa_output(fname) \
                    for fname in filepaths}

        return land_surf

    def get_spa_filepaths(self, dir_path):
        """
        Given the top directory, find all subfolders that contain
        outputs from the SPA 1 model.

        regex_remove: controls the files that are ignored in the scan
        regex_take: controls the subfolders that are search for SPA outputs

        Returns file_list, which contains the full paths for each output
        for each SPA simulation folder.
        """
        # Files to be ignored in the search
        files_remove = ['canopy','DS_Store','ci','drivers','energy','iceprop', \
                        'parcheck','power','waterfluxes','test','solar','daily']
        # Turn list into a regex argument
        regex_remove = r'^((?!{0}).)*$'.format("|".join(files_remove))
        # Subfolders to be searched
        regex_take = r'(\binputs\b)|(\boutputs\b)$'

        # Number of subfolders
        n_subfold = len(regex_take.split("|"))

        # Walk down through the directory path looking for SPA files
        raw_list = [[os.path.join(dp, f) for f in fn \
                    if re.search(regex_remove, f)] \
                    for (dp, dn, fn) in os.walk(dir_path, followlinks=False) \
                    if re.search(regex_take, dp)]

        # Exit program if no files found
        if len(raw_list) == 0:
            # Echo user
            sys.stderr.write('No output files found. Check your FilePath. \n')
            sys.exit(-1)
        else:
            # Attach the two subfolders in each simulation together
            file_list = [natsort.natsorted(raw_list[i] + raw_list[i+1]) \
                        for i in np.arange(0, len(raw_list), n_subfold)]
            # Pass back to user
            return natsort.natsorted(file_list)

    def integrated_mean_soil(self, soil_matrix, soil_depths):

        soil_int = lambda y, x: integrate.trapz(y=y, x=x)/max(x)
        int_soil_prof = soil_matrix.apply(soil_int, axis=1, args=(soil_depths,))

        return int_soil_prof

    def process_outputs(self, fpack):
        """
        Given a list of files from a simulation folder load them using pandas
        and process them into netcdf4 files
        """
        fpath_list, nc_fout = fpack
        # Static soil profile path
        profpaths = [fp for fp in fpath_list if \
                     re.search(r'_soils', os.path.basename(fp))]
        # Paths for soil profile water fluxes
        soilpaths = [fp for fp in fpath_list \
                     if re.search(r'soil(?!s)|upfrac', os.path.basename(fp))]
        # Canopy layer outputs
        canopaths = [fp for fp in fpath_list \
                     if re.search(r'l.*[0-9]+.csv$', os.path.basename(fp))]
        # Paths for land-surface fluxes
        landpaths = [fp for fp in fpath_list \
                     if re.search(r'^((?!l\d+|soil|upfrac|temp|phen).)*$', os.path.basename(fp))]

        # Canopy layer outputs
        canopy10 = self.load_gasex_raw(canopaths)
        # Land-surface and meteorology fluxes
        landsurf = self.load_hourly_raw(landpaths)
        # Below-ground water fluxes
        watrprof = self.load_hourly_raw(soilpaths)
        # soil layer thicknesses
        soil_depths = self.import_soil_profile(profpaths[0])
        # Determine the integrated average soil state properties
        intProf_swc = self.integrated_mean_soil(watrprof['soilwater'].ix[:, 1:-1], soil_depths)

        print("> writing netCDF file: {0}".format(nc_fout))
        #import ipdb; ipdb.set_trace()

        # Open a NCDF4 file for SPA simulation outputs
        nc_file = nc.Dataset(nc_fout, 'w', format='NETCDF4')

        # Assign attributes
        ncdf_attr = spa_netCDF4()
        ncdf_attr.assign_variables(nc_file)
        ncdf_attr.assign_units(nc_file)
        ncdf_attr.assign_longNames(nc_file)

        # Assign values to variables
        tseries = pd.timedelta_range(0, periods=len(landsurf['hourly']), freq="1800s") \
                    .astype('timedelta64[s]')

        # Get time from netcdf driver file
        nc_file.variables['time'][:] = tseries.values
        nc_file.variables['soildepth'][:] = soil_depths
        # [Land-surface fluxes]
        nc_file.variables['GPP'][:] = landsurf['hourly']['gpp'].values
        nc_file.variables['Qle'][:] = landsurf['hourly']['lemod'].values
        nc_file.variables['TVeg'][:] = landsurf['hourly']['transle'].values
        nc_file.variables['Esoil'][:] = landsurf['hourly']['soille'].values
        nc_file.variables['Ecanop'][:] = landsurf['hourly']['wetle'].values
        nc_file.variables['AutoResp'][:] = landsurf['hourly']['resp'].values
        # [Vegetation fluxes]
        nc_file.variables['Atree'][:] = canopy10['trees_Ag'].values
        nc_file.variables['Agrass'][:] = canopy10['grass_Ag'].values
        nc_file.variables['Rtree'][:] = canopy10['trees_Rd'].values
        nc_file.variables['Rgrass'][:] = canopy10['grass_Rd'].values
        nc_file.variables['Etree'][:] = canopy10['trees_Et'].values
        nc_file.variables['Egrass'][:] = canopy10['grass_Et'].values
        nc_file.variables['Gtree'][:] = canopy10['trees_gs'].values
        nc_file.variables['Ggrass'][:] = canopy10['grass_gs'].values
        nc_file.variables['LAItree'][:] = canopy10['trees_LAI'].values
        nc_file.variables['LAIgrass'][:] = canopy10['grass_LAI'].values
        # [Soil Profile]
        nc_file.variables['SWC20'][:] = watrprof['soilwater']['w1'].values
        nc_file.variables['SWC80'][:] = watrprof['soilwater']['w4'].values
        nc_file.variables['IntSWC'][:] = intProf_swc.values
        nc_file.variables['IntSWP'][:] = watrprof['soilwater']['w_swp'].values
        nc_file.variables['SoilMoist'][:] = watrprof['soilwater'].ix[:, 1:-1].values

        # Close file
        nc_file.close()

