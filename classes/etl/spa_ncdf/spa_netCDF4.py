#!/usr/bin/env python2

class spa_netCDF4(object):
    def __init__(self):
        self.start_date = "2001-01-01 00:00:00"

    def assign_variables(self, nc_obj):
        # CREATE DIMENSIONS
        nc_obj.createDimension('x', 1)
        nc_obj.createDimension('y', 1)
        nc_obj.createDimension('z', 1)
        nc_obj.createDimension('soil', 20)
        nc_obj.createDimension('time', None)
        # CREATE VARIABLES
        nc_obj.createVariable('x', 'f8', ('x'))
        nc_obj.createVariable('y', 'f8', ('y'))
        nc_obj.createVariable('latitude', 'f8', ('x', 'y'))
        nc_obj.createVariable('longitude', 'f8', ('x', 'y'))
        nc_obj.createVariable('time', 'f8', ('time'))
        nc_obj.createVariable('soildepth', 'f8', ('soil'))
        # >> Meteorology
        nc_obj.createVariable('SWdown', 'f8', ('time', 'x', 'y'))
        nc_obj.createVariable('Tair', 'f8', ('time', 'x', 'y'))
        nc_obj.createVariable('VPD', 'f8', ('time', 'x', 'y'))
        nc_obj.createVariable('Cair', 'f8', ('time', 'x', 'y'))
        nc_obj.createVariable('Rainfall', 'f8', ('time', 'x', 'y'))
        # >> Land-surface fluxes
        nc_obj.createVariable('GPP', 'f8', ('time', 'x', 'y'))
        nc_obj.createVariable('AutoResp', 'f8', ('time', 'x', 'y'))
        nc_obj.createVariable('NPP', 'f8', ('time', 'x', 'y'))
        nc_obj.createVariable('Qle', 'f8', ('time', 'x', 'y'))
        nc_obj.createVariable('TVeg', 'f8', ('time', 'x', 'y'))
        nc_obj.createVariable('Esoil', 'f8', ('time', 'x', 'y'))
        nc_obj.createVariable('Ecanop', 'f8', ('time', 'x', 'y'))
        # >> Vegetation
        nc_obj.createVariable('Etree', 'f8', ('time', 'x', 'y'))
        nc_obj.createVariable('Egrass', 'f8', ('time', 'x', 'y'))
        nc_obj.createVariable('Atree', 'f8', ('time', 'x', 'y'))
        nc_obj.createVariable('Agrass', 'f8', ('time', 'x', 'y'))
        nc_obj.createVariable('Rtree', 'f8', ('time', 'x', 'y'))
        nc_obj.createVariable('Rgrass', 'f8', ('time', 'x', 'y'))
        nc_obj.createVariable('Gtree', 'f8', ('time', 'x', 'y'))
        nc_obj.createVariable('Ggrass', 'f8', ('time', 'x', 'y'))
        nc_obj.createVariable('LAItree', 'f8', ('time', 'x', 'y'))
        nc_obj.createVariable('LAIgrass', 'f8', ('time', 'x', 'y'))
        # >> Soils
        nc_obj.createVariable('SWC20', 'f8', ('time', 'x', 'y'))
        nc_obj.createVariable('SWC80', 'f8', ('time', 'x', 'y'))
        nc_obj.createVariable('IntSWC', 'f8', ('time', 'x', 'y'))
        nc_obj.createVariable('IntSWP', 'f8', ('time', 'x', 'y'))
        #
        nc_obj.createVariable('SoilMoist', 'f8', ('time', 'soil', 'y'))
        return None

    def assign_units(self, nc_obj):
        # ASSIGN UNITS
        # >> [Dimensions]
        nc_obj.variables['x'].units = ""
        nc_obj.variables['y'].units = ""
        nc_obj.variables['latitude'].units = "degrees_north"
        nc_obj.variables['longitude'].units = "degrees_east"
        nc_obj.variables['time'].units = "seconds since " + self.start_date
        nc_obj.variables['soildepth'].units = "m"
        # >> [Time-varying values]
        # Land-surface
        nc_obj.variables['NPP'].units = "umol/m^2/s"
        nc_obj.variables['GPP'].units = "umol/m^2/s"
        nc_obj.variables['AutoResp'].units = "umol/m^2/s"
        nc_obj.variables['Qle'].units = "W/m^2"
        nc_obj.variables['TVeg'].units = "W/m^2"
        nc_obj.variables['Esoil'].units = "W/m^2"
        nc_obj.variables['Ecanop'].units = "W/m^2"
        # Vegetation
        nc_obj.variables['Etree'].units = "mmol/m^2/s"
        nc_obj.variables['Egrass'].units = "mmol/m^2/s"
        nc_obj.variables['Atree'].units = "umol/m^2/s"
        nc_obj.variables['Agrass'].units = "umol/m^2/s"
        nc_obj.variables['Rtree'].units = "umol/m^2/s"
        nc_obj.variables['Rgrass'].units = "umol/m^2/s"
        nc_obj.variables['Gtree'].units = "mmol/m^2/s"
        nc_obj.variables['Ggrass'].units = "mmol/m^2/s"
        nc_obj.variables['LAItree'].units = "m^2/m^2"
        nc_obj.variables['LAIgrass'].units = "m^2/m^2"
        # Soil Profile
        nc_obj.variables['SWC20'].units = "m^3/m^-3"
        nc_obj.variables['SWC80'].units = "m^3/m^-3"
        nc_obj.variables['IntSWC'].units = "m^3/m^-3"
        nc_obj.variables['IntSWP'].units = "MPa"
        nc_obj.variables['SoilMoist'].units = "m^3/m^-3"
        return None

    def assign_longNames(self, nc_obj):
        # LONG NAMES
        nc_obj.variables['NPP'].longname = "Net primary productivity"
        nc_obj.variables['GPP'].longname = "Gross primary productivity"
        nc_obj.variables['AutoResp'].longname = "Autotrophic respiration"
        nc_obj.variables['Qle'].longname = "Latent heat flux"
        nc_obj.variables['TVeg'].longname = "Vegetation transpiration"
        nc_obj.variables['Esoil'].longname = "Soil evaporation"
        nc_obj.variables['Ecanop'].longname = "Wet canopy evaporation"
        # Vegetation
        nc_obj.variables['Etree'].longname = "C3 transpiration"
        nc_obj.variables['Egrass'].longname = "C4 transpiration"
        nc_obj.variables['Atree'].longname = "C3 photosynthesis"
        nc_obj.variables['Agrass'].longname = "C4 photosynthesis"
        nc_obj.variables['Rtree'].longname = "C3 respiration"
        nc_obj.variables['Rgrass'].longname = "C4 respiration"
        nc_obj.variables['Gtree'].longname = "C3 stomatal conductance"
        nc_obj.variables['Ggrass'].longname = "C4 stomatal conductance"
        nc_obj.variables['LAItree'].longname = "C3 leaf area index"
        nc_obj.variables['LAIgrass'].longname = "C4 leaf area index"
        # Soil Profile
        nc_obj.variables['SWC20'].longname = "Soil water content at 20 cm depth"
        nc_obj.variables['SWC80'].longname = "Soil water content at 80 cm depth"
        nc_obj.variables['IntSWC'].longname = "Integrated average soil water content"
        nc_obj.variables['IntSWP'].longname = "Integrated average soil water potential"
        nc_obj.variables['SoilMoist'].longname = "Soil water content of the profile"
        return None
