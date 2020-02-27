# Useful functions to work with radiosonde data measured during EUREC4A 

import numpy as np
from netCDF4 import Dataset
import glob
import xarray
import os
import matplotlib.pyplot as plt
import typhon
import copy
from scipy.interpolate import interp1d

def get_filenames(dates, start_time, end_time, platforms, branches, datapath):
    """
    Returns list of filenames of all available soundings for the specified days, time periods and platforms.
    """
    filenames = []
        
    for platform in platforms:
        if platform == 'ATL':
            time_ind = slice(-9, -5)
        else:
            time_ind = slice(-7, -3)
        for date in dates:
            if 'ascending' in branches:
                filenames_all = glob.glob(os.path.join(datapath, f'{platform}_SoundingAscentProfile_*_{date}_*.nc'))
                if filenames_all:
                    filenames_time = [f for f in filenames_all if int(f[time_ind]) > int(start_time) and int(f[time_ind]) < int(end_time)]
                    filenames.extend(filenames_time)
            if 'descending' in branches:
                filenames_all = glob.glob(os.path.join(datapath, f'{platform}_SoundingDescentProfile_*_{date}_*.nc'))
                if filenames_all:
                    filenames_time = [f for f in filenames_all if int(f[time_ind]) > int(start_time) and int(f[time_ind]) < int(end_time)]
                    filenames.extend(filenames_time)
    return filenames

def profiles_from_netcdf(path2file, variables):
    """
    Returns a dictionary with the specified variables for one sounding. 
    """
    profiles = {} 
    coords = ['pressure', 'latitude', 'longitude', 'flight_time']
    attrs = ['time_of_launch_HHmmss', 'date_YYYYMMDD', 'platform_name']
    filename = os.path.split(path2file)[1]
    platform = filename[0:3]
    va = variables.copy()
    ds = xarray.open_dataset(path2file)
    for variable in va:
        if variable in coords:
            profiles[variable] = ds.coords[variable].data[0]
        elif variable in attrs:
            profiles[variable] = ds.attrs[variable]
        else:
            profiles[variable] = ds.variables[variable].data[0]
    return profiles

def interpolate_profiles(profiles, variables, height):
    """
    Returns a dictionary containing the specified variables, interpolated to a given height vector. 
    """
    profiles_interp = copy.deepcopy(profiles)
    for i, profile in enumerate(profiles):
        for var in variables:
            var_interp = interp1d(profile['altitude'], profile[var], fill_value='extrapolate', bounds_error=False)(height)
            profiles_interp[i][var] = var_interp
        profiles_interp[i]['altitude'] = height
    
    return profiles_interp

def calc_potential_temperature(temperature, pressure):
    """
    Returns potential temperature [K].
    """
    R = typhon.constants.gas_constant_dry_air
    cp = typhon.constants.isobaric_mass_heat_capacity
    exponent = R / cp
    pot_temp = temperature * (1e5 / pressure) ** exponent
    
    return pot_temp

def calc_virtual_potential_temperature(temperature, mixing_ratio, pressure):
    """
    Returns virtual potential temperature [K]
    """
    potential_temp = calc_potential_temperature(temperature, pressure)
    virtual_potential_temp = potential_temp * (1 + 0.61 * mixing_ratio)
    
    return virtual_potential_temp

def calc_lower_tropospheric_stability(temperature, pressure):
    """
    Returns lower tropospheric stability (LTS) [K] calculated following Wood and Hartmann (2006).
    """
    if pressure[-1] < pressure[0]:
        surf_ind = 0
    else:
        surf_ind = -1
        
    potential_temperature = calc_potential_temperature(temperature, pressure)
    pot_temp_interp = interp1d(pressure, potential_temperature, fill_value='extrapolate', bounds_error=False)
    lts = pot_temp_interp(700e2) - potential_temperature[surf_ind]
    return lts

def calc_integrated_water_vapor(vmr, temperature, pressure, altitude):
    """
    Returns integrated water vapour [kg/m**-3]
    """
    nan_mask = np.isnan(vmr) + np.isnan(temperature) + np.isnan(pressure) + np.isnan(altitude)
    nan_ind = np.where(~nan_mask)
    
    if pressure[-1] > pressure[0]:
        vmr = np.flipud(vmr[nan_ind])
        temp = np.flipud(temperature[nan_ind])
        pres = np.flipud(pressure[nan_ind])
        alt = np.flipud(altitude[nan_ind])
    else:
        vmr = vmr[nan_ind]
        temp = temperature[nan_ind]
        pres = pressure[nan_ind]
        alt = altitude[nan_ind]
    
    iwv = typhon.physics.atmosphere.integrate_water_vapor(
        vmr=vmr, 
        p=pres, 
        T=temp, 
        z=alt
    )
    
    return iwv

def calc_precipitable_water(vmr, temperature, pressure, altitude):
    """
    Returns precipitable water [mm]
    """
    water_density = 997
    iwv = calc_integrated_water_vapor(vmr, temperature, pressure, altitude)
    pw = iwv / water_density * 1e3
    
    return pw

def calc_lcl_simple(relative_humidity, altitude):
    """
    Returns lifting condensation level (LCL) [m] calculated using a rough rule of thumb.
    """
    lcl_simple = altitude[0] + (100 - relative_humidity[0] * 1e2) * 25
    return lcl_simple

def calc_lcl(temperature, relative_humidity, altitude, levels=0):
    """
    Retruns lifting condensation level (LCL) [m] calculated according to Bolton (1980).
    """
    z0 = altitude[levels]
    cp = typhon.constants.isobaric_mass_heat_capacity
    g = typhon.constants.g
    tlcl = 1 / ((1 / (temperature[levels] - 55)) - (np.log(relative_humidity[levels]) / 2840.)) + 55
    zlcl = z0 - (cp * (tlcl - temperature[levels]) / g)
    mean_zlcl = np.mean(zlcl)
    
    return mean_zlcl

def calc_lcl_bco(temperature, dew_point, altitude, levels=0):
    """
    Returns lifting condensation level (LCL) [m] calculated as in the BCO quicklooks.
    """
    zlcl = ((temperature[levels] - dew_point[levels]) / 8. + altitude[levels] / 1000.) * 1000.
    mean_zlcl = np.mean(zlcl)
    return mean_zlcl

def calc_eis(pressure, temperature, relative_humidity, altitude):
    """
    Returns estimated inversion strength (EIS) [K] calculated according to Wood and Bretherton (2006).
    """
    moist_lapse_rate = typhon.physics.moist_lapse_rate(pressure, temperature)
    moist_lapse_rate_850hPa = interp1d(pressure, moist_lapse_rate)(850e2)
    
    lts = calc_lower_tropospheric_stability(temperature, pressure)
    lcl = calc_lcl(temperature, relative_humidity, altitude)
    z700hPa = interp1d(pressure, altitude)(700e2)
    eis = lts - moist_lapse_rate_850hPa * (z700hPa - lcl)
    
    return eis

def get_inversion_props(temperature, mixing_ratio, pressure, height, min_pres=700e2, virtual=True):
    """
    Returns inversion height defined as the maximum in the vertical gradient of virtual 
    potential temprature.
    """
    if virtual:
        potential_temperature = calc_virtual_potential_temperature(temperature, mixing_ratio, pressure)
    else:
        potential_temperature = calc_potential_temperature(temperature, pressure)
    # smooth curves
    gradient_pot_temp = np.diff(potential_temperature) / np.diff(height)
    gradient_pot_temp = np.append(gradient_pot_temp, gradient_pot_temp[-1])
    lower_tropo = np.logical_and(pressure > min_pres, pressure < 1000e2)
    inversion_ind = np.argmax(gradient_pot_temp[lower_tropo])
    inversion_height = height[lower_tropo][inversion_ind]
    
    return inversion_height, gradient_pot_temp

def smooth_profiles(profiles, variables, n_smooth):
    """ 
    Returns dictionary containing smoothed profiles of specified variables.
    """
    
    profiles_smoothed = copy.deepcopy(profiles)
    for i, profile in enumerate(profiles):
        for var in variables:
            var_smoothed = smooth_profile(profile[var], n_smooth)
            profiles_smoothed[i][var] = var_smoothed
    return profiles_smoothed
    
def smooth_profile(profile, box_pts):
    box = np.ones(box_pts)/box_pts
    y_smooth = np.convolve(profile, box, mode='valid')
    return y_smooth

def calc_mean_profiles(profiles, heights, variables):
    """ Returns mean profiles of specified variables from different soundings.
    Only works with profiles that are already interpolated to a common vertical grid.
    """
    mean_profiles = {}
    for var in variables:
        var_arr = np.zeros((len(profiles), len(heights)))
        for i in range(len(profiles)):
            var_arr[i] = profiles[i][var]
        var_mean = np.mean(var_arr, axis=0)
        mean_profiles[var] = var_mean
    return mean_profiles
    
def layer_mean_speed(wind_speed, altitude, z):
    """
    Returns mean over layer 250 m around given height 
    over given list of regridded wind profiles.
    """
    layer_ind = np.logical_and(altitude > z - 125, altitude < z + 150)
    layer_depth = np.diff(altitude)
    wind_speed_mean = np.sum(wind_speed[layer_ind] * layer_depth[layer_ind[:-1]]) / np.sum(layer_depth[layer_ind[:-1]])
    
    return wind_speed_mean

def layer_mean_dir(wind_direction, wind_speed, altitude, z):
    """
    Returns mean over layer 250 m around given height 
    """
    u, v = calc_wind_components(wind_direction, wind_speed)
    layer_ind = np.logical_and(altitude > z - 125, altitude < z + 150)
    layer_depth = np.diff(altitude)
    u_mean = np.sum(u[layer_ind] * layer_depth[layer_ind[:-1]]) / np.sum(layer_depth[layer_ind[:-1]])
    v_mean = np.sum(v[layer_ind] * layer_depth[layer_ind[:-1]]) / np.sum(layer_depth[layer_ind[:-1]])
    
    return calc_wind_direction(u_mean, v_mean)

def calc_wind_components(wind_direction, wind_speed):
    """ Returns wind components U and V.
    """
    u = -np.sin(wind_direction * np.pi / 180.) * wind_speed
    v = -np.cos(wind_direction * np.pi / 180.) * wind_speed
    
    return u, v

def calc_wind_direction(u, v):
    """ Returns wind direction.
    """
    direction = (180. + np.arctan2(u, v) * 180. / np.pi)
    return direction

def calc_wind_speed(u, v):
    """ Returns wind speed.
    """
    speed = np.sqrt(u**2 + v**2)
    return speed
