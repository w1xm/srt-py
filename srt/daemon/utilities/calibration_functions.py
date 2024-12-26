"""calibration_functions.py

dsheen 2024/11/30

Now with covariance matrix fits files

Calibration Mathhematics for new SRT calibration scheme

reimplements basic cold sky cal implemented directly in gnuradio 
and adds additional capabilities for more advanced telescopes.
"""

import numpy as np
import numpy.polynomial.polynomial as poly
from astropy.io import fits

def calibration_command_parameters(cal_type, num_channels=1, cal_duration=10, valid_states=[0,1]):
    '''
    return basic control parameters fromm the daemon that are not easily executed from this file
    if you have a system with anything interesting going on you'll probably need to modify this.

    Inputs
    ------

    num_channels : integer
        number of radio channels in use
    cal_type : string
        type of calibration to perform
    cal_duration : integer
        nebulously corresponds to the number of integration cycles the cal sequence should run for
    valid_states :
        list of valid calibrator states
        assumes valid_states[0] is all off and valid_states[-1], for middle states likely need to customize for a specific telescope

    Returns
    -------

    wait_cycles : list of integers
        list of calibrator periods for which to wait in between antenna commands
    cal_states : list of integers
        calibrator state to set prior to each wait period 
    '''

    if cal_type=='COLD_SKY':
        #no active calibrators involved
        #just take a single measurement at the position the telescope is pointed at
        wait_cycles = [cal_duration]
        cal_states = [0]

    elif cal_type=='NOISE_DIODE':
        #sequence on off measurements of active noise calibrator without phase calibration
        #assumes crosscoupling is low enough to not matter much but for 2 channel cal also mostly cancels it out anyway (only like 0.2K error on W1XMBIGDISH regardless)
        wait_cycles = [cal_duration]*4

        if num_channels==2: #there's a cute crosstalk cancelling sequence to be had for this case
            cal_states = [0,1,2,3] #half overlap calibrator pulses
        else:
            cal_states = []
            for i in len(wait_cycles):
                cal_states.append(valid_states[0] if i%2==0 else valid_states[-1])

        #pad the end to account for calibrator control latency (unavoidable due to integration time)
        wait_cycles.append(3)
        cal_states.append(cal_states[-1])

    elif cal_type=='REFL_PHASE':
        #ONLY VALID for dual pol feeds, Hyperspecific to implementation
        #no attempt has been made to generalize this for different channel counts
        cycle_time = 3 #just lock it in at this speed
        num_cycles = max(int(cal_duration/(cycle_time*3)),3)

        wait_cycles=[cal_duration]*num_cycles
        cal_states =[0,1,2]*int(num_cycles/3)

        #pad the end to account for calibrator control latency (unavoidable due to integration time)
        wait_cycles.append(3)
        cal_states.append(cal_states[-1])

    else:
        raise ValueError(f"Bad cal_type: {cal_type} is not a recognized calibration type")

    return wait_cycles, cal_states



def get_averaged_spectrum(fits_file):
    """
    open fits file, reconstruct complex array, and average all included spectra together
    
    returns the average covariance matrix as a complex array (or a single spectrum if using only one channel)
    """
    spectrum_file = fits.open(fits_file)
    
    average_spectrum = np.zeros(np.shape(spectrum_file[0].data)[0:-1],dtype=np.complex64) 
    #use array of same dimesions as fits file except for last axis (real vs complex) 
    #print(f'spectrum data shape: {np.shape(spectrum_file[0].data)}')
    #print(f'spectrum file shape: {np.shape(average_spectrum)}')

    num_spectra = len(spectrum_file)
    
    #average all integration periods together
    for i in range(0,num_spectra):
        spectrum=spectrum_file[i]
        average_spectrum += spectrum.data[:,:,:,0]+1j*spectrum.data[:,:,:,1]

    average_spectrum /= num_spectra

    return average_spectrum


def get_fits_data(fits_file):
    """
    open fits file and return numpy array of data plus list of metadata

    Inputs
    ------
    fits_file : path to fits file

    Returns
    -------
    data : numpy array dtype=complex64
        data from fits file as complex numpy array
    metadata : list
        list of metadata from fits file
    """

    hdul = fits.open(fits_file)

    #get metadata and rf data out into a more useful form for my purpose
    fits_data = []
    fits_metadata = []

    for i, hdu in enumerate(hdul):
        fits_metadata.append(json.loads(hdu.header["METADATA"]))
        fits_data.append(np.array(hdu.data[:,:,:,0] +1j*hdu.data[:,:,:,1]))

    hdul.close()

    #make data into nice big numpy array
    fits_data=np.array(fits_data)

    return fits_data, fits_metadata



def calculate_calibration_corrections(ref_file, cal_type, tsys=np.array([300]), tref=np.array([300]), num_channels=1, valid_states=range(len(2))):

    """
    takes in a bunch of parameters plus a fits file with recorded calibration data and returns calibration corrections for the telescope. 
    accuracy and complexity dependent on calibration type

    Returns
    -------

    cal_values : numpy array
        complex calibration correctionn matrix to be applied to data coming out of the radio
    """

    #internal variables
    polynomial_order=20

    #start by pulling in fits file data and metadata
    fits_data, fits_metadata = get_fits_data(fits_file)
    #and create a reference axis for fitting data
    relative_freq_values = np.linspace(-1, 1, len(fits_data[0,0,0]))

    if cal_type=="COLD_SKY":

        #just average across the whole data set and try to correct for estimated total temperature
        average_spectra = np.mean(fits_data,axis=0)
        amplitude_correction_mat = np.ones_like(average_spectra)

        #compute corections for diagonal of covariance matrix

        for i in range(num_channels):
            poly_fit = poly.Polynomial.fit(relative_freq_values, np.real(average_spectra[i,i]), polynomial_order)
            amplitude_correction_mat[i,i] = (tref[i]+tsys[i])/poly_fit(relative_freq_values)

        #compute corections for diagonal of covariance matrix

    elif cal_type=="NOISE_DIODE":

    elif cal_type=='REFL_PHASE':

    else:
        raise ValueError(f"Bad cal_type: {cal_type} is not a recognized calibration type")

    return cal_values






def basic_cold_sky_calibration_fit(cold_sky_reference_filepath, t_sys=np.array([300]), t_cal=np.array([300]), num_channels=1, polynomial_order=20):
    """
    very basic calibration for single point temperature reference measurement. 
    calculates a polynomial fit for the spectrum and appropriately normalizes it
    only accounts for amplitude and assumes noise covariance between channels is zero for reference observation
    """

    average_cold_sky_spectrum = get_averaged_spectrum(cold_sky_reference_filepath) 
    #average_cold_sky_spectrum_real = average_cold_sky_spectrum[:,:,:,0] #only get real values for now
    average_cold_sky_spectrum_real = np.real(average_cold_sky_spectrum)
    relative_freq_values = np.linspace(-1, 1, np.shape(average_cold_sky_spectrum_real)[2])

    smoothed_cold_sky_spectrum =np.ones_like(average_cold_sky_spectrum_real[0]) #drop a dimension to save my sanity here. 
    #this is ONLY calculating the diagonal elements


    for i in range(num_channels):
        polynomial_fit = poly.Polynomial.fit(relative_freq_values, average_cold_sky_spectrum_real[i,i], polynomial_order,)
        smoothed_cold_sky_spectrum[i] = polynomial_fit(relative_freq_values)

    #calculate gain corrections for the diagonal terms

    average_value = np.mean(smoothed_cold_sky_spectrum, axis=1)
    #following nonsense is needed because of handling pultiple channels
    normalized_gain_spectrum = smoothed_cold_sky_spectrum/(average_value*np.ones_like(smoothed_cold_sky_spectrum).transpose()).transpose()
    average_gain_correction = average_value/(t_sys+t_cal)

    full_normalized_spectra = np.ones_like(average_cold_sky_spectrum)
    cal_coefficients = np.ones_like(average_cold_sky_spectrum)
    full_average_gains = np.ones((num_channels,num_channels))

    #infer gain corrections for the off-diagonal terms from the individual channel gains

    for i in range(num_channels):
        for j in range(num_channels):
            if i==j:
                full_normalized_spectra[i,j] = normalized_gain_spectrum[i]
                full_average_gains[i,j] = average_gain_correction[i]
                cal_coefficients[i,j] = 1.0/(full_normalized_spectra[i,j]*full_average_gains[i,j])

            else: #gain correction is just product of corrections of the path gains
                full_normalized_spectra[i,j] = np.sqrt(normalized_gain_spectrum[i]*normalized_gain_spectrum[j]) 
                full_average_gains[i,j] = np.sqrt(average_gain_correction[i]*average_gain_correction[j])
                cal_coefficients[i,j] = 1.0/(full_normalized_spectra[i,j]*full_average_gains[i,j])

    return cal_coefficients.reshape(num_channels**2,-1), full_average_gains.reshape(num_channels**2)

    
def additive_noise_calibration_fit(cold_sky_reference_filepath, calibrator_reference_filepath, t_sys=np.array([300]), t_cal=np.array([300]), num_channels=1, polynomial_order=20):

    """calibration using injected noise calibrator added to background signal
    only accounts for amplitude and assumes noise covariance between channels is zero for now
    """

    average_cold_sky_spectrum = get_averaged_spectrum(cold_sky_reference_filepath)
    #average_cold_sky_spectrum_real = average_cold_sky_spectrum[:,:,:,0] #only get real values for now
    average_cold_sky_spectrum_real = np.real(average_cold_sky_spectrum)
    average_calibrator_plus_sky_spectrum = get_averaged_spectrum(calibrator_reference_filepath)
    #average_calibrator_plus_sky_spectrum_real = average_calibrator_plus_sky_spectrum[:,:,:,0] #only get real values for now
    average_calibrator_plus_sky_spectrum_real = np.real(average_calibrator_plus_sky_spectrum)
    average_calibrator_spectrum = average_calibrator_plus_sky_spectrum_real - average_cold_sky_spectrum_real
    

    smoothed_calibrator_spectrum =np.ones_like(average_calibrator_spectrum[0]) #collapse to a single dimension corresponding to the diagonal

    relative_freq_values = np.linspace(-1, 1, np.shape(average_cold_sky_spectrum)[2])
    #print(f'calibrator spectrum shape {np.shape(average_calibrator_spectrum)}')

    #compute corrections on the diagonal

    for i in range(num_channels):
        polynomial_fit = poly.Polynomial.fit(relative_freq_values, average_calibrator_spectrum[i,i], polynomial_order,)
        smoothed_calibrator_spectrum[i] = polynomial_fit(relative_freq_values)

    average_value = np.mean(smoothed_calibrator_spectrum,axis=1)
    #following nonsense is needed because of handling pultiple channels
    normalized_gain_spectrum = smoothed_calibrator_spectrum/(average_value*np.ones_like(smoothed_calibrator_spectrum).transpose()).transpose()
    average_gain_correction = average_value/t_cal

    full_normalized_spectra = np.ones_like(average_cold_sky_spectrum)
    cal_coefficients = np.ones_like(average_cold_sky_spectrum)
    full_average_gains = np.ones((num_channels,num_channels))

    #infer gain corrections for the off-diagonal terms from the individual channel gains

    for i in range(num_channels):
        for j in range(num_channels):
            if i==j:
                full_normalized_spectra[i,j] = normalized_gain_spectrum[i]
                full_average_gains[i,j] = average_gain_correction[i]
                cal_coefficients[i,j] = 1.0/(full_normalized_spectra[i,j]*full_average_gains[i,j])

            else: #gain correction is just product of corrections of the path gains
                full_normalized_spectra[i,j] = np.sqrt(normalized_gain_spectrum[i]*normalized_gain_spectrum[j])
                full_average_gains[i,j] = np.sqrt(average_gain_correction[i]*average_gain_correction[j])
                cal_coefficients[i,j] = 1.0/(full_normalized_spectra[i,j]*full_average_gains[i,j])
    
    return cal_coefficients.reshape(num_channels**2,-1), full_average_gains.reshape(num_channels**2)







