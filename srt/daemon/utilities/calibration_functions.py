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

def get_averaged_spectrum(fits_file):
    """
    open fits file and average all included spectra together
    """
    spectrum_file = fits.open(fits_file)
    average_spectrum = np.zeros(np.shape(spectrum_file[0].data),dtype=np.float64)
    #print(f'spectrum data shape: {np.shape(spectrum_file[0].data)}')
    #print(f'spectrum file shape: {np.shape(average_spectrum)}')

    num_spectra = len(spectrum_file)
    for i in range(0,num_spectra):
        spectrum=spectrum_file[i]
        average_spectrum += spectrum.data

    average_spectrum /= num_spectra

    return average_spectrum


def basic_cold_sky_calibration_fit(cold_sky_reference_filepath, t_sys=300, t_cal=300, num_channels=1, polynomial_order=20):
    """
    very basic calibration for single point temperature reference measurement. 
    calculates a polynomial fit for the spectrum and appropriately normalizes it
    only accounts for amplitude and assumes noise covariance between channels is zero for reference observation
    """

    average_cold_sky_spectrum = get_averaged_spectrum(cold_sky_reference_filepath) 
    average_cold_sky_spectrum_real = average_cold_sky_spectrum[:,:,:,0] #only get real values for now
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
    full_average_gains = np.ones((num_channels,num_channels))

    #infer gain corrections for the off-diagonal terms

    for i in range(num_channels):
        for j in range(num_channels):
            if i==j:
                full_normalized_spectra[i,j] = normalized_gain_spectrum[i]
                full_average_gains[i,j] = average_gain_correction[i]

            else: #gain correction is just product of corrections of the path gains
                full_normalized_spectra[i,j] = np.sqrt(normalized_gain_spectrum[i]*normalized_gain_spectrum[j])
                full_average_gains[i,j] = np.sqrt(average_gain_correction[i]*average_gain_correction[j])


    return full_normalized_spectra.reshape(num_channels**2,-1), full_average_gains.reshape(num_channels**2)

    
def additive_noise_calibration_fit(cold_sky_reference_filepath, calibrator_reference_filepath, t_sys=300, t_cal=300, num_channels=1, polynomial_order=20):

    """calibration using injected noise calibrator added to background signal
    only accounts for amplitude and assumes noise covariance between channels is zero for now
    """

    average_cold_sky_spectrum = np.abs(get_averaged_spectrum(cold_sky_reference_filepath))
    average_cold_sky_spectrum_real = average_cold_sky_spectrum[:,:,:,0] #only get real values for now
    average_calibrator_plus_sky_spectrum = np.abs(get_averaged_spectrum(calibrator_reference_filepath))
    average_calibrator_plus_sky_spectrum_real = average_calibrator_plus_sky_spectrum[:,:,:,0] #only get real values for now
    average_calibrator_spectrum = average_calibrator_plus_sky_spectrum_real - average_cold_sky_spectrum_real

    smoothed_calibrator_spectrum =np.ones_like(average_calibrator_spectrum[0]) #collapse to a single dimension corresponding to the diagonal

    relative_freq_values = np.linspace(-1, 1, np.shape(average_cold_sky_spectrum)[2])

    #compute corrections on the diagonal

    for i in range(num_channels):
        polynomial_fit = poly.Polynomial.fit(relative_freq_values, average_calibrator_spectrum[i,i], polynomial_order,)
        smoothed_calibrator_spectrum[i] = polynomial_fit(relative_freq_values)

    average_value = np.mean(smoothed_calibrator_spectrum,axis=1)
    #following nonsense is needed because of handling pultiple channels
    normalized_gain_spectrum = smoothed_calibrator_spectrum/(average_value*np.ones_like(smoothed_calibrator_spectrum).transpose()).transpose()
    average_gain_correction = average_value/t_cal

    full_normalized_spectra = np.ones_like(average_cold_sky_spectrum)
    full_average_gains = np.ones((num_channels,num_channels))

    #infer gain corrections for the off-diagonal terms

    for i in range(num_channels):
        for j in range(num_channels):
            if i==j:
                full_normalized_spectra[i,j] = normalized_gain_spectrum[i]
                full_average_gains[i,j] = average_gain_correction[i]

            else: #gain correction is just product of corrections of the path gains
                full_normalized_spectra[i,j] = np.sqrt(normalized_gain_spectrum[i]*normalized_gain_spectrum[j])
                full_average_gains[i,j] = np.sqrt(average_gain_correction[i]*average_gain_correction[j])


    return full_normalized_spectra.reshape(num_channels**2,-1), full_average_gains.reshape(num_channels**2)







