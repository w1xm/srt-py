"""calibration_functions.py

dsheen 2024/11/30

Now with covariance matrix fits files

Calibration Mathhematics for new SRT calibration scheme

reimplements basic cold sky cal implemented directly in gnuradio 
and adds additional capabilities for more advanced telescopes.
"""

import numpy as np
import numpy.polynomial.polynomial as poly
import scipy.stats as stats
import json
from astropy.io import fits

def calibration_command_parameters(cal_type, num_channels=1, cal_duration=10, valid_masks=[0,1]):
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
    valid_masks :
        list of valid calibrator masks that cal_on can take, the index into this is what actually needs to be sent back to the daemon
        assumes valid_masks[0] is all off and valid_masks[-1], for middle states likely need to customize for a specific telescope

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
                cal_states.append(0 if i%2==0 else len(valid_masks)-1)

        #pad the end to account for calibrator control latency (unavoidable due to integration time)
        wait_cycles.append(3)
        cal_states.append(cal_states[-1])

    elif cal_type=='REFL_PHASE':
        #ONLY VALID for dual pol feeds, Hyperspecific to implementation
        #no attempt has been made to generalize this for different channel counts
        if num_channels !=2:
            raise ValueError(f"cal_type: 'REFL_PHASE' only valid for dual pol feeds")

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



def calculate_calibration_corrections(ref_file, cal_type, tsys=np.array([300]), tref=np.array([300]), num_channels=1, valid_masks=range(2)):

    """
    takes in a bunch of parameters plus a fits file with recorded calibration data and returns calibration corrections for the telescope. 
    accuracy and complexity dependent on calibration type

    Returns
    -------

    correction_mat : numpy array
        complex calibration correctionn matrix to be applied to data coming out of the radio
    """

    #internal variables
    polynomial_order=20

    #start by pulling in fits file data and metadata
    fits_data, fits_metadata = get_fits_data(ref_file)
    #and create a reference axis for fitting data
    relative_freq_values = np.linspace(-1, 1, len(fits_data[0,0,0]))

    if cal_type=="COLD_SKY":

        #just average across the whole data set and try to correct for estimated total temperature
        average_spectra = np.mean(fits_data,axis=0)
        correction_mat = np.ones_like(average_spectra)

        #compute corections for diagonal of covariance matrix

        for i in range(num_channels):
            poly_fit = poly.Polynomial.fit(relative_freq_values, np.real(average_spectra[i,i]), polynomial_order)
            correction_mat[i,i] = (tref[i]+tsys[i])/poly_fit(relative_freq_values)

        #propagate to off-diagonal terms in the matrix

        for i in range(num_channels):
            for j in range(num_channels):
                if i!=j:
                    correction_mat[i,j] = np.sqrt(correction_mat[i,i]*correction_mat[j,j])

    elif cal_type=="NOISE_DIODE":

        #select out and average specific calibration states based on metadata
        state_averages = np.zeros((len(valid_masks),np.shape(fits_data[0])[0],np.shape(fits_data[0])[1],np.shape(fits_data[0])[2]),dtype=np.complex64)
        correction_mat = np.ones_like(fits_data[0])

        for i, cal_state in enumerate(valid_masks):
            count = 0
            temp_array = np.zeros_like(fits_data[0])

            for j in range(len(fits_metadata)):
                if fits_metadata[j]["cal_on"] == cal_state:
                    temp_array = temp_array + fits_data[j]
                    count +=1

            try:
                state_averages[i] = temp_array/count
            except:
                state_averages[i] = state_averages[i]

        #baseline subtraction for difference estimation

        if num_channels==2: #handle the clever overlap from above
            calibrator_diag_spectra = [np.mean(np.array([state_averages[1,0,0],state_averages[3,0,0]]),axis=0)-np.mean(np.array([state_averages[0,0,0],state_averages[2,0,0]]),axis=0),
                                        np.mean(np.array([state_averages[2,1,1],state_averages[3,1,1]]),axis=0)-np.mean(np.array([state_averages[0,1,1],state_averages[1,1,1]]),axis=0)]
        else:
            calibrator_diag_spectra = [state_averages[-1,i,i]-state_averages[0,i,i] for i in range(num_channels)] #assume first is off state and last is all on

        #compute corections for diagonal of covariance matrix

        for i in range(num_channels):
            poly_fit = poly.Polynomial.fit(relative_freq_values, np.real(calibrator_diag_spectra[i]), polynomial_order)
            correction_mat[i,i] = tref[i]/poly_fit(relative_freq_values)

        #propagate to off-diagonal terms in the matrix

        for i in range(num_channels):
            for j in range(num_channels):
                if i!=j:
                    correction_mat[i,j] = np.sqrt(correction_mat[i,i]*correction_mat[j,j])

    elif cal_type=='REFL_PHASE':
        #just assume this needs to be customized for a given telescope
        if num_channels !=2:
            raise ValueError(f"cal_type: 'REFL_PHASE' only valid for dual pol feeds")

        #select out and average specific calibration states based on metadata
        state_averages = np.zeros((len(valid_masks),np.shape(fits_data[0])[0],np.shape(fits_data[0])[1],np.shape(fits_data[0])[2]),dtype=np.complex64)
        amplitude_correction_mat = np.ones_like(fits_data[0])
        phase_correction_mat = np.ones_like(fits_data[0])

        for i, cal_state in enumerate(valid_masks):
            count = 0
            temp_array = np.zeros_like(fits_data[0])

            for j in range(len(fits_metadata)):
                if fits_metadata[j]["cal_on"] == cal_state:
                    temp_array = temp_array + fits_data[j]
                    count +=1
            try:
                state_averages[i] = temp_array/count
            except:
                state_averages[i] = state_averages[i]

        #baseline subtraction HARD CODED FOR WR66

        cal_1_subtracted = state_averages[1] - state_averages[0]
        cal_2_subtracted = state_averages[2] - state_averages[0]

        #values to feed into amplitude cal matrix

        diag_spectra = [cal_1_subtracted[0,0],cal_1_subtracted[1,1]]

        #compute diagonal of amplitude correction matrix

        for i in range(num_channels):
            poly_fit = poly.Polynomial.fit(relative_freq_values, np.real(diag_spectra[i]), polynomial_order)
            amplitude_correction_mat[i,i] = tref[i]/poly_fit(relative_freq_values)

        #propagate to off-diagonal terms in the matrix

        for i in range(num_channels):
            for j in range(num_channels):
                if i!=j:
                    amplitude_correction_mat[i,j] = np.sqrt(amplitude_correction_mat[i,i]*amplitude_correction_mat[j,j])

        #Phase Cal
        cal_1_phase = np.unwrap(np.angle(cal_1_subtracted))
        cal_2_phase = np.unwrap(np.angle(cal_2_subtracted))
        phase_error = 0.5*(cal_1_phase[0,1] + cal_2_phase[0,1])

        phasefit = stats.linregress(relative_freq_values,phase_error)
        print(f'r value = {phasefit.rvalue}')
        print(f'p value = {phasefit.pvalue}')

        fitphase = phasefit.intercept*np.ones_like(relative_freq_values) + phasefit.slope*relative_freq_values

        phase_correction_mat[0,1] = np.exp(-1j*fitphase)
        phase_correction_mat[1,0] = np.exp(1j*fitphase)

        correction_mat = amplitude_correction_mat * phase_correction_mat

    else:
        raise ValueError(f"Bad cal_type: {cal_type} is not a recognized calibration type")

    average_gains = 1/np.mean(np.abs(correction_mat),axis=2)

    return correction_mat.reshape(num_channels**2,-1), average_gains.reshape(num_channels**2)