"""
Block to save covariance spectra to a fits file
"""

import numpy as np
from gnuradio import gr
import pmt
import json

import pathlib
from datetime import datetime, timezone
from astropy.io import fits


class blk(gr.sync_block):
    """Embedded Python Block - Saving """

    def __init__(
        self, directory=".", filename="test.fits", spectrum_len=4096, num_channels=2
    ):  # only default arguments here
        """arguments to this function show up as parameters in GRC"""
        gr.sync_block.__init__(
            self,
            name="Embedded Python Block",  # will show up in GRC
            in_sig=[(np.complex64, spectrum_len * num_channels**2)] ,
            #in_sig=[(np.float32, spectrum_len) for i in range(num_channels)],
            out_sig=None,
        )
        # if an attribute with the same name as a parameter is found,
        # a callback is registered (properties work, too).
        self.directory = directory
        self.filename = filename
        self.spectrum_len = spectrum_len
        self.num_channels = num_channels

    def work(self, input_items, output_items):
        """Saving Spectrum Data to a FITS File"""
        # we're just going to assume the inputs are the same length because they will be, and for now assume they share the same metadata 
        #not too worried babout getting this perfect because I'll need to rewrite this later anyway
        file_path = pathlib.Path(self.directory, self.filename)

        with open(file_path, "ab+") as file:
            for input_array in input_items[0]:

                tags_0 = self.get_tags_in_window(0, 0, len(input_items[0]))
                #tags_1 = self.get_tags_in_window(0, 0, len(input_items[1]))
                tags_dict_0 = {pmt.to_python(tag.key): pmt.to_python(tag.value) for tag in tags_0}
                #tags_dict_1 = {pmt.to_python(tag.key): pmt.to_python(tag.value) for tag in tags_1}

                time_since_epoch = tags_dict_0["rx_time"][0] + tags_dict_0["rx_time"][1]
                date = datetime.fromtimestamp(time_since_epoch, timezone.utc)
                metadata = tags_dict_0["metadata"]
                samp_rate = metadata["samp_rate"]
                num_integrations = metadata["num_integrations"]
                freq = metadata["freq"]
                num_bins = metadata["num_bins"]
                soutrack = metadata["soutrack"]

                hdr = fits.Header()
                hdr["BUNIT"] = "K"
                hdr["CTYPE4"] = "Channel 0"
                hdr["CTYPE3"] = "Channel 1"
                hdr["CTYPE2"] = "Freq"
                hdr["CTYPE1"] = "COMPLEX"
                hdr["CRPIX4"] = 0 #referenced to channel 0 at coordinate 0
                hdr["CRPIX3"] = 0 #referenced to channel 0 at coordinate 0
                hdr["CRPIX2"] = num_bins / float(2)  # Reference pixel (center)
                hdr["CRPIX1"] = 0
                hdr["CRVAL4"] = 0 #channel 0
                hdr["CRVAL3"] = 0 #channel 0
                hdr["CRVAL2"] = freq  # Center, USRP, frequency
                hdr["CRVAL1"] = 0
                hdr["CDELT4"] = 1
                hdr["CDELT3"] = 1
                hdr["CDELT2"] = samp_rate / (1 * num_bins)  # Channel width
                hdr["CDELT1"] = 1
                hdr["CUNIT4"] = "Radio Channel"
                hdr["CUNIT3"] = "Radio Channel"
                hdr["CUNIT2"] = "Hz"

                #hdr["TELESCOP"] = "SmallRadioTelescope"
                hdr["TELESCOP"] = "MediumRadioTelescope"
                hdr["OBJECT"] = soutrack
                hdr["OBSTIME"] = (num_bins * num_integrations) / samp_rate

                hdr["DATE-OBS"] = date.strftime("%Y-%m-%d")
                hdr["UTC"] = date.strftime("%H:%M:00%s")
                hdr["METADATA"] = json.dumps(metadata)

                #need to add an axis that separately includes real and complex parts of the data
                covariances = input_array.reshape(self.num_channels,self.num_channels,self.spectrum_len)
                float_array = np.moveaxis(np.array([np.real(covariances),np.imag(covariances)]),0,-1) #.swapaxes(0,3)

                #append neatly reshaped input containing covariance matrix data
                fits.append(file, float_array.reshape(self.num_channels,self.num_channels,self.spectrum_len,2), hdr) #need to explicitly reshape inline to force it to save array in correct shape
                #fits.append(file, combined_data, hdr) #append both spectra.
                #file.close()
                # p = np.sum(input_array)
                # a = len(input_array)
                # pwr = (tsys + tcal) * p / (a * calpwr)
                # ppwr = pwr - tsys
        return len(input_items[0])
