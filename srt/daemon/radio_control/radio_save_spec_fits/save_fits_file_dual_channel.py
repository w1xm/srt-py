"""
Embedded Python Blocks:

Each time this file is saved, GRC will instantiate the first class it finds
to get ports and parameters of your block. The arguments to __init__  will
be the parameters. All of them are required to have default values!
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
        self, directory=".", filename="test.fits", vec_length=4096, num_channels=2
    ):  # only default arguments here
        """arguments to this function show up as parameters in GRC"""
        gr.sync_block.__init__(
            self,
            name="Embedded Python Block",  # will show up in GRC
            in_sig=[(np.float32, vec_length) for i in range(num_channels)],
            out_sig=None,
        )
        # if an attribute with the same name as a parameter is found,
        # a callback is registered (properties work, too).
        self.directory = directory
        self.filename = filename
        self.vec_length = vec_length

    def work(self, input_items, output_items):
        """Saving Spectrum Data to a FITS File"""
        # we're just going to assume the inputs are the same length because they will be, and for now assume they share the same metadata 
        #not too worried babout getting this perfect because I'll need to rewrite this later anyway
        file_path = pathlib.Path(self.directory, self.filename)
        #for i, input_array in enumerate(input_items[0]):
        for input_array_0, input_array_1 in zip(input_items[0],input_items[1]): #idk why enumerate was involved here. not needed
            file = open(file_path, "ab+")
            tags_0 = self.get_tags_in_window(0, 0, len(input_items[0]))
            tags_1 = self.get_tags_in_window(0, 0, len(input_items[1]))
            tags_dict_0 = {pmt.to_python(tag.key): pmt.to_python(tag.value) for tag in tags_0}
            tags_dict_1 = {pmt.to_python(tag.key): pmt.to_python(tag.value) for tag in tags_1}

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
            hdr["CTYPE1"] = "Freq"
            hdr["CRPIX1"] = num_bins / float(2)  # Reference pixel (center)
            hdr["CRVAL1"] = freq  # Center, USRP, frequency
            hdr["CDELT1"] = samp_rate / (1 * num_bins)  # Channel width
            hdr["CUNIT1"] = "Hz"

            hdr["TELESCOP"] = "SmallRadioTelescope"
            hdr["OBJECT"] = soutrack
            hdr["OBSTIME"] = (num_bins * num_integrations) / samp_rate

            hdr["DATE-OBS"] = date.strftime("%Y-%m-%d")
            hdr["UTC"] = date.strftime("%H:%M:00%s")
            hdr["METADATA"] = json.dumps(metadata)

            fits.append(file, [input_array_0,input_array_1], hdr) #append both spectra.
            file.close()
            # p = np.sum(input_array)
            # a = len(input_array)
            # pwr = (tsys + tcal) * p / (a * calpwr)
            # ppwr = pwr - tsys
        return len(input_items[0])
