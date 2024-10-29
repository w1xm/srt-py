# -*- coding: utf-8 -*-

#
# SPDX-License-Identifier: GPL-3.0
#
# GNU Radio Python Flow Graph
# Title: filter and integration
# GNU Radio version: 3.10.10.0

from gnuradio import blocks
from gnuradio import fft
from gnuradio.fft import window
from gnuradio import filter
from gnuradio import gr
from gnuradio.filter import firdes
import sys
import signal
import math
import numpy as np







class filter_integrate(gr.hier_block2):
    def __init__(self, fft_window, num_bins=256, num_integrations=100000):
        gr.hier_block2.__init__(
            self, "filter and integration",
                gr.io_signature(1, 1, gr.sizeof_gr_complex*1),
                gr.io_signature(1, 1, gr.sizeof_float*num_bins),
        )

        ##################################################
        # Parameters
        ##################################################
        self.fft_window = fft_window
        self.num_bins = num_bins
        self.num_integrations = num_integrations

        ##################################################
        # Variables
        ##################################################
        self.sinc_sample_locations = sinc_sample_locations = np.arange(-np.pi*4/2.0, np.pi*4/2.0, np.pi/num_bins)
        self.sinc_samples = sinc_samples = np.sinc(sinc_sample_locations/np.pi)
        self.custom_window = custom_window = sinc_samples*np.hamming(4*num_bins)

        ##################################################
        # Blocks
        ##################################################

        self.fft_vxx_0 = fft.fft_vcc(num_bins, True, fft_window, True, 3)
        self.dc_blocker_xx_0 = filter.dc_blocker_cc((num_bins*num_integrations), False)
        self.blocks_stream_to_vector_0_2 = blocks.stream_to_vector(gr.sizeof_gr_complex*1, num_bins)
        self.blocks_stream_to_vector_0_1 = blocks.stream_to_vector(gr.sizeof_gr_complex*1, num_bins)
        self.blocks_stream_to_vector_0_0 = blocks.stream_to_vector(gr.sizeof_gr_complex*1, num_bins)
        self.blocks_stream_to_vector_0 = blocks.stream_to_vector(gr.sizeof_gr_complex*1, num_bins)
        self.blocks_skiphead_0 = blocks.skiphead(gr.sizeof_gr_complex*1, (num_bins*num_integrations))
        self.blocks_multiply_const_xx_0 = blocks.multiply_const_ff(1.0/float(num_integrations), num_bins)
        self.blocks_multiply_const_vxx_0_0_0_0 = blocks.multiply_const_vcc(custom_window[0:num_bins])
        self.blocks_multiply_const_vxx_0_0_0 = blocks.multiply_const_vcc(custom_window[num_bins:2*num_bins])
        self.blocks_multiply_const_vxx_0_0 = blocks.multiply_const_vcc(custom_window[2*num_bins:3*num_bins])
        self.blocks_multiply_const_vxx_0 = blocks.multiply_const_vcc(custom_window[-num_bins:])
        self.blocks_integrate_xx_0 = blocks.integrate_ff(num_integrations, num_bins)
        self.blocks_delay_0_1 = blocks.delay(gr.sizeof_gr_complex*1, num_bins)
        self.blocks_delay_0_0 = blocks.delay(gr.sizeof_gr_complex*1, (num_bins*2))
        self.blocks_delay_0 = blocks.delay(gr.sizeof_gr_complex*1, (num_bins*3))
        self.blocks_complex_to_mag_squared_0 = blocks.complex_to_mag_squared(num_bins)
        self.blocks_add_xx_0 = blocks.add_vcc(num_bins)


        ##################################################
        # Connections
        ##################################################
        self.connect((self.blocks_add_xx_0, 0), (self.fft_vxx_0, 0))
        self.connect((self.blocks_complex_to_mag_squared_0, 0), (self.blocks_integrate_xx_0, 0))
        self.connect((self.blocks_delay_0, 0), (self.blocks_stream_to_vector_0_2, 0))
        self.connect((self.blocks_delay_0_0, 0), (self.blocks_stream_to_vector_0_0, 0))
        self.connect((self.blocks_delay_0_1, 0), (self.blocks_stream_to_vector_0_1, 0))
        self.connect((self.blocks_integrate_xx_0, 0), (self.blocks_multiply_const_xx_0, 0))
        self.connect((self.blocks_multiply_const_vxx_0, 0), (self.blocks_add_xx_0, 0))
        self.connect((self.blocks_multiply_const_vxx_0_0, 0), (self.blocks_add_xx_0, 1))
        self.connect((self.blocks_multiply_const_vxx_0_0_0, 0), (self.blocks_add_xx_0, 2))
        self.connect((self.blocks_multiply_const_vxx_0_0_0_0, 0), (self.blocks_add_xx_0, 3))
        self.connect((self.blocks_multiply_const_xx_0, 0), (self, 0))
        self.connect((self.blocks_skiphead_0, 0), (self.blocks_delay_0, 0))
        self.connect((self.blocks_skiphead_0, 0), (self.blocks_delay_0_0, 0))
        self.connect((self.blocks_skiphead_0, 0), (self.blocks_delay_0_1, 0))
        self.connect((self.blocks_skiphead_0, 0), (self.blocks_stream_to_vector_0, 0))
        self.connect((self.blocks_stream_to_vector_0, 0), (self.blocks_multiply_const_vxx_0, 0))
        self.connect((self.blocks_stream_to_vector_0_0, 0), (self.blocks_multiply_const_vxx_0_0_0, 0))
        self.connect((self.blocks_stream_to_vector_0_1, 0), (self.blocks_multiply_const_vxx_0_0, 0))
        self.connect((self.blocks_stream_to_vector_0_2, 0), (self.blocks_multiply_const_vxx_0_0_0_0, 0))
        self.connect((self.dc_blocker_xx_0, 0), (self.blocks_skiphead_0, 0))
        self.connect((self.fft_vxx_0, 0), (self.blocks_complex_to_mag_squared_0, 0))
        self.connect((self, 0), (self.dc_blocker_xx_0, 0))


    def get_fft_window(self):
        return self.fft_window

    def set_fft_window(self, fft_window):
        self.fft_window = fft_window

    def get_num_bins(self):
        return self.num_bins

    def set_num_bins(self, num_bins):
        self.num_bins = num_bins
        self.set_custom_window(self.sinc_samples*np.hamming(4*self.num_bins))
        self.set_sinc_sample_locations(np.arange(-np.pi*4/2.0, np.pi*4/2.0, np.pi/self.num_bins))
        self.blocks_delay_0.set_dly(int((self.num_bins*3)))
        self.blocks_delay_0_0.set_dly(int((self.num_bins*2)))
        self.blocks_delay_0_1.set_dly(int(self.num_bins))
        self.blocks_multiply_const_vxx_0.set_k(self.custom_window[-self.num_bins:])
        self.blocks_multiply_const_vxx_0_0.set_k(self.custom_window[2*self.num_bins:3*self.num_bins])
        self.blocks_multiply_const_vxx_0_0_0.set_k(self.custom_window[self.num_bins:2*self.num_bins])
        self.blocks_multiply_const_vxx_0_0_0_0.set_k(self.custom_window[0:self.num_bins])

    def get_num_integrations(self):
        return self.num_integrations

    def set_num_integrations(self, num_integrations):
        self.num_integrations = num_integrations
        self.blocks_multiply_const_xx_0.set_k(1.0/float(self.num_integrations))

    def get_sinc_sample_locations(self):
        return self.sinc_sample_locations

    def set_sinc_sample_locations(self, sinc_sample_locations):
        self.sinc_sample_locations = sinc_sample_locations
        self.set_sinc_samples(np.sinc(self.sinc_sample_locations/np.pi))

    def get_sinc_samples(self):
        return self.sinc_samples

    def set_sinc_samples(self, sinc_samples):
        self.sinc_samples = sinc_samples
        self.set_custom_window(self.sinc_samples*np.hamming(4*self.num_bins))

    def get_custom_window(self):
        return self.custom_window

    def set_custom_window(self, custom_window):
        self.custom_window = custom_window
        self.blocks_multiply_const_vxx_0.set_k(self.custom_window[-self.num_bins:])
        self.blocks_multiply_const_vxx_0_0.set_k(self.custom_window[2*self.num_bins:3*self.num_bins])
        self.blocks_multiply_const_vxx_0_0_0.set_k(self.custom_window[self.num_bins:2*self.num_bins])
        self.blocks_multiply_const_vxx_0_0_0_0.set_k(self.custom_window[0:self.num_bins])

