#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#
# SPDX-License-Identifier: GPL-3.0
#
# GNU Radio Python Flow Graph
# Title: radio_process_dual_channel
# GNU Radio version: 3.10.9.2

import os
import sys
#sys.path.append(os.environ.get('GRC_HIER_PATH', os.path.expanduser('~/.grc_gnuradio')))


from gnuradio import blocks
import pmt
from gnuradio import gr
from gnuradio.filter import firdes
from gnuradio.fft import window
import signal
from argparse import ArgumentParser
from gnuradio.eng_arg import eng_float, intx
from gnuradio import eng_notation
from gnuradio import uhd
import time
from gnuradio import zeromq
from xmlrpc.server import SimpleXMLRPCServer
import threading
import math
import numpy as np
from . import add_clock_tags
from . import covariance_matrix
from . import weighted_overlap_fft  # grc-generated hier_block manually relocated to directory
#from . import calibrator_control_strobe
from . import calibrator_timestamping_block_with_full_time_synchronization as calibrator_control_strobe


class radio_process_dual_channel(gr.top_block):

    def __init__(self, num_bins=256, num_integrations=100000):
        gr.top_block.__init__(self, "radio_process_dual_channel", catch_exceptions=True)

        ##################################################
        # Parameters
        ##################################################
        self.num_bins = num_bins
        self.num_integrations = num_integrations
        self.num_channels = num_channels = 2

        ##################################################
        # Variables
        ##################################################
        self.sinc_sample_locations = sinc_sample_locations = np.arange(-np.pi*4/2.0, np.pi*4/2.0, np.pi/num_bins)
        self.sinc_samples = sinc_samples = np.sinc(sinc_sample_locations/np.pi)
        self.freq = freq = 1420000000
        self.vlsr = vlsr = np.nan
        self.tsys = tsys = np.array([171]*num_channels)
        self.tcal = tcal = np.array([290]*num_channels)
        self.tag_period = tag_period = num_bins*num_integrations
        self.soutrack = soutrack = "at_stow"
        self.samp_rate = samp_rate = 2000000
        self.rf_gain = rf_gain = 20
        self.rf_freq = rf_freq = freq
        self.motor_el = motor_el = np.nan
        self.motor_az = motor_az = np.nan
        self.is_running = is_running = False
        self.glon = glon = np.nan
        self.glat = glat = np.nan
        self.fft_window = fft_window = window.blackmanharris(num_bins)
        self.custom_window = custom_window = sinc_samples*np.hamming(4*num_bins)
        self.calibrator_mask = calibrator_mask = 0b000000000011
        self.cal_values_real = cal_values_real = np.ones((num_channels**2,num_bins))
        self.cal_values_imag = cal_values_imag = np.zeros((num_channels**2,num_bins))
        self.cal_values = cal_values = cal_values_real+1j*cal_values_imag
        self.cal_pwr = cal_pwr = np.array([1]*num_channels**2)
        self.cal_on = cal_on = 0
        self.beam_switch = beam_switch = 0

        ##################################################
        # Blocks
        ##################################################

        self.zeromq_pub_sink_2_0 = zeromq.pub_sink(gr.sizeof_gr_complex, (4*num_bins), 'tcp://127.0.0.1:5561', 100, False, (-1), '', True)
        self.zeromq_pub_sink_2 = zeromq.pub_sink(gr.sizeof_gr_complex, (4*num_bins), 'tcp://127.0.0.1:5560', 100, True, (-1), '', True)
        self.zeromq_pub_sink_1_0 = zeromq.pub_sink(gr.sizeof_gr_complex, (4*num_bins), 'tcp://127.0.0.1:5562', 100, True, (-1), '', True)
        self.zeromq_pub_sink_1 = zeromq.pub_sink(gr.sizeof_gr_complex, (4*num_bins), 'tcp://127.0.0.1:5563', 100, False, (-1), '', True)
        self.zeromq_pub_sink_0_0 = zeromq.pub_sink(gr.sizeof_gr_complex, 2, 'tcp://127.0.0.1:5559', 100, False, (-1), '', True)
        self.zeromq_pub_sink_0 = zeromq.pub_sink(gr.sizeof_gr_complex, 2, 'tcp://127.0.0.1:5558', 100, True, (-1), '', True)
        self.xmlrpc_server_0 = SimpleXMLRPCServer(('localhost', 5557), allow_none=True)
        self.xmlrpc_server_0.register_instance(self)
        self.xmlrpc_server_0_thread = threading.Thread(target=self.xmlrpc_server_0.serve_forever)
        self.xmlrpc_server_0_thread.daemon = True
        self.xmlrpc_server_0_thread.start()

        #blocks_tags_strobe blocks need to come before slow radio startup commands for some silly reason
        #self.blocks_tags_strobe_0_0 = blocks.tags_strobe(gr.sizeof_gr_complex*1, pmt.to_pmt({"num_bins": num_bins, "samp_rate": samp_rate, "num_integrations": num_integrations, "motor_az": motor_az, "motor_el": motor_el, "freq": freq, "tsys": [float(n) for n in tsys], "tcal": [float(n) for n in tcal], "cal_pwr": [float(n) for n in cal_pwr], "vlsr": vlsr, "glat": glat, "glon": glon, "soutrack": soutrack, "bsw": beam_switch, "cal_on":cal_on}), tag_period, pmt.intern("metadata"))
        self.blocks_tags_strobe_0_0 = blocks.tags_strobe(gr.sizeof_gr_complex*1, pmt.to_pmt({"num_bins": num_bins, "samp_rate": samp_rate, "num_integrations": num_integrations, "motor_az": motor_az, "motor_el": motor_el, "freq": freq, "tsys": [float(n) for n in tsys], "tcal": [float(n) for n in tcal], "cal_pwr": [float(n) for n in cal_pwr], "vlsr": vlsr, "glat": glat, "glon": glon, "soutrack": soutrack, "bsw": beam_switch}), tag_period, pmt.intern("metadata"))
        self.blocks_tags_strobe_0 = blocks.tags_strobe(gr.sizeof_gr_complex*1, pmt.to_pmt(float(freq)), tag_period, pmt.intern("rx_freq"))


        self.uhd_usrp_source_1 = uhd.usrp_source(
            ",".join(("addr=172.25.14.11", '')),
            uhd.stream_args(
                cpu_format="fc32",
                args='',
                channels=list(range(0,2)),
            ),
        )

        
        self.uhd_usrp_source_1.set_clock_source("external")
        self.uhd_usrp_source_1.set_time_source("external")
        self.uhd_usrp_source_1.set_samp_rate(samp_rate)
        _last_pps_time = self.uhd_usrp_source_1.get_time_last_pps().get_real_secs()
        # Poll get_time_last_pps() every 50 ms until a change is seen
        while(self.uhd_usrp_source_1.get_time_last_pps().get_real_secs() == _last_pps_time):
            time.sleep(0.05)
        # Set the time to PC time on next PPS
        self.uhd_usrp_source_1.set_time_next_pps(uhd.time_spec(int(time.time()) + 1.0))
        # Sleep 1 second to ensure next PPS has come
        time.sleep(1)


        ###### initial USRP channel Setup
        self.uhd_usrp_source_1.set_center_freq(rf_freq, 0)
        self.uhd_usrp_source_1.set_antenna("RX2", 0)
        self.uhd_usrp_source_1.set_bandwidth(samp_rate, 0)
        self.uhd_usrp_source_1.set_gain(rf_gain, 0)
        self.uhd_usrp_source_1.set_auto_dc_offset(True, 0)

        self.uhd_usrp_source_1.set_center_freq(rf_freq, 1)
        self.uhd_usrp_source_1.set_antenna("RX2", 1)
        self.uhd_usrp_source_1.set_bandwidth(samp_rate, 1)
        self.uhd_usrp_source_1.set_gain(rf_gain, 1)
        self.uhd_usrp_source_1.set_auto_dc_offset(True, 1)

        ##### Manually Configure USRP GPIO
        self.uhd_usrp_source_1.set_gpio_attr('FP0A', 'CTRL', 0x000, 0xFFF ^ calibrator_mask)  #set pins 2 and 3 manual
        self.uhd_usrp_source_1.set_gpio_attr('FP0A', 'DDR', 0xFFF, calibrator_mask) #set pins 2 and 3 as output
        self.uhd_usrp_source_1.set_gpio_attr('FP0A', 'OUT', 0x000 , calibrator_mask)

        ##### configure LO sharing

        #self.uhd_usrp_source_1.set_lo_source('internal', uhd.ALL_LOS, 0)
        #self.uhd_usrp_source_1.set_lo_export_enabled(True, uhd.ALL_LOS, 0)
        #self.uhd_usrp_source_1.set_lo_source('external', uhd.ALL_LOS, 1)
        #self.uhd_usrp_source_1.set_lo_export_enabled(False, uhd.ALL_LOS, 1)

        ##### timed tuning command 

        self.uhd_usrp_source_1.clear_command_time()
        now_time = self.uhd_usrp_source_1.get_time_last_pps()
        self.uhd_usrp_source_1.set_command_time(now_time + uhd.time_spec(1.0)) #occur at next second or ASAP
        
        #self.uhd_usrp_source_1.set_center_freq(self.rf_freq, 0)
        self.uhd_usrp_source_1.set_center_freq(uhd.tune_request(self.rf_freq,self.samp_rate*0.6), 0)
        #self.uhd_usrp_source_1.set_center_freq(self.rf_freq, 1)
        self.uhd_usrp_source_1.set_center_freq(uhd.tune_request(self.rf_freq,self.samp_rate*0.6), 1)

        self.uhd_usrp_source_1.clear_command_time()




        self.calibrator_control_strobe_0 = calibrator_control_strobe.blk(cal_mask=calibrator_mask, cal_state=cal_on, cal_interval=tag_period/samp_rate, samp_rate=samp_rate)
        #self.calibrator_control_strobe = calibrator_control_strobe.msg_blk(calibrator_mask=calibrator_mask, cal_state=cal_on)
        self.blocks_vector_to_streams_0 = blocks.vector_to_streams(gr.sizeof_gr_complex*num_bins, (num_channels**2))
        self.blocks_streams_to_vector_1 = blocks.streams_to_vector(gr.sizeof_gr_complex*num_bins, (num_channels**2))
        self.blocks_streams_to_vector_0_0_0 = blocks.streams_to_vector(gr.sizeof_gr_complex*num_bins, 4)
        self.blocks_streams_to_vector_0 = blocks.streams_to_vector(gr.sizeof_gr_complex*1, 2)
        self.blocks_multiply_const_xx_0_0_0_0 = blocks.multiply_const_cc(1.0/float(num_integrations), (num_bins*(num_channels**2)))
        self.blocks_multiply_const_vxx_1_0_0_0 = blocks.multiply_const_vcc(cal_values[2])
        self.blocks_multiply_const_vxx_1_0_0 = blocks.multiply_const_vcc(cal_values[1])
        self.blocks_multiply_const_vxx_1_0 = blocks.multiply_const_vcc(cal_values[3])
        self.blocks_multiply_const_vxx_1 = blocks.multiply_const_vcc(cal_values[0])
        #self.blocks_message_strobe_0 = blocks.message_strobe(pmt.to_pmt(is_running), int(tag_period/samp_rate*1000))
        self.blocks_integrate_xx_0 = blocks.integrate_cc(num_integrations, (num_bins*(num_channels**2)))
        self.blocks_add_xx_0_0_0 = blocks.add_vcc(1)
        self.blocks_add_xx_0_0 = blocks.add_vcc(1)
        #self.add_clock_tags_0 = add_clock_tags.clk(nsamps=tag_period)
        #self.add_clock_tags = add_clock_tags.clk(nsamps=tag_period)

        self.weighted_overlap_fft_0_0 = weighted_overlap_fft.weighted_overlap_fft(
            fft_window=fft_window,
            num_bins=num_bins,
            num_integrations=num_integrations,
        )
        self.weighted_overlap_fft_0 = weighted_overlap_fft.weighted_overlap_fft(
            fft_window=fft_window,
            num_bins=num_bins,
            num_integrations=num_integrations,
        )
        self.covariance_matrix_1 = covariance_matrix.covariance_matrix_block(num_channels=num_channels, vec_length=num_bins)

        ##################################################
        # Connections
        ##################################################
        #self.msg_connect((self.blocks_message_strobe_0, 'strobe'), (self.calibrator_control_strobe, 'strobe'))
        self.msg_connect((self.calibrator_control_strobe_0, 'command'), (self.uhd_usrp_source_1, 'command'))
        self.connect((self.covariance_matrix_1, 2), (self.blocks_streams_to_vector_1, 2))
        self.connect((self.covariance_matrix_1, 3), (self.blocks_streams_to_vector_1, 3))
        self.connect((self.covariance_matrix_1, 1), (self.blocks_streams_to_vector_1, 1))
        self.connect((self.covariance_matrix_1, 0), (self.blocks_streams_to_vector_1, 0))
        self.connect((self.weighted_overlap_fft_0, 0), (self.covariance_matrix_1, 0))
        self.connect((self.weighted_overlap_fft_0_0, 0), (self.covariance_matrix_1, 1))
        #self.connect((self.add_clock_tags, 0), (self.blocks_add_xx_0_0, 1))
        #self.connect((self.add_clock_tags_0, 0), (self.blocks_add_xx_0_0_0, 1))

        self.connect((self.blocks_add_xx_0_0, 0), (self.weighted_overlap_fft_0, 0))
        self.connect((self.blocks_add_xx_0_0, 0), (self.blocks_streams_to_vector_0, 0))
        self.connect((self.blocks_add_xx_0_0_0, 0), (self.weighted_overlap_fft_0_0, 0))
        self.connect((self.blocks_add_xx_0_0_0, 0), (self.blocks_streams_to_vector_0, 1))

        self.connect((self.blocks_integrate_xx_0, 0), (self.blocks_multiply_const_xx_0_0_0_0, 0))
        self.connect((self.blocks_multiply_const_vxx_1, 0), (self.blocks_streams_to_vector_0_0_0, 0))
        self.connect((self.blocks_multiply_const_vxx_1_0, 0), (self.blocks_streams_to_vector_0_0_0, 3))
        self.connect((self.blocks_multiply_const_vxx_1_0_0, 0), (self.blocks_streams_to_vector_0_0_0, 1))
        self.connect((self.blocks_multiply_const_vxx_1_0_0_0, 0), (self.blocks_streams_to_vector_0_0_0, 2))
        self.connect((self.blocks_multiply_const_xx_0_0_0_0, 0), (self.blocks_vector_to_streams_0, 0))
        self.connect((self.blocks_multiply_const_xx_0_0_0_0, 0), (self.zeromq_pub_sink_2, 0))
        self.connect((self.blocks_multiply_const_xx_0_0_0_0, 0), (self.zeromq_pub_sink_2_0, 0))
        self.connect((self.blocks_streams_to_vector_0, 0), (self.zeromq_pub_sink_0, 0))
        self.connect((self.blocks_streams_to_vector_0, 0), (self.zeromq_pub_sink_0_0, 0))
        self.connect((self.blocks_streams_to_vector_0_0_0, 0), (self.zeromq_pub_sink_1, 0))
        self.connect((self.blocks_streams_to_vector_0_0_0, 0), (self.zeromq_pub_sink_1_0, 0))
        self.connect((self.blocks_streams_to_vector_1, 0), (self.blocks_integrate_xx_0, 0))
        self.connect((self.blocks_tags_strobe_0, 0), (self.blocks_add_xx_0_0, 0))
        self.connect((self.blocks_tags_strobe_0, 0), (self.blocks_add_xx_0_0_0, 0))
        self.connect((self.blocks_tags_strobe_0_0, 0), (self.blocks_add_xx_0_0, 2))
        self.connect((self.blocks_tags_strobe_0_0, 0), (self.blocks_add_xx_0_0_0, 2))
        self.connect((self.blocks_vector_to_streams_0, 0), (self.blocks_multiply_const_vxx_1, 0))
        self.connect((self.blocks_vector_to_streams_0, 3), (self.blocks_multiply_const_vxx_1_0, 0))
        self.connect((self.blocks_vector_to_streams_0, 1), (self.blocks_multiply_const_vxx_1_0_0, 0))
        self.connect((self.blocks_vector_to_streams_0, 2), (self.blocks_multiply_const_vxx_1_0_0_0, 0))
        #self.connect((self.uhd_usrp_source_1, 0), (self.add_clock_tags, 0))
        #self.connect((self.uhd_usrp_source_1, 1), (self.add_clock_tags_0, 0))
        self.connect((self.calibrator_control_strobe_0, 0), (self.blocks_add_xx_0_0, 1))
        self.connect((self.calibrator_control_strobe_0, 1), (self.blocks_add_xx_0_0_0, 1))
        self.connect((self.uhd_usrp_source_1, 0), (self.calibrator_control_strobe_0, 0))
        self.connect((self.uhd_usrp_source_1, 1), (self.calibrator_control_strobe_0, 1))


    def get_num_bins(self):
        return self.num_bins

    def set_num_bins(self, num_bins):
        self.num_bins = num_bins
        self.set_cal_values_imag(np.zeros((self.num_channels**2,self.num_bins)))
        self.set_cal_values_real(np.ones((self.num_channels**2,self.num_bins)))
        #self.set_cal_values(np.ones((self.num_channels**2,self.num_bins))+1j*np.zeros((self.num_channels**2,self.num_bins)))
        self.set_custom_window(self.sinc_samples*np.hamming(4*self.num_bins))
        self.set_fft_window(window.blackmanharris(self.num_bins))
        self.set_sinc_sample_locations(np.arange(-np.pi*4/2.0, np.pi*4/2.0, np.pi/self.num_bins))
        self.set_tag_period(self.num_bins*self.num_integrations)
        self.covariance_matrix_1.vec_length = self.num_bins
        self.weighted_overlap_fft_0.set_num_bins(self.num_bins)
        self.weighted_overlap_fft_0_0.set_num_bins(self.num_bins)
        #self.blocks_tags_strobe_0_0.set_value(pmt.to_pmt({"num_bins": self.num_bins, "samp_rate": self.samp_rate, "num_integrations": self.num_integrations, "motor_az": self.motor_az, "motor_el": self.motor_el, "freq": self.freq, "tsys": [float(n) for n in self.tsys], "tcal": [float(n) for n in self.tcal], "cal_pwr": [float(n) for n in self.cal_pwr], "vlsr": self.vlsr, "glat": self.glat, "glon": self.glon, "soutrack": self.soutrack, "bsw": self.beam_switch, "cal_on":self.cal_on}))
        self.blocks_tags_strobe_0_0.set_value(pmt.to_pmt({"num_bins": self.num_bins, "samp_rate": self.samp_rate, "num_integrations": self.num_integrations, "motor_az": self.motor_az, "motor_el": self.motor_el, "freq": self.freq, "tsys": [float(n) for n in self.tsys], "tcal": [float(n) for n in self.tcal], "cal_pwr": [float(n) for n in self.cal_pwr], "vlsr": self.vlsr, "glat": self.glat, "glon": self.glon, "soutrack": self.soutrack, "bsw": self.beam_switch}))


    def get_num_integrations(self):
        return self.num_integrations

    def set_num_integrations(self, num_integrations):
        self.num_integrations = num_integrations
        self.set_tag_period(self.num_bins*self.num_integrations)
        #self.blocks_tags_strobe_0_0.set_value(pmt.to_pmt({"num_bins": self.num_bins, "samp_rate": self.samp_rate, "num_integrations": self.num_integrations, "motor_az": self.motor_az, "motor_el": self.motor_el, "freq": self.freq, "tsys": [float(n) for n in self.tsys], "tcal": [float(n) for n in self.tcal], "cal_pwr": [float(n) for n in self.cal_pwr], "vlsr": self.vlsr, "glat": self.glat, "glon": self.glon, "soutrack": self.soutrack, "bsw": self.beam_switch, "cal_on":self.cal_on}))
        self.blocks_tags_strobe_0_0.set_value(pmt.to_pmt({"num_bins": self.num_bins, "samp_rate": self.samp_rate, "num_integrations": self.num_integrations, "motor_az": self.motor_az, "motor_el": self.motor_el, "freq": self.freq, "tsys": [float(n) for n in self.tsys], "tcal": [float(n) for n in self.tcal], "cal_pwr": [float(n) for n in self.cal_pwr], "vlsr": self.vlsr, "glat": self.glat, "glon": self.glon, "soutrack": self.soutrack, "bsw": self.beam_switch}))
        self.weighted_overlap_fft_0.set_num_integrations(self.num_integrations)
        self.weighted_overlap_fft_0_0.set_num_integrations(self.num_integrations)
        self.blocks_multiply_const_xx_0_0_0_0.set_k(1.0/float(self.num_integrations))

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

    def get_freq(self):
        return self.freq

    def set_freq(self, freq):
        self.freq = freq
        self.set_rf_freq(self.freq)
        self.blocks_tags_strobe_0.set_value(pmt.to_pmt(float(self.freq)))
        #self.blocks_tags_strobe_0_0.set_value(pmt.to_pmt({"num_bins": self.num_bins, "samp_rate": self.samp_rate, "num_integrations": self.num_integrations, "motor_az": self.motor_az, "motor_el": self.motor_el, "freq": self.freq, "tsys": [float(n) for n in self.tsys], "tcal": [float(n) for n in self.tcal], "cal_pwr": [float(n) for n in self.cal_pwr], "vlsr": self.vlsr, "glat": self.glat, "glon": self.glon, "soutrack": self.soutrack, "bsw": self.beam_switch, "cal_on":self.cal_on}))
        self.blocks_tags_strobe_0_0.set_value(pmt.to_pmt({"num_bins": self.num_bins, "samp_rate": self.samp_rate, "num_integrations": self.num_integrations, "motor_az": self.motor_az, "motor_el": self.motor_el, "freq": self.freq, "tsys": [float(n) for n in self.tsys], "tcal": [float(n) for n in self.tcal], "cal_pwr": [float(n) for n in self.cal_pwr], "vlsr": self.vlsr, "glat": self.glat, "glon": self.glon, "soutrack": self.soutrack, "bsw": self.beam_switch}))

    def get_vlsr(self):
        return self.vlsr

    def set_vlsr(self, vlsr):
        self.vlsr = vlsr
        #self.blocks_tags_strobe_0_0.set_value(pmt.to_pmt({"num_bins": self.num_bins, "samp_rate": self.samp_rate, "num_integrations": self.num_integrations, "motor_az": self.motor_az, "motor_el": self.motor_el, "freq": self.freq, "tsys": [float(n) for n in self.tsys], "tcal": [float(n) for n in self.tcal], "cal_pwr": [float(n) for n in self.cal_pwr], "vlsr": self.vlsr, "glat": self.glat, "glon": self.glon, "soutrack": self.soutrack, "bsw": self.beam_switch, "cal_on":self.cal_on}))
        self.blocks_tags_strobe_0_0.set_value(pmt.to_pmt({"num_bins": self.num_bins, "samp_rate": self.samp_rate, "num_integrations": self.num_integrations, "motor_az": self.motor_az, "motor_el": self.motor_el, "freq": self.freq, "tsys": [float(n) for n in self.tsys], "tcal": [float(n) for n in self.tcal], "cal_pwr": [float(n) for n in self.cal_pwr], "vlsr": self.vlsr, "glat": self.glat, "glon": self.glon, "soutrack": self.soutrack, "bsw": self.beam_switch}))

    def get_tsys(self):
        return self.tsys

    def set_tsys(self, tsys):
        self.tsys = tsys
        #self.blocks_tags_strobe_0_0.set_value(pmt.to_pmt({"num_bins": self.num_bins, "samp_rate": self.samp_rate, "num_integrations": self.num_integrations, "motor_az": self.motor_az, "motor_el": self.motor_el, "freq": self.freq, "tsys": [float(n) for n in self.tsys], "tcal": [float(n) for n in self.tcal], "cal_pwr": [float(n) for n in self.cal_pwr], "vlsr": self.vlsr, "glat": self.glat, "glon": self.glon, "soutrack": self.soutrack, "bsw": self.beam_switch, "cal_on":self.cal_on}))
        self.blocks_tags_strobe_0_0.set_value(pmt.to_pmt({"num_bins": self.num_bins, "samp_rate": self.samp_rate, "num_integrations": self.num_integrations, "motor_az": self.motor_az, "motor_el": self.motor_el, "freq": self.freq, "tsys": [float(n) for n in self.tsys], "tcal": [float(n) for n in self.tcal], "cal_pwr": [float(n) for n in self.cal_pwr], "vlsr": self.vlsr, "glat": self.glat, "glon": self.glon, "soutrack": self.soutrack, "bsw": self.beam_switch}))

    def get_tcal(self):
        return self.tcal

    def set_tcal(self, tcal):
        self.tcal = tcal
        #self.blocks_tags_strobe_0_0.set_value(pmt.to_pmt({"num_bins": self.num_bins, "samp_rate": self.samp_rate, "num_integrations": self.num_integrations, "motor_az": self.motor_az, "motor_el": self.motor_el, "freq": self.freq, "tsys": [float(n) for n in self.tsys], "tcal": [float(n) for n in self.tcal], "cal_pwr": [float(n) for n in self.cal_pwr], "vlsr": self.vlsr, "glat": self.glat, "glon": self.glon, "soutrack": self.soutrack, "bsw": self.beam_switch, "cal_on":self.cal_on}))
        self.blocks_tags_strobe_0_0.set_value(pmt.to_pmt({"num_bins": self.num_bins, "samp_rate": self.samp_rate, "num_integrations": self.num_integrations, "motor_az": self.motor_az, "motor_el": self.motor_el, "freq": self.freq, "tsys": [float(n) for n in self.tsys], "tcal": [float(n) for n in self.tcal], "cal_pwr": [float(n) for n in self.cal_pwr], "vlsr": self.vlsr, "glat": self.glat, "glon": self.glon, "soutrack": self.soutrack, "bsw": self.beam_switch}))

    def get_tag_period(self):
        return self.tag_period

    def set_tag_period(self, tag_period):
        self.tag_period = tag_period
        #self.add_clock_tags.nsamps = self.tag_period
        #self.add_clock_tags_0.nsamps = self.tag_period
        #self.blocks_message_strobe_0.set_period((int(self.tag_period/self.samp_rate*1000)))
        self.blocks_tags_strobe_0.set_nsamps(self.tag_period)
        self.blocks_tags_strobe_0_0.set_nsamps(self.tag_period)
        self.calibrator_control_strobe_0.cal_interval = self.tag_period/self.samp_rate

    def get_soutrack(self):
        return self.soutrack

    def set_soutrack(self, soutrack):
        self.soutrack = soutrack
        #self.blocks_tags_strobe_0_0.set_value(pmt.to_pmt({"num_bins": self.num_bins, "samp_rate": self.samp_rate, "num_integrations": self.num_integrations, "motor_az": self.motor_az, "motor_el": self.motor_el, "freq": self.freq, "tsys": [float(n) for n in self.tsys], "tcal": [float(n) for n in self.tcal], "cal_pwr": [float(n) for n in self.cal_pwr], "vlsr": self.vlsr, "glat": self.glat, "glon": self.glon, "soutrack": self.soutrack, "bsw": self.beam_switch, "cal_on":self.cal_on}))
        self.blocks_tags_strobe_0_0.set_value(pmt.to_pmt({"num_bins": self.num_bins, "samp_rate": self.samp_rate, "num_integrations": self.num_integrations, "motor_az": self.motor_az, "motor_el": self.motor_el, "freq": self.freq, "tsys": [float(n) for n in self.tsys], "tcal": [float(n) for n in self.tcal], "cal_pwr": [float(n) for n in self.cal_pwr], "vlsr": self.vlsr, "glat": self.glat, "glon": self.glon, "soutrack": self.soutrack, "bsw": self.beam_switch}))

    def get_samp_rate(self):
        return self.samp_rate

    def set_samp_rate(self, samp_rate):
        #note that we are not yet implementing phase locking of the channels here. X300 with UBX cards can do that so need to return to this
        self.samp_rate = samp_rate
        #self.blocks_tags_strobe_0_0.set_value(pmt.to_pmt({"num_bins": self.num_bins, "samp_rate": self.samp_rate, "num_integrations": self.num_integrations, "motor_az": self.motor_az, "motor_el": self.motor_el, "freq": self.freq, "tsys": [float(n) for n in self.tsys], "tcal": [float(n) for n in self.tcal], "cal_pwr": [float(n) for n in self.cal_pwr], "vlsr": self.vlsr, "glat": self.glat, "glon": self.glon, "soutrack": self.soutrack, "bsw": self.beam_switch, "cal_on":self.cal_on}))
        self.blocks_tags_strobe_0_0.set_value(pmt.to_pmt({"num_bins": self.num_bins, "samp_rate": self.samp_rate, "num_integrations": self.num_integrations, "motor_az": self.motor_az, "motor_el": self.motor_el, "freq": self.freq, "tsys": [float(n) for n in self.tsys], "tcal": [float(n) for n in self.tcal], "cal_pwr": [float(n) for n in self.cal_pwr], "vlsr": self.vlsr, "glat": self.glat, "glon": self.glon, "soutrack": self.soutrack, "bsw": self.beam_switch}))
        self.uhd_usrp_source_1.set_samp_rate(self.samp_rate)
        #self.blocks_message_strobe_0.set_period((int(self.tag_period/self.samp_rate*1000)))
        self.calibrator_control_strobe_0.cal_interval = self.tag_period/self.samp_rate
        self.calibrator_control_strobe_0.samp_rate = self.samp_rate

        ##### timed tuning command 

        self.uhd_usrp_source_1.clear_command_time()
        now_time = self.uhd_usrp_source_1.get_time_last_pps()
        self.uhd_usrp_source_1.set_command_time(now_time + uhd.time_spec(1.0)) 

        self.uhd_usrp_source_1.set_bandwidth(self.samp_rate, 0)
        self.uhd_usrp_source_1.set_bandwidth(self.samp_rate, 1)
        self.uhd_usrp_source_1.set_center_freq(uhd.tune_request(self.rf_freq,self.samp_rate*0.6), 0)
        self.uhd_usrp_source_1.set_center_freq(uhd.tune_request(self.rf_freq,self.samp_rate*0.6), 1)

        self.uhd_usrp_source_1.clear_command_time()

    def get_rf_gain(self):
        return self.rf_gain

    def set_rf_gain(self, rf_gain):
        self.rf_gain = rf_gain
        self.uhd_usrp_source_1.set_gain(self.rf_gain, 0)
        self.uhd_usrp_source_1.set_gain(self.rf_gain, 1)

    def get_rf_freq(self):
        return self.rf_freq

    def set_rf_freq(self, rf_freq):
        self.rf_freq = rf_freq

        ##### timed tuning command 

        self.uhd_usrp_source_1.clear_command_time()
        now_time = self.uhd_usrp_source_1.get_time_last_pps()
        self.uhd_usrp_source_1.set_command_time(now_time + uhd.time_spec(1.0)) #occur at next second or ASAP

        self.rf_freq = rf_freq
        #self.uhd_usrp_source_1.set_center_freq(self.rf_freq, 0)
        self.uhd_usrp_source_1.set_center_freq(uhd.tune_request(self.rf_freq,self.samp_rate*0.6), 0)
        #self.uhd_usrp_source_1.set_center_freq(self.rf_freq, 1)
        self.uhd_usrp_source_1.set_center_freq(uhd.tune_request(self.rf_freq,self.samp_rate*0.6), 1)

        self.uhd_usrp_source_1.clear_command_time()

    def get_num_channels(self):
        return self.num_channels

    def set_num_channels(self, num_channels):
        self.num_channels = num_channels
        self.set_cal_pwr(np.array([1]*self.num_channels**2))
        self.set_cal_values_imag(np.zeros((self.num_channels**2,self.num_bins)))
        self.set_cal_values_real(np.ones((self.num_channels**2,self.num_bins)))
        #self.set_cal_values(np.ones((self.num_channels**2,self.num_bins))+1j*np.zeros((self.num_channels**2,self.num_bins)))
        self.set_tcal(np.array([290]*self.num_channels))
        self.set_tsys(np.array([171]*self.num_channels))
        self.covariance_matrix_1.num_channels = self.num_channels

    def get_motor_el(self):
        return self.motor_el

    def set_motor_el(self, motor_el):
        self.motor_el = motor_el
        #self.blocks_tags_strobe_0_0.set_value(pmt.to_pmt({"num_bins": self.num_bins, "samp_rate": self.samp_rate, "num_integrations": self.num_integrations, "motor_az": self.motor_az, "motor_el": self.motor_el, "freq": self.freq, "tsys": [float(n) for n in self.tsys], "tcal": [float(n) for n in self.tcal], "cal_pwr": [float(n) for n in self.cal_pwr], "vlsr": self.vlsr, "glat": self.glat, "glon": self.glon, "soutrack": self.soutrack, "bsw": self.beam_switch, "cal_on":self.cal_on}))
        self.blocks_tags_strobe_0_0.set_value(pmt.to_pmt({"num_bins": self.num_bins, "samp_rate": self.samp_rate, "num_integrations": self.num_integrations, "motor_az": self.motor_az, "motor_el": self.motor_el, "freq": self.freq, "tsys": [float(n) for n in self.tsys], "tcal": [float(n) for n in self.tcal], "cal_pwr": [float(n) for n in self.cal_pwr], "vlsr": self.vlsr, "glat": self.glat, "glon": self.glon, "soutrack": self.soutrack, "bsw": self.beam_switch}))

    def get_motor_az(self):
        return self.motor_az

    def set_motor_az(self, motor_az):
        self.motor_az = motor_az
        #self.blocks_tags_strobe_0_0.set_value(pmt.to_pmt({"num_bins": self.num_bins, "samp_rate": self.samp_rate, "num_integrations": self.num_integrations, "motor_az": self.motor_az, "motor_el": self.motor_el, "freq": self.freq, "tsys": [float(n) for n in self.tsys], "tcal": [float(n) for n in self.tcal], "cal_pwr": [float(n) for n in self.cal_pwr], "vlsr": self.vlsr, "glat": self.glat, "glon": self.glon, "soutrack": self.soutrack, "bsw": self.beam_switch, "cal_on":self.cal_on}))
        self.blocks_tags_strobe_0_0.set_value(pmt.to_pmt({"num_bins": self.num_bins, "samp_rate": self.samp_rate, "num_integrations": self.num_integrations, "motor_az": self.motor_az, "motor_el": self.motor_el, "freq": self.freq, "tsys": [float(n) for n in self.tsys], "tcal": [float(n) for n in self.tcal], "cal_pwr": [float(n) for n in self.cal_pwr], "vlsr": self.vlsr, "glat": self.glat, "glon": self.glon, "soutrack": self.soutrack, "bsw": self.beam_switch}))

    def get_is_running(self):
        return self.is_running

    def set_is_running(self, is_running):
        self.is_running = is_running
        self.blocks_message_strobe_0.set_msg(pmt.to_pmt(self.is_running))

    def get_glon(self):
        return self.glon

    def set_glon(self, glon):
        self.glon = glon
        #self.blocks_tags_strobe_0_0.set_value(pmt.to_pmt({"num_bins": self.num_bins, "samp_rate": self.samp_rate, "num_integrations": self.num_integrations, "motor_az": self.motor_az, "motor_el": self.motor_el, "freq": self.freq, "tsys": [float(n) for n in self.tsys], "tcal": [float(n) for n in self.tcal], "cal_pwr": [float(n) for n in self.cal_pwr], "vlsr": self.vlsr, "glat": self.glat, "glon": self.glon, "soutrack": self.soutrack, "bsw": self.beam_switch, "cal_on":self.cal_on}))
        self.blocks_tags_strobe_0_0.set_value(pmt.to_pmt({"num_bins": self.num_bins, "samp_rate": self.samp_rate, "num_integrations": self.num_integrations, "motor_az": self.motor_az, "motor_el": self.motor_el, "freq": self.freq, "tsys": [float(n) for n in self.tsys], "tcal": [float(n) for n in self.tcal], "cal_pwr": [float(n) for n in self.cal_pwr], "vlsr": self.vlsr, "glat": self.glat, "glon": self.glon, "soutrack": self.soutrack, "bsw": self.beam_switch}))

    def get_glat(self):
        return self.glat

    def set_glat(self, glat):
        self.glat = glat
        #self.blocks_tags_strobe_0_0.set_value(pmt.to_pmt({"num_bins": self.num_bins, "samp_rate": self.samp_rate, "num_integrations": self.num_integrations, "motor_az": self.motor_az, "motor_el": self.motor_el, "freq": self.freq, "tsys": [float(n) for n in self.tsys], "tcal": [float(n) for n in self.tcal], "cal_pwr": [float(n) for n in self.cal_pwr], "vlsr": self.vlsr, "glat": self.glat, "glon": self.glon, "soutrack": self.soutrack, "bsw": self.beam_switch, "cal_on":self.cal_on}))
        self.blocks_tags_strobe_0_0.set_value(pmt.to_pmt({"num_bins": self.num_bins, "samp_rate": self.samp_rate, "num_integrations": self.num_integrations, "motor_az": self.motor_az, "motor_el": self.motor_el, "freq": self.freq, "tsys": [float(n) for n in self.tsys], "tcal": [float(n) for n in self.tcal], "cal_pwr": [float(n) for n in self.cal_pwr], "vlsr": self.vlsr, "glat": self.glat, "glon": self.glon, "soutrack": self.soutrack, "bsw": self.beam_switch}))

    def get_fft_window(self):
        return self.fft_window

    def set_fft_window(self, fft_window):
        self.fft_window = fft_window
        self.weighted_overlap_fft_0.set_fft_window(self.fft_window)
        self.weighted_overlap_fft_0_0.set_fft_window(self.fft_window)

    def get_custom_window(self):
        return self.custom_window

    def set_custom_window(self, custom_window):
        self.custom_window = custom_window

    def get_calibrator_mask(self):
        return self.calibrator_mask

    def set_calibrator_mask(self, calibrator_mask):
        self.calibrator_mask = calibrator_mask
        #self.calibrator_control_strobe.calibrator_mask = self.calibrator_mask
        self.calibrator_control_strobe_0.cal_mask = self.calibrator_mask
        ##### Configure USRP GPIO (not on the fly though, that's silly)
        #self.uhd_usrp_source_1.set_gpio_attr('FP0A', 'CTRL', 0x000, 0xFFF ^ calibrator_mask)  #set pins 2 and 3 manual
        #self.uhd_usrp_source_1.set_gpio_attr('FP0A', 'DDR', 0xFFF, calibrator_mask) #set pins 2 and 3 as output
        #self.uhd_usrp_source_1.set_gpio_attr('FP0A', 'OUT', 0x000 , calibrator_mask)
        
    def get_cal_values_real(self):
        return self.cal_values_real

    def set_cal_values_real(self, cal_values_real):
        self.cal_values_real = np.array(cal_values_real)
        self.set_cal_values(self.cal_values_real+1j*self.cal_values_imag)

    def get_cal_values_imag(self):
        return self.cal_values_imag

    def set_cal_values_imag(self, cal_values_imag):
        self.cal_values_imag = np.array(cal_values_imag)
        self.set_cal_values(self.cal_values_real+1j*self.cal_values_imag)

    def get_cal_values(self):
        return self.cal_values

    def set_cal_values(self, cal_values):
        self.cal_values = np.array(cal_values)
        self.blocks_multiply_const_vxx_1.set_k(self.cal_values[0])
        self.blocks_multiply_const_vxx_1_0.set_k(self.cal_values[3])
        self.blocks_multiply_const_vxx_1_0_0.set_k(self.cal_values[1] )
        self.blocks_multiply_const_vxx_1_0_0_0.set_k(self.cal_values[2])

    def get_cal_pwr(self):
        return self.cal_pwr

    def set_cal_pwr(self, cal_pwr):
        self.cal_pwr = cal_pwr
        #self.blocks_tags_strobe_0_0.set_value(pmt.to_pmt({"num_bins": self.num_bins, "samp_rate": self.samp_rate, "num_integrations": self.num_integrations, "motor_az": self.motor_az, "motor_el": self.motor_el, "freq": self.freq, "tsys": [float(n) for n in self.tsys], "tcal": [float(n) for n in self.tcal], "cal_pwr": [float(n) for n in self.cal_pwr], "vlsr": self.vlsr, "glat": self.glat, "glon": self.glon, "soutrack": self.soutrack, "bsw": self.beam_switch, "cal_on":self.cal_on}))
        self.blocks_tags_strobe_0_0.set_value(pmt.to_pmt({"num_bins": self.num_bins, "samp_rate": self.samp_rate, "num_integrations": self.num_integrations, "motor_az": self.motor_az, "motor_el": self.motor_el, "freq": self.freq, "tsys": [float(n) for n in self.tsys], "tcal": [float(n) for n in self.tcal], "cal_pwr": [float(n) for n in self.cal_pwr], "vlsr": self.vlsr, "glat": self.glat, "glon": self.glon, "soutrack": self.soutrack, "bsw": self.beam_switch}))

    def get_cal_on(self):
        return self.cal_on

    def set_cal_on(self, cal_on):
        self.cal_on = cal_on
        #self.blocks_tags_strobe_0_0.set_value(pmt.to_pmt({"num_bins": self.num_bins, "samp_rate": self.samp_rate, "num_integrations": self.num_integrations, "motor_az": self.motor_az, "motor_el": self.motor_el, "freq": self.freq, "tsys": [float(n) for n in self.tsys], "tcal": [float(n) for n in self.tcal], "cal_pwr": [float(n) for n in self.cal_pwr], "vlsr": self.vlsr, "glat": self.glat, "glon": self.glon, "soutrack": self.soutrack, "bsw": self.beam_switch, "cal_on":self.cal_on}))
        self.blocks_tags_strobe_0_0.set_value(pmt.to_pmt({"num_bins": self.num_bins, "samp_rate": self.samp_rate, "num_integrations": self.num_integrations, "motor_az": self.motor_az, "motor_el": self.motor_el, "freq": self.freq, "tsys": [float(n) for n in self.tsys], "tcal": [float(n) for n in self.tcal], "cal_pwr": [float(n) for n in self.cal_pwr], "vlsr": self.vlsr, "glat": self.glat, "glon": self.glon, "soutrack": self.soutrack, "bsw": self.beam_switch}))
        #self.calibrator_control_strobe.cal_state = self.cal_on
        self.calibrator_control_strobe_0.cal_state = self.cal_on

    def get_beam_switch(self):
        return self.beam_switch

    def set_beam_switch(self, beam_switch):
        self.beam_switch = beam_switch
        #self.blocks_tags_strobe_0_0.set_value(pmt.to_pmt({"num_bins": self.num_bins, "samp_rate": self.samp_rate, "num_integrations": self.num_integrations, "motor_az": self.motor_az, "motor_el": self.motor_el, "freq": self.freq, "tsys": [float(n) for n in self.tsys], "tcal": [float(n) for n in self.tcal], "cal_pwr": [float(n) for n in self.cal_pwr], "vlsr": self.vlsr, "glat": self.glat, "glon": self.glon, "soutrack": self.soutrack, "bsw": self.beam_switch, "cal_on":self.cal_on}))
        self.blocks_tags_strobe_0_0.set_value(pmt.to_pmt({"num_bins": self.num_bins, "samp_rate": self.samp_rate, "num_integrations": self.num_integrations, "motor_az": self.motor_az, "motor_el": self.motor_el, "freq": self.freq, "tsys": [float(n) for n in self.tsys], "tcal": [float(n) for n in self.tcal], "cal_pwr": [float(n) for n in self.cal_pwr], "vlsr": self.vlsr, "glat": self.glat, "glon": self.glon, "soutrack": self.soutrack, "bsw": self.beam_switch}))



def argument_parser():
    parser = ArgumentParser()
    parser.add_argument(
        "--num-bins", dest="num_bins", type=intx, default=256,
        help="Set num_bins [default=%(default)r]")
    parser.add_argument(
        "--num-integrations", dest="num_integrations", type=intx, default=100000,
        help="Set num_integrations [default=%(default)r]")
    return parser


def main(top_block_cls=radio_process_dual_channel, options=None):
    if options is None:
        options = argument_parser().parse_args()
    if gr.enable_realtime_scheduling() != gr.RT_OK:
        gr.logger("realtime").warn("Error: failed to enable real-time scheduling.")
    tb = top_block_cls(num_bins=options.num_bins, num_integrations=options.num_integrations)

    def sig_handler(sig=None, frame=None):
        tb.stop()
        tb.wait()

        sys.exit(0)

    signal.signal(signal.SIGINT, sig_handler)
    signal.signal(signal.SIGTERM, sig_handler)

    tb.start()

    tb.wait()


if __name__ == '__main__':
    main()
