"""
block to generate calibrator control commands for the X300 radio. 
note that this block assumes a latency lower than the calibrator cycle time
and will break if that condition is not met. modified to add timestamps on 2024/12/23

Embedded Python Blocks:

Each time this file is saved, GRC will instantiate the first class it finds
to get ports and parameters of your block. The arguments to __init__  will
be the parameters. All of them are required to have default values!
"""
import numpy as np
from gnuradio import gr 
import pmt
import time
import uhd

def make_time_pair(t):
    return pmt.make_tuple(
        pmt.to_pmt(int(np.trunc(t))), pmt.to_pmt(t - int(np.trunc(t)))
    )

class tagging_and_ctl(gr.sync_block):
    def __init__(self, num_channels=2, cal_mask=0xFFF, cal_state=0, cal_interval=1.0, samp_rate=32e3, center_frequency=1.42e9, metadata_pmt=pmt.to_pmt({"num_bins":512})):
        gr.sync_block.__init__(
            self,
            name="metadata_tagging_and_control",
            in_sig=[np.complex64]*num_channels,
            out_sig=[np.complex64]*num_channels
        )

        #input parameters
        self.cal_mask = cal_mask
        self.cal_interval = cal_interval
        self.cal_state = cal_state
        self.samp_rate = samp_rate
        self.num_channels = num_channels
        self.metadata_pmt = metadata_pmt
        self.center_frequency=center_frequency

        #fixed derived variables

        self.calibrator_sample_interval = int(self.samp_rate * self.cal_interval)

        self.last_cal_state = False
        self.rx_time = pmt.to_python(pmt.cons(pmt.from_uint64(int(0)),pmt.from_double(0)))

        #self.message_port_register_in(pmt.intern('get_gpio'))
        self.message_port_register_out(pmt.intern('command'))
        self.message_port_register_out(pmt.intern('time_reference'))
        #self.set_msg_handler(pmt.intern('gpio_command'), self.handle_msg)

    def work(self, input_items, output_items):


        #when SDR first starts capture its timestamp off the first sample. use to time all subsequent events

        tags = self.get_tags_in_window(0, 0, len(input_items[0]))

        for tag in tags:
            key = pmt.to_python(tag.key) # convert from PMT to python string
            if key == "rx_time":
                self.rx_time = pmt.to_python(tag.value) # Note that the type(value) can be several things, it depends what PMT type it was

                msg = pmt.make_dict()
                msg = pmt.dict_add(msg, pmt.to_pmt('radio_start_time'), pmt.to_python(self.rx_time))
                self.message_port_pub(pmt.intern('time_reference'), msg) #issue message
                #print('key entry:', key)
                #print('value:', self.rx_time[0],self.rx_time[1], type(self.rx_time))
                #print('')

        #determine what the sample number of the last sample in the input is

        n_last_sample = (self.nitems_written(0) + len(input_items[0])) % self.calibrator_sample_interval

        #determine calibrator state at samples being recieved

        #check if we are seeing a sample we are interested in adding a tag to
        if (n_last_sample-len(input_items[0])) <= 0:

            writeindex = len(input_items[0]) - n_last_sample

            #generate tags to be applied to data (pmt.cons does not work for metadata here, needs to be dict)
            #we take in all the radio state from an external metadata constructor EXCEPT for cal state 
            #since we really want that to line up with the transition.

            key = pmt.intern('metadata')
            value = self.metadata_pmt
            value = pmt.dict_add(value, pmt.to_pmt('cal_on'),pmt.to_pmt(int(self.last_cal_state)))
            #value = pmt.cons(pmt.to_pmt('cal_on'), pmt.from_bool(self.last_cal_state))

            #apply tags

            for i in range(self.num_channels):
                self.add_item_tag(i, self.nitems_written(0) + writeindex,key,value)
                self.add_item_tag(i, self.nitems_written(0) + writeindex, pmt.intern("rx_time"), make_time_pair(time.time()))
                self.add_item_tag(i, self.nitems_written(0) + writeindex, pmt.intern("rx_freq"), pmt.to_pmt(float(self.center_frequency)))


            if self.last_cal_state != self.cal_state:

                ########################################
                #issue command to usrp for next state of calibrator, 
                #needs to be a timed command so it ends up synced with the integration periods
                #only do this if we are changing things
                #######################################

                #set command time for approx 1 cycle hence (winds up being less when recieved at SDR)
                #I probably need to fix this to actually match the time as recorded by the SDR

                command_time = pmt.cons(pmt.from_uint64(int((self.nitems_written(0)+len(input_items[0]))/self.calibrator_sample_interval+self.cal_interval+self.rx_time[0])),pmt.from_double(self.rx_time[1]))
                msg = pmt.make_dict()
                msg = pmt.dict_add(msg, pmt.to_pmt('time'), command_time)

                self.message_port_pub(pmt.intern('command'), msg) #issue message

                #issue command to toggle gpio

                set_gpio = pmt.make_dict()
                set_gpio = pmt.dict_add(set_gpio, pmt.to_pmt('bank'), pmt.to_pmt('FP0A'))
                set_gpio = pmt.dict_add(set_gpio, pmt.to_pmt('attr'), pmt.to_pmt('OUT'))
                set_gpio = pmt.dict_add(set_gpio, pmt.to_pmt('value'), pmt.from_double(self.cal_state))
                set_gpio = pmt.dict_add(set_gpio, pmt.to_pmt('mask'), pmt.from_double(self.cal_mask))

                msg = pmt.make_dict()
                msg = pmt.dict_add(msg, pmt.to_pmt('gpio'), set_gpio)

                self.message_port_pub(pmt.intern('command'), msg) #issue message


                #clear command time 

                msg = pmt.make_dict()
                msg = pmt.dict_add(msg, pmt.to_pmt('time'), pmt.PMT_NIL)

                self.message_port_pub(pmt.intern('command'), msg) #issue message

                #self.message_port_pub(pmt.intern('command'), pmt.cons(pmt.to_pmt('time'), pmt.PMT_NIL))

                self.last_cal_state = self.cal_state

        for i in range(self.num_channels):
            output_items[i][:] = input_items[i]
        return len(output_items[0])
        