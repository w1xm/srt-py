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

class blk(gr.sync_block):
    def __init__(self, cal_mask=0xFFF, cal_state=0, cal_interval=1.0, samp_rate=32e3):
        gr.sync_block.__init__(
            self,
            name="calibrator_control_and_timestamp",
            in_sig=[np.complex64, np.complex64],
            out_sig=[np.complex64, np.complex64]
        )

        #input parameters
        self.cal_mask = cal_mask
        self.cal_interval = cal_interval
        self.cal_state = cal_state
        self.samp_rate = samp_rate

        #fixed derived variables

        self.calibrator_sample_interval = int(self.samp_rate * self.cal_interval)

        self.last_cal_state = False
        self.rx_time = pmt.to_python(pmt.cons(pmt.from_uint64(int(0)),pmt.from_double(0)))

        #self.message_port_register_in(pmt.intern('get_gpio'))
        self.message_port_register_out(pmt.intern('command'))
        #self.set_msg_handler(pmt.intern('gpio_command'), self.handle_msg)

    def work(self, input_items, output_items):


        #when SDR first starts capture its timestamp off the first sample. use to time all subsequent events

        tags = self.get_tags_in_window(0, 0, len(input_items[0]))

        for tag in tags:
            key = pmt.to_python(tag.key) # convert from PMT to python string
            if key == "rx_time":
                self.rx_time = pmt.to_python(tag.value) # Note that the type(value) can be several things, it depends what PMT type it was
                #print('key entry:', key)
                #print('value:', self.rx_time[0],self.rx_time[1], type(self.rx_time))
                #print('')

        #determine what the sample number of the last sample in the input is

        n_last_sample = (self.nitems_written(0) + len(input_items[0])) % self.calibrator_sample_interval

        #determine calibrator state at samples being recieved

        #check if we are seeing a sample we are interested in adding a tag to
        if (n_last_sample-len(input_items[0])) <= 0:

            writeindex = len(input_items[0]) - n_last_sample

            #generate tag to be applied to data (pmt.cons does not work here, needs to be dict)

            key = pmt.intern('metadata')
            value = pmt.make_dict()
            value = pmt.dict_add(value, pmt.to_pmt('cal_on'),pmt.from_bool(self.last_cal_state))
            #value = pmt.cons(pmt.to_pmt('cal_on'), pmt.from_bool(self.last_cal_state))

            #apply tag

            self.add_item_tag(0, self.nitems_written(0) + writeindex,key,value)
            self.add_item_tag(1, self.nitems_written(0) + writeindex,key,value)

            #apply time

            self.add_item_tag(0, self.nitems_written(0) + writeindex, pmt.intern("rx_time"), make_time_pair(time.time()))
            self.add_item_tag(1, self.nitems_written(0) + writeindex, pmt.intern("rx_time"), make_time_pair(time.time()))

            #send message to define next calibrator state change

            ########################################
            #issue command to usrp for next flip of calibrator, 
            #needs to be a timed command
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


        output_items[0][:] = input_items[0]
        output_items[1][:] = input_items[1]
        return len(output_items[0])
        