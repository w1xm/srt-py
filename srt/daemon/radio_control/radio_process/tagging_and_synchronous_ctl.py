"""

This block generates synchronous calibrator control commands timed to the edge of
integration periodsb to control the X300 radio GPIO. 
It also manages all metadata tagging to enable an applied lag for the calibrator state information

note that this block assumes a latency lower than the calibrator cycle time
and will break if that condition is not met. modified to add timestamps on 2024/12/23

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
    def __init__(self, num_channels=2, cal_mask=0xFFF, cal_state=0, integration_time=1.0, samp_rate=32e3, center_frequency=1.42e9, metadata_pmt=pmt.to_pmt({"num_bins":512})):
        gr.sync_block.__init__(
            self,
            name="metadata_tagging_and_control",
            in_sig=[np.complex64]*num_channels,
            out_sig=[np.complex64]*num_channels
        )

        #input parameters
        self.cal_mask = cal_mask
        self.integration_time = integration_time
        self.cal_state = cal_state
        self.samp_rate = samp_rate
        self.num_channels = num_channels
        self.metadata_pmt = metadata_pmt
        self.center_frequency=center_frequency

        #fixed derived variables

        self.calibrator_sample_interval = int(self.samp_rate * self.integration_time)

        self.last_cal_state = False
        self.rx_time = None
        self.next_cal_time = None

        self.offset = 0
        
        #self.message_port_register_in(pmt.intern('get_gpio'))
        self.message_port_register_out(pmt.intern('command'))
        self.message_port_register_out(pmt.intern('time_reference'))
        #self.set_msg_handler(pmt.intern('gpio_command'), self.handle_msg)

    def work(self, input_items, output_items):


        #when SDR first starts, capture its internal timestamp off the first sample. use to time all subsequent events
        #ONLY accept radio timestamp once. It gets resent upon tuning commands and thoroughly borks things

        if self.rx_time == None:

            tags = self.get_tags_in_window(0, 0, len(input_items[0]))

            for tag in tags:
                key = pmt.to_python(tag.key) # convert from PMT to python string
                if key == "rx_time":
                    self.rx_time = pmt.to_python(tag.value) 

                    rx_time_float = self.rx_time[0]+self.rx_time[1]
                    msg = pmt.cons(pmt.string_to_symbol('radio_start_time'),pmt.to_pmt(float(self.rx_time[0]+self.rx_time[1])))
                    self.message_port_pub(pmt.intern('time_reference'), msg) #issue message
                    #print('key entry:', key)
                    #print('rx_time:', self.rx_time[0],self.rx_time[1], type(self.rx_time))
                    #print('')

        ################################################
        # Main Work Function
        ################################################

        else:

            #determine what the sample number of the last sample in the input is
            nitems = len(input_items[0]) + self.nitems_written(0)
            #n_last_sample = (self.nitems_written(0) + len(input_items[0])) % self.calibrator_sample_interval

            #while there are integration period boundaries present
            while (nitems - self.offset) > self.calibrator_sample_interval:
                
                self.offset += self.calibrator_sample_interval


                current_rx_time = float(self.rx_time[0]+self.rx_time[1]) + float(self.offset)/self.samp_rate

                if self.next_cal_time:
                    if current_rx_time == self.next_cal_time: #actually only want this flag a full period after the calibrator switches
                        self.last_cal_state = self.cal_state
                        self.next_cal_time = None

                #generate tags to be applied to data (pmt.cons does not work for metadata here, needs to be dict)
                #we take in all the radio state from an external metadata constructor EXCEPT for cal state 
                #since we really want that to line up with the transition.


                key = pmt.intern('metadata')
                value = self.metadata_pmt
                value = pmt.dict_add(value, pmt.to_pmt('cal_on'),pmt.to_pmt(int(self.last_cal_state)))

                #apply tags

                for i in range(self.num_channels):
                    self.add_item_tag(i, self.offset,key,value)
                    #self.add_item_tag(i, self.offset, pmt.intern("rx_time"), make_time_pair(time.time()))
                    self.add_item_tag(i, self.offset, pmt.intern("rx_time"), make_time_pair(current_rx_time)) #true to radio  timestamp
                    self.add_item_tag(i, self.offset, pmt.intern("rx_freq"), pmt.to_pmt(float(self.center_frequency)))


                if self.last_cal_state != self.cal_state:

                    ########################################
                    #issue command to usrp for next state of calibrator, 
                    #needs to be a timed command so it ends up synced with the integration periods
                    #only do this if we are changing things
                    #######################################

                    #set command time for approx 1 cycle hence (winds up being less when recieved at SDR)
                    #I probably need to fix this to actually match the time as recorded by the SDR

                    rftime = time.time() - float(self.rx_time[0]+self.rx_time[1])  #get the actual exact time since the radio started sampling
                    current_num_integration_cycles = int((rftime+0.01)/self.integration_time) #number of cycles that have been completed before now with just a little padding to ensure theres time to send it
                    self.next_cal_time = float(self.rx_time[0]+self.rx_time[1]) + (current_num_integration_cycles+1)*self.integration_time #one cycle ahead of now

                    command_time = pmt.cons(pmt.from_uint64(int(self.next_cal_time)),pmt.from_double(self.next_cal_time-int(self.next_cal_time)))
                    #command_time = make_time_pair(self.next_cal_time)
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

                    #self.last_cal_state = self.cal_state

        for i in range(self.num_channels):
            output_items[i][:] = input_items[i]
        return len(output_items[0])
        