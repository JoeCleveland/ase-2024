use std::iter::Filter;

use crate::ring_buffer::RingBuffer;

pub struct CombFilter {
    pub buffer: RingBuffer<f32>,
    filter_type: FilterType,
    gain: f32,
    sample_rate_hz: f32
}

#[derive(Debug, Clone, Copy)]
pub enum FilterType {
    FIR,
    IIR,
}

#[derive(Debug, Clone, Copy)]
pub enum FilterParam {
    Gain,
    Delay,
}

#[derive(Debug, Clone)]
pub enum Error {
    InvalidValue { param: FilterParam, value: f32 }
}

impl CombFilter {
    pub fn new(filter_type: FilterType, max_delay_secs: f32, sample_rate_hz: f32, num_channels: usize) -> Self {
        CombFilter {
            filter_type: filter_type,
            buffer: RingBuffer::new((max_delay_secs * sample_rate_hz) as usize),
            gain: 0.0,
            sample_rate_hz: sample_rate_hz
        }
    }

    pub fn reset(&mut self) {
        self.buffer.reset()
    }

    pub fn process(&mut self, input: &[&[f32]], output: &mut [&mut [f32]]) {
        for channel in 0..input.len() {
            for sample in 0..input[channel].len() {
                if matches!(self.filter_type, FilterType::FIR) {
                    self.buffer.push(input[channel][sample]);
                    output[channel][sample] = self.buffer.pop() * self.gain + input[channel][sample];
                }
                if matches!(self.filter_type, FilterType::IIR) {
                    output[channel][sample] = self.buffer.pop() * self.gain + input[channel][sample];
                    self.buffer.push(output[channel][sample]);
                } 
            }
        }
    }

    pub fn set_param(&mut self, param: FilterParam, value: f32) -> Result<(), Error> {
        if matches!(param, FilterParam::Gain) {
            self.gain = value;
            Result::Ok(())
        } else if matches!(param, FilterParam::Delay) {
            self.buffer.set_write_index((value * self.sample_rate_hz) as usize);
            Result::Ok(())
        } else {
            Result::Err(Error::InvalidValue { param: param, value: value })
        }
    }

    pub fn get_param(&self, param: FilterParam) -> f32 {
        if matches!(param, FilterParam::Gain) {
            self.gain
        } else {
            (self.buffer.get_write_index() - self.buffer.get_read_index()) as f32
        }

    }

    // TODO: feel free to define other functions for your own use
}

// TODO: feel free to define other types (here or in other modules) for your own use
