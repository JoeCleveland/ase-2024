use crate::ring_buffer::RingBuffer;
use crate::lfo::LFO;

pub struct Vibrato {
    buffer: RingBuffer<f32>,
    lfo: LFO,
    delay: f32,
    sample_rate: usize
}

impl Vibrato {
    pub fn new(freq: f32, max_delay_secs: f32, sample_rate: usize) -> Self {
        Vibrato {
            buffer: RingBuffer::new((max_delay_secs * sample_rate as f32) as usize),
            lfo: LFO::new(freq, sample_rate),
            delay: 0.0,
            sample_rate: sample_rate,
        }
    }

    pub fn set_delay(&mut self, delay_in_secs: f32) {
        if delay_in_secs < 0.0 {
            panic!("Delay must be positive");
        }
        self.delay = delay_in_secs * self.sample_rate as f32;
    }

    pub fn set_freq(&mut self, freq: f32) {
        if freq < 0.0 {
            panic!("Vibrato freq must be positive");
        }
        self.lfo = LFO::new(freq, self.sample_rate);
    }

    pub fn process_block(&mut self, input: &[f32], output: &mut [f32]) {
        let mut lfo_samples = vec![0.0; input.len()];
        self.lfo.get_block(lfo_samples.as_mut_slice());

        for i in 0..input.len() {
            self.buffer.set_read_index(self.buffer.get_write_index() as i32 - self.delay as i32);
            self.buffer.push(input[i]);
            let offset = (lfo_samples[i]) * self.delay;
            output[i] = self.buffer.get_frac(offset);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const EPSILON: f32 = 0.00001;

    fn sine_gen(freq: f32, Fs: f32, length: usize) -> Vec<f32> {
        let mut output = vec![0f32; length];
        for i in 0..length {
            output[i] = f32::sin((i as f32) * 2.0 * freq * std::f32::consts::PI / Fs);
        }
        return output;
    }

    #[test]
    fn output_equals_delayed_input() {

        //NOTE: The way I've implemented my vibrato, the amplitude of the effect is determined by the delay parameter
        // Therefore, to test the no-amplitude case, the delay is set to 0 and the input and output samples at the same index are identical
        let length = 16000;
        let delay = 0.0;
        let sample_rate = 16000;

        let mut vibrato = Vibrato::new(1.0, 1.0, sample_rate);
        vibrato.set_delay(delay);
        let input = sine_gen(1.0, 16000.0, length);
        let mut output = vec![0.0; length];

        vibrato.process_block(&input, &mut output);

        let delay_in_samples = (length as f32 * delay) as usize; //This will be zero
        let check_range = length - delay_in_samples;
        for i in 0..check_range {
            assert!(f32::abs(input[i] - output[i + delay_in_samples]) <= EPSILON);
        }
    }

    #[test]
    fn dc_in_eq_dc_out() {
        let length = 16000;
        let delay = 0.314159;
        let sample_rate = 16000;

        let mut vibrato = Vibrato::new(1.0, 1.0, sample_rate);
        vibrato.set_delay(delay);
        let input = vec![1.0; length];
        let mut output = vec![0.0; length];

        vibrato.process_block(&input, &mut output);

        // The direct DC output doesn't 'warm up' until we've passed the number of delay samples
        let start_sample = (sample_rate as f32 * delay) as usize;
        for i in start_sample..length {
            assert!(output[i] == 1.0);
        }
    }

    #[test]
    fn test_vary_block() {
        //Repeating the above DC to DC test with a change in block size
        let delay = 0.314159;
        let sample_rate = 16000;

        //We will run two blocks through, here are the two sizes
        let length_one = 16000;
        let length_two = 8000;

        let mut vibrato = Vibrato::new(1.0, 1.0, sample_rate);
        vibrato.set_delay(delay);

        let input_one = vec![1.0; length_one];
        let mut output_one = vec![0.0; length_one];
        let input_two = vec![1.0; length_two];
        let mut output_two = vec![0.0; length_two];

        vibrato.process_block(&input_one, &mut output_one);
        vibrato.process_block(&input_two, &mut output_two);

        let start_sample = (sample_rate as f32 * delay) as usize;
        for i in start_sample..(length_one + length_two) {
            if i >= length_one {
                assert!(output_two[i - length_one] == 1.0);
            } else {
                assert!(output_one[i] == 1.0);
            }
        }
    }

    #[test]
    fn test_zero_input() {
        let length = 16000;
        let delay = 0.314159;
        let sample_rate = 16000;

        let mut vibrato = Vibrato::new(1.0, 1.0, sample_rate);
        vibrato.set_delay(delay);
        let input = vec![0.0; length];
        let mut output = vec![0.0; length];

        vibrato.process_block(&input, &mut output);

        for i in 0..length {
            assert!(output[i] == 0.0);
        }
    }
}