use rustfft::{num_complex::Complex32, FftPlanner};

use crate::ring_buffer::RingBuffer;

struct FastConvolver {
    mode: ConvolutionMode,
    //Used for time domain convolution:
    impulse_response: Vec<f32>,
    previous_block: Vec<f32>,
    //Used for freq domain convolution:
    freq_domain_delay_line: Vec<Complex32>,
    time_domain_input_buffer: Vec<f32>,
    filter_bank: Vec<Vec<Complex32>>,
    block_size: usize,
}

#[derive(Debug, Clone, Copy)]
pub enum ConvolutionMode {
    TimeDomain,
    FrequencyDomain { block_size: usize },
}

fn zero_pad<T: Default + Clone + Copy>(input: &[T], num_zeros: usize) -> Vec<T> {
    let mut v = vec![T::default(); input.len() + num_zeros];
    for i in 0..input.len() {
        v[i] = input[i];
    }
    return v;
}

fn build_filter_bank(impulse_response: Vec<f32>, block_size: usize) -> Vec<Vec<Complex32>> {
    let mut planner = FftPlanner::new();
    let fft = planner.plan_fft_forward(block_size * 2);
    let impulse_response: Vec<Complex32> = impulse_response.into_iter().map(|e| Complex32::from(e)).collect();

    let mut output = Vec::new();
    for i in 0..impulse_response.len()/block_size {
        let mut filter: Vec<Complex32> = vec![Complex32::new(0.0, 0.0); block_size * 2];
        filter[..block_size].clone_from_slice(&impulse_response[i*block_size..(i+1)*block_size]);
        fft.process(&mut filter);
        output.push(filter[..block_size+1].to_vec());
    }
    return output;
}

impl FastConvolver {
    pub fn new(impulse_response: &[f32], mode: ConvolutionMode) -> Self {
        let mut block_size = 1;
        let mut n_filters = 1;

        let vec_impulse: Vec<f32> = match mode {
            ConvolutionMode::FrequencyDomain { block_size: b_s } => {
                block_size = b_s;
                //Add some zeros to make impulse response a nice multiple of block size
                let nzeros = block_size - impulse_response.len() % block_size;
                n_filters = (impulse_response.len() + nzeros) / block_size;
                zero_pad(impulse_response, nzeros)
            },
            _ => impulse_response.to_vec(),
        };

        let impulse_len = vec_impulse.len();
        FastConvolver {
            mode,

            impulse_response: vec_impulse.clone(),
            previous_block: vec![0f32; impulse_len],

            freq_domain_delay_line: vec![Complex32::new(0.0, 0.0); (block_size + 1) * n_filters],
            time_domain_input_buffer: vec![0.0; block_size * 2],
            filter_bank: build_filter_bank(vec_impulse, block_size),
            block_size
        }
    }

    pub fn reset(&mut self) {
        todo!("implement")
    }

    fn process_linear(&mut self, input: &[f32], output: &mut [f32]) {
        let iL: i32 = self.impulse_response.len() as i32;
        let pbL: i32 = self.previous_block.len() as i32;
        for i in 0..input.len() as i32 {
            let mut accum = 0.0;
            for j in 0..self.impulse_response.len() as i32 {
                if i + j < iL - 1 {
                    accum += self.impulse_response[j as usize] * self.previous_block[(i + j + 1) as usize];
                } else {
                    accum += self.impulse_response[j as usize] * input[(i + j - iL + 1) as usize];
                }
            }
            output[i as usize] = accum;
        }
        //Shift previous block back by input size
        //Then put input at the end
        for i in 0..pbL - input.len() as i32 {
            self.previous_block[i as usize] = self.previous_block[(i + input.len() as i32) as usize]
        }
        for i in 0..input.len() as i32 {
            self.previous_block[(pbL - input.len() as i32 + i) as usize] = input[i as usize];
        }
    }

    fn process_freq_domain(&mut self, input: &[f32], output: &mut [f32]) {
        //shift input buffer and store input
        for i in 0..self.block_size {
            self.time_domain_input_buffer[i] = self.time_domain_input_buffer[i+self.block_size];
            self.time_domain_input_buffer[i+self.block_size] = input[i];
        }
        //shift the freq-domain delay-line
        for i in (self.block_size+1..self.freq_domain_delay_line.len()).rev() {
            self.freq_domain_delay_line[i] = self.freq_domain_delay_line[i - (self.block_size+1)];
        }
        //Compute FFT of the input buffer
        let mut planner = FftPlanner::new();
        let fft = planner.plan_fft_forward(self.block_size * 2);
        let mut freq_domain_input_buffer: Vec<Complex32>  = self.time_domain_input_buffer.clone().into_iter().map(|e| Complex32::from(e)).collect();
        fft.process(freq_domain_input_buffer.as_mut_slice());
        //Discard right half of FFT
        freq_domain_input_buffer = freq_domain_input_buffer[..self.block_size+1].to_vec();

        //Store result in the delay line
        self.freq_domain_delay_line[..self.block_size+1].copy_from_slice(&freq_domain_input_buffer);

        //Multiply each filter in the filter bank by each block of the delay-line
        //And accumulate the output
        let mut freq_domain_output = vec![Complex32::new(0.0, 0.0); self.block_size+1];
        for i in 0..self.filter_bank.len() {
            let delay_chunk = &self.freq_domain_delay_line[i*(self.block_size+1)..(i+1)*(self.block_size+1)];
            let result: Vec<Complex32> = self.filter_bank[i].iter().zip(delay_chunk).map(|(&F, &D)| F * D).collect();
            freq_domain_output = freq_domain_output.iter().zip(result).map(|(&O, R)| O + R).collect();
        }

        //Compute the inverse FFT
        //Input is B+1 but we need B*2 so mirror that guy on the mid axis
        let ifft = planner.plan_fft_inverse(self.block_size * 2);
        let mut ifft_output = vec![Complex32::new(0.0, 0.0); self.block_size * 2];
        for i in 0..self.block_size {
            ifft_output[i] = freq_domain_output[i];
            ifft_output[(self.block_size * 2) - i] = freq_domain_output[i];
        }
        ifft.process(&mut ifft_output);

        //Output the right half of the output buffer 
        let mut real_output: Vec<f32> = ifft_output.into_iter().map(|e| e.re).collect();
        output.clone_from_slice(&real_output[self.block_size..]);
 
    }

    pub fn process(&mut self, input: &[f32], output: &mut [f32]) {
        if matches!(self.mode, ConvolutionMode::TimeDomain) {
            self.process_linear(input, output);
        }
        if matches!(self.mode, ConvolutionMode::FrequencyDomain { .. }) {
            self.process_freq_domain(input, output);
        }
    }

    pub fn flush(&mut self, output: &mut [f32]) {
        //Run linear process on an input of zeros 
        let input = vec![0.0; self.previous_block.len()];
        self.process_linear(&input, output);
    }

    pub fn get_flush_length(&mut self) -> usize {
        return self.previous_block.len();
    }

}

#[cfg(test)]
mod tests {
    use std::cmp::{max, min};

    use super::*;
    use rand::Rng;

    /** 
     * identity: generate a random IR of length 51 samples and an impulse as input signal at sample index 3; make the input signal 10 samples long. 
     * Check the correct values of the 10 samples of the output signal.
    * flush: use the same signals as in identity, but now get the full tail of the impulse response 
    and check for correctness
    * blocksize: for an input signal of length 10000 samples, implement a test similar to identity with a succession of 
    different input/output block sizes (1, 13, 1023, 2048,1,17, 5000, 1897)
     * **/
     
    #[test]
    fn test_identity() {
        let mut rng = rand::thread_rng();
        let input_signal: Vec<f32> = vec![0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0];
        let mut output_signal: Vec<f32> = vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0];
        let mut impulse_response: Vec<f32> = vec![0f32; 51];
        let mut fast_convolver = FastConvolver::new(&impulse_response, ConvolutionMode::TimeDomain);
        for i in impulse_response.iter_mut() {
            *i = rng.gen_range(-1.0..=1.0);
        }
        fast_convolver.process(&input_signal, &mut output_signal);
        for i in 0..input_signal.len() {
            assert_eq!(input_signal[i], output_signal[i]);
        }
    }

    #[test]
    // fn test_flush() {
    //     let mut rng = rand::thread_rng();
    //     let input_signal: Vec<f32> = vec![0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0];
    //     let mut output_signal: Vec<f32> = vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0];
    //     let mut impulse_response: Vec<f32> = vec![0f32; 51];
    //     let mut fast_convolver = FastConvolver::new(&impulse_response, ConvolutionMode::TimeDomain);
    //     for i in impulse_response.iter_mut() {
    //         *i = rng.gen_range(-1.0..=1.0);
    //     }
    //     fast_convolver.process(&input_signal, &mut output_signal);
    //     for i in 0..input_signal.len() {
    //         assert_eq!(input_signal[i], output_signal[i]);
    //     }
    // }

    #[test]
    fn test_blocksize () {
        let mut rng = rand::thread_rng();
        let mut input_signal: Vec<f32> = vec![0f32; 10000];
        input_signal[3] = 1.0;
        let mut output_signal: Vec<f32> = vec![0f32; 10000];
        let mut impulse_response: Vec<f32> = vec![0f32; 51];
        let block_sizes = [1, 13, 1023, 2048,1,17, 5000, 1897];
        let mut fast_convolver = FastConvolver::new(&impulse_response, ConvolutionMode::TimeDomain);

        for block_size in block_sizes {
            for block_index in 0..=input_signal.len() / block_size {
                let block_start_index = block_index * block_size;
                let block_end_index = min(input_signal.len() - 1, (block_index + 1) * block_size);
                // println!("block_start_index: {}, block_end_index: {}", block_start_index, block_end_index);
                fast_convolver.process(&input_signal[block_start_index ..block_end_index], &mut output_signal[block_start_index..block_end_index]);
            }
            for i in 0..input_signal.len() {
                assert_eq!(input_signal[i], output_signal[i]);
            }
        }
    }

    fn test_identity_frequency_domain() {
        let mut rng = rand::thread_rng();
        let input_signal: Vec<f32> = vec![0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0];
        let mut output_signal: Vec<f32> = vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0];
        let mut impulse_response: Vec<f32> = vec![0f32; 51];
        let mut fast_convolver = FastConvolver::new(&impulse_response, ConvolutionMode::FrequencyDomain { block_size: (1024) }); // What should block size be?
        for i in impulse_response.iter_mut() {
            *i = rng.gen_range(-1.0..=1.0);
        }
        fast_convolver.process(&input_signal, &mut output_signal);
        for i in 0..input_signal.len() {
            assert_eq!(input_signal[i], output_signal[i]);
        }
    }

    fn test_blocksize_frequency_domain() {
        let mut rng = rand::thread_rng();
        let mut input_signal: Vec<f32> = vec![0f32; 10000];
        input_signal[3] = 1.0;
        let mut output_signal: Vec<f32> = vec![0f32; 10000];
        let mut impulse_response: Vec<f32> = vec![0f32; 51];
        let block_sizes = [1, 13, 1023, 2048,1,17, 5000, 1897];
        let mut fast_convolver = FastConvolver::new(&impulse_response, ConvolutionMode::FrequencyDomain { block_size: (1024) }); // What should block size be?

        for block_size in block_sizes {
            for block_index in 0..=input_signal.len() / block_size {
                let block_start_index = block_index * block_size;
                let block_end_index = min(input_signal.len() - 1, (block_index + 1) * block_size);
                // println!("block_start_index: {}, block_end_index: {}", block_start_index, block_end_index);
                fast_convolver.process(&input_signal[block_start_index ..block_end_index], &mut output_signal[block_start_index..block_end_index]);
            }
            for i in 0..input_signal.len() {
                assert_eq!(input_signal[i], output_signal[i]);
            }
        }
    }
}
