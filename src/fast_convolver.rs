struct FastConvolver {
    mode: ConvolutionMode,
    impulse_response: Vec<f32>,
    previous_block: Vec<f32>
}

#[derive(Debug, Clone, Copy)]
pub enum ConvolutionMode {
    TimeDomain,
    FrequencyDomain { block_size: usize },
}

impl FastConvolver {
    pub fn new(impulse_response: &[f32], mode: ConvolutionMode) -> Self {
        FastConvolver {
            mode,
            impulse_response: impulse_response.to_vec(),
            previous_block: vec![0f32; impulse_response.len()]
        }
    }

    pub fn reset(&mut self) {
        todo!("implement")
    }

    pub fn process(&mut self, input: &[f32], output: &mut [f32]) {
        if matches!(self.mode, ConvolutionMode::TimeDomain) {
            let iL = self.impulse_response.len();
            let pbL = self.previous_block.len();
            for i in 0..input.len() - iL {
                let mut accum = 0.0;
                for j in 0..self.impulse_response.len() {
                    if i + j < iL - 1 {
                        accum += self.impulse_response[j] * self.previous_block[pbL - (iL - 2) + i + j];
                    } else {
                        accum += self.impulse_response[j] * input[i + j - iL + 1];
                    }
                }
                output[i] = accum;
            }
            self.previous_block = output.to_vec();
        }
    }

    pub fn flush(&mut self, output: &mut [f32]) {
        let iL = self.impulse_response.len();
        for i in 0..self.previous_block.len() - iL {
            let mut accum = 0.0;
            for j in 0..self.impulse_response.len() {
                accum += self.previous_block[i + j] * self.impulse_response[j];
            }
            output[i] = accum;
        }
    }

    pub fn get_flush_length(&mut self) -> usize {
        return self.previous_block.len() - (self.impulse_response.len() - 1);
    }

}

