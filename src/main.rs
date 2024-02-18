use std::{f32::consts::PI, fs::File, io::Write};

use comb_filter::CombFilter;

mod comb_filter;
mod ring_buffer;

fn show_info() {
    eprintln!("MUSI-6106 Assignment Executable");
    eprintln!("(c) 2024 Joseph Cleveland");
}

fn main() {
   show_info();

    // Parse command line arguments
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 3 {
        eprintln!("Usage: {} <input wave filename> <output text filename>", args[0]);
        return
    }

    // Open the input wave file
    let mut reader = hound::WavReader::open(&args[1]).unwrap();
    let spec = reader.spec();
    let channels: usize = spec.channels as usize;

    let block_size = 1024;
    let mut filter = CombFilter::new(comb_filter::FilterType::FIR, 1f32, spec.sample_rate as f32, channels as usize);
    filter.set_param(comb_filter::FilterParam::Delay, args[3].parse().unwrap());
    filter.set_param(comb_filter::FilterParam::Gain, args[4].parse().unwrap());

    // Read audio data and write it to the output text file (one column per channel)
    let mut input_buffer: Vec<Vec<f32>> = Vec::new();
    for i in 0..channels {
        input_buffer.push(vec![0f32; block_size]);
    }
    let mut out = File::create(&args[2]).expect("Unable to create file");
    for (i, sample) in reader.samples::<i16>().enumerate() {
        let sample = sample.unwrap() as f32 / (1 << 15) as f32;
        input_buffer[i % channels as usize].push(sample);
        if (i + 1) / channels % block_size == 0 {
            let mut input_slices: Vec<&[f32]> = Vec::new();
            let mut output_buffer: Vec<Vec<f32>> = vec![vec![0f32; block_size]; channels];

            for c in 0..channels {
                input_slices.push(&input_buffer[c][i+1-block_size..i]);
            }
            let mut slice2d: Vec<&mut [f32]> = output_buffer.iter_mut().map(|v: &mut Vec<f32>| v.as_mut_slice()).collect();
            filter.process(&input_slices, &mut slice2d);

            for j in 0..channels {
                for k in 0..block_size {
                    write!(out, "{}{}", output_buffer[j][k], if i % channels as usize == (channels - 1).into() { "\n" } else { " " }).unwrap();
                }
            }
        }
    }
}

fn sine_gen(freq: f32, Fs: f32, length: usize) -> Vec<f32> {
    let mut output = vec![0f32; length];
    for i in 0..length {
        output[i] = f32::sin((i as f32) * 2.0 * freq * PI / Fs);
    }
    return output;
}

#[test]
fn fir_output_zero() {
    const FS: f32 = 44100.0;
    let mut filter = CombFilter::new(comb_filter::FilterType::FIR, 1.0, FS, 1);
    filter.set_param(comb_filter::FilterParam::Delay, 1.0 / 200.0);
    filter.set_param(comb_filter::FilterParam::Gain, 1.0);
    let input = sine_gen(100.0, FS, 2*FS as usize);
    let mut output = [0f32; 2*FS as usize];

    filter.process([input.as_slice()].as_slice(), [output.as_mut_slice()].as_mut_slice());

    for i in FS as usize..2*FS as usize{
        assert!(f32::abs(output[i]) < 0.01);
    }
}

#[test]
fn iir_magnitude_increase() {
    const FS: f32 = 44100.0;
    let mut filter = CombFilter::new(comb_filter::FilterType::IIR, 1.0, FS, 1);
    filter.set_param(comb_filter::FilterParam::Delay, 1.0 / 100.0);
    filter.set_param(comb_filter::FilterParam::Gain, 0.5);
    let input = sine_gen(100.0, FS, 4*FS as usize);
    let mut output = [0f32; 4*FS as usize];

    filter.process([input.as_slice()].as_slice(), [output.as_mut_slice()].as_mut_slice());
    
    let mut max = 0.0;
    for i in FS as usize..4*FS as usize{
        if f32::abs(output[i]) > max {
            max = f32::abs(output[i])
        }
    }
    assert!(f32::abs(max - 2.0) <= 0.001);
}

#[test]
fn vary_block_size() {
    const FS: f32 = 44100.0;
    let mut filter = CombFilter::new(comb_filter::FilterType::FIR, 1.0, FS, 1);
    filter.set_param(comb_filter::FilterParam::Delay, 1.0 / 100.0);
    filter.set_param(comb_filter::FilterParam::Gain, 1.0);
    let input = sine_gen(100.0, FS, 10*FS as usize);
    let mut output = [0f32; 10*FS as usize];

    for i in 0..10 {
        let block_size = if i % 2 == 0 { 6300} else {3150};
        for block_index in (i*FS as usize..(i+1)*FS as usize).step_by(block_size) {
            filter.process([&input[block_index..block_index+block_size]].as_slice(), [&mut output[block_index..block_index+block_size]].as_mut_slice());
        }
    }

    for i in FS as usize..10*FS as usize{
        assert!(f32::abs(output[i] - 2.0 * input[i]) < 0.001);
    }
}

#[test]
fn test_reset() {
    const FS: f32 = 44100.0;
    let mut filter = CombFilter::new(comb_filter::FilterType::IIR, 1.0, FS, 1);
    filter.set_param(comb_filter::FilterParam::Delay, 1.0 / 100.0);
    filter.set_param(comb_filter::FilterParam::Gain, 0.5);
    let input = sine_gen(100.0, FS, 4*FS as usize);
    let mut output = [0f32; 4*FS as usize];

    filter.process([input.as_slice()].as_slice(), [output.as_mut_slice()].as_mut_slice());

    filter.reset();

    for i in 0..filter.buffer.capacity() {
        assert!(f32::abs(filter.buffer.get(0)) <= 0.0001);
    }
}