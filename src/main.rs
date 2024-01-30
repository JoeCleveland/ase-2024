use std::{fs::File, io::Write};

fn show_info() {
    eprintln!("MUSI-6106 Assignment Executable");
    eprintln!("(c) 2024 Stephen Garrett & Ian Clester");
}

fn main() -> std::io::Result<()> {
   show_info();

    // Parse command line arguments
    // First argument is input .wav file, second argument is output text file.
    let args: Vec<String> = std::env::args().collect();
    let in_file = &args[1];
    let out_file = &args[2];
    dbg!(&in_file);

    // Open the input wave file and determine number of channels
    let mut wav_file = hound::WavReader::open(in_file).unwrap();
    let mut txt_file = File::create(out_file).unwrap();
    let spec = wav_file.spec();
    dbg!((&wav_file).len());
    let mut channel = 0;
    for sample in wav_file.samples::<i16>() {
        let fval = sample.unwrap() as f32;
        if channel == spec.channels {
            writeln!(&mut txt_file)?;
            channel = 0;
        }
        write!(&mut txt_file, "{} ", fval / 65535.0)?;
        channel += 1;
    }

    Ok(())
}
