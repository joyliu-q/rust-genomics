use rust_genomics::{Sequence, FASTA};

fn main() {
    let file_path = "data/haha-1.fasta";

    // Generate an instance of the FASTA struct by reading data from the given file
    let fasta_result = FASTA::rayon_read_fasta(file_path);

    // The FASTA instance has the Display trait, allowing it to be visualized in a clean fashion
    println!("{}", fasta_result);

    // 
}