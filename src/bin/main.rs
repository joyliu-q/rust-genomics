use rust_genomics::{read_fasta, Sequence};
use std::fs;

fn main() {
    let seq = Sequence::gen_random_seq(1000);
    seq.pack_into_codons();
    //read_fasta("data/sars_cov2_snip.fasta");
}