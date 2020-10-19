use rust_genomics::{read_fasta};
use std::fs;

use std::time::{Instant};

fn main() {
    read_fasta("data/sars_cov2_snip.fasta");
}