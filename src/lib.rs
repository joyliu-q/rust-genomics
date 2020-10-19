use rand::Rng;
use std::io::{Error, ErrorKind};

// Reading fasta/gff3 files
use std::io::{Write, BufReader, BufRead};
use std::fs;
use std::fmt;

use std::time::{Instant};

pub const SEQUENCE_LEN: i32 = 20;
pub const NUCLEOTIDE: [char;4] = ['A', 'T', 'C', 'G'];

pub struct FASTA {
    pub name: String,
    pub content: Vec<FastaRecord>,
}
impl FASTA {
    pub fn new(name: String, content: Vec<FastaRecord>) -> FASTA {
        FASTA{name, content}
    }
}
impl fmt::Display for FASTA {
    //TODO
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "Name: {}", self.name)
    }
}

pub struct FastaRecord {
    pub header: String,
    pub sequence: String,
}
impl FastaRecord {
    pub fn new(header: String, sequence: String) -> FastaRecord {
        FastaRecord{header, sequence}
    }
}
impl fmt::Display for FastaRecord {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "Header: {}\nSequence: {}", self.header, self.sequence)
    }    
}

//TODO: RIGHT NOW NOT USING SEQUENCE TYPE
pub struct Sequence {
    pub seq: String,
}
impl Sequence {
    pub fn new(seq: String) -> Sequence {
        Sequence{seq}
    }
    // TODO: check sequence validity (e.g. is ATCG)
    pub fn check(&self) -> bool {
        true
    }
    pub fn find_at_index(&self, index: usize) -> Result<&str, Error> {
        match self.seq.get(index..index+1) {
            Some(n) => {
                return Ok(n);
            },
            None => {
                return Err(Error::new(ErrorKind::Other, format!["Nucleotide not found at {}.", &index]));
            }
        }
    }
    /// TODO: There are 3 reading frames for codons: let's call them shift 0, 1, and 2.
    /// For example, for sequence ATCAGGCAT
    ///     Reading Frame 0: [ATC][AGG][CAT]...
    ///     Reading Frame 1: A[TCA][GGC]AT...
    ///     Reading Frame 2: AT[CAG][GCA]T...
    /// This function takes the current sequence and returns a vector of codons at these 3 frames
    pub fn pack_into_codons(&self) -> Vec<Vec<char>>{
        let codon_list = vec![Vec::new(), Vec::new(), Vec::new()];
        return codon_list;
    }
    pub fn gen_random_seq(len: i32) -> Sequence {    
        let mut rng = rand::thread_rng();
        let mut seq: String = String::new();
    
        for _ in 0..len {
            let i = rng.gen_range(0, 4);
            seq.push(NUCLEOTIDE[i]);
        }
        Sequence{seq}
    }
}

pub fn read_fasta(path: &str) {
    //let now = Instant::now();

    // Speedy Way
    let data = fs::read_to_string(path).unwrap();
    let data: Vec<&str> = data.split('>').collect();

    let mut records: Vec<FastaRecord> = Vec::new();

    for entry in data {
        if entry.is_empty() {continue}
        let mut entry: Vec<&str> = entry.split("\n").collect();
        let header = entry.remove(0);
        let mut sequence: String = entry.into_iter().collect();
        sequence = sequence.replace("\n", "");
        sequence = sequence.replace("\r", "");

        records.push(FastaRecord::new(header.to_string(), sequence));
    }

    // More controlled way
    /*let file = fs::File::open(path).expect("path to file not found");
    let reader = BufReader::new(file);

    let mut data = Vec::new();
    let mut temp_header = "".to_string();
    let mut temp_seq = "".to_string();

    for (index, line) in reader.lines().enumerate() {
        let line = line.unwrap();
        if line.is_empty() {continue}
        if line.contains('>') {
            // push all prev record, don't push for 1st line
            if index > 0 {
                data.push(FastaRecord::new(temp_header, temp_seq.to_owned()));
            }
            // start new record
            temp_header = line;
            continue;
        }
        temp_seq.push_str(&line);
    }
    // push final record
    data.push(FastaRecord::new(temp_header, temp_seq.to_owned()));*/

    //println!("{}", now.elapsed().subsec_nanos());
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn gen_correct_length() {
        let sequence = Sequence::gen_random_seq(1000);
        assert_eq!(sequence.seq.len(), 1000);
    }

    #[test]
    fn find_correct_index() {
        let sequence = Sequence::gen_random_seq(1000);
        assert![sequence.find_at_index(10).is_ok()];
        assert![sequence.find_at_index(2000).is_err()];
    }
}