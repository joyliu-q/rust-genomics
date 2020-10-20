use rand::Rng;
use std::io::{Error, ErrorKind};

// Reading fasta/gff3 files
use std::fs;
use std::fmt;
use std::cmp;
use std::io::{Write, BufReader, BufRead};
use std::time::{Instant};


pub const SEQUENCE_LEN: i32 = 20;
pub const NUCLEOTIDE: [char;4] = ['A', 'T', 'C', 'G'];

type StartIndex = usize;
type StopIndex = usize; 

/// Longest open reading frame (start, stop). Can have one or multiple
pub enum LORF {
    One((StartIndex, StopIndex)),
    Multiple(Vec<(StartIndex, StopIndex)>),
} 

//TODO: RIGHT NOW NOT USING SEQUENCE TYPE
pub struct Sequence {
    pub seq: String,
}
impl fmt::Display for Sequence {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.seq)
    }
}
impl Sequence {
    pub fn new(seq: String) -> Sequence {
        Sequence{seq}
    }
    // TODO: check sequence validity (e.g. is ATCG)
    fn check(&self) -> bool {
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
    /// Takes the current sequence and returns a vector of codons at 3 reading frames.
    /// # Example
    /// For sequence ATCAGGCAT, there are 3 frames. By calling this function on that sequence, it returns 
    /// ```
    /// [
    ///     ["ATC", "AGG", "CAT"],
    ///     ["A", "TCA", "GGC", "AT"],
    ///     ["AT", "CAG", "GCA", "T"]
    /// ]
    /// ```
    pub fn return_reading_frames(&self) -> Vec<Vec<&str>> {
        let mut reading_frame = vec![Vec::new(), Vec::new(), Vec::new()];
        for i in 0..3 {
            let mut seq = &self.seq[..];
            if i > 0 {reading_frame[i].push(&seq[0..i]);}
            seq.to_string().replace_range(0..i, "");
            while !seq.is_empty() {
                let (codon, remaining_seq) = seq.split_at(cmp::min(3, seq.len()));
                reading_frame[i].push(codon);
                seq = remaining_seq;
            }
        }
        println!("{:?}", &reading_frame);
        reading_frame
    }
    fn return_start_stop_positions(codon_list: Vec<&str>) -> (Vec<usize>, Vec<usize>){
        let mut start_positons = Vec::new();
        let mut stop_positions = Vec::new();

        for (index, codon) in codon_list.into_iter().enumerate() {
            match codon {
                "ATG" => start_positons.push(index),
                "TAA" | "TAG" | "TGA" => stop_positions.push(index),
                _ => continue
            }
        }
        (start_positons, stop_positions)
    }
    // TODO: Find longest Open Reading Frame in sequence
    pub fn find_lorf(&self) -> LORF {
        let reading_frames = self.return_reading_frames();
        for frame in reading_frames {
            let (start_positons, stop_positions) = Sequence::return_start_stop_positions(frame);
            //TODO: find ORF based on start and stop lol
            /*
                ORF:
                1. start_index < stop_index
                2. no stops between current start and stop
                LORF: max length in ORFs
            */
        }
        LORF::One((1,2))
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

        let sequence = Sequence::new(sequence);
        records.push(FastaRecord::new(header.to_string(), sequence));
    }

    // The other way
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

pub struct FastaRecord {
    pub header: String,
    pub sequence: Sequence,
}
impl FastaRecord {
    pub fn new(header: String, sequence: Sequence) -> FastaRecord {
        FastaRecord{header, sequence}
    }
}
impl fmt::Display for FastaRecord {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "Header: {}\nSequence: {}", self.header, self.sequence)
    }    
}

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

    #[test]
    fn test_codon_packing() {
        let seq = Sequence::gen_random_seq(1000);
        seq.return_reading_frames();
    }
}