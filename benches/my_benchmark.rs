use criterion::{black_box, criterion_group, criterion_main, Criterion};
use rust_genomics::{Sequence};

pub fn criterion_benchmark(c: &mut Criterion) {
    c.bench_function("gen seq", |b| b.iter(|| {
        let long_sequence = Sequence::gen_random_seq(black_box(10000));
    }));
    c.bench_function("normal lorf", |b| b.iter(|| {
        let long_sequence = Sequence::gen_random_seq(black_box(10000));
        let lorf = long_sequence.find_lorf();
    }));
    c.bench_function("concurrent lorf", |b| b.iter(|| {
        let long_sequence = Sequence::gen_random_seq(black_box(10000));
        let lorf_concurrent = long_sequence.concurrent_find_lorf();
    }));
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
