use criterion::{black_box, criterion_group, criterion_main, Criterion, BenchmarkId};
use rust_genomics::{Sequence};

pub fn criterion_benchmark(c: &mut Criterion) {
    /*c.bench_function("gen seq", |b| b.iter(|| {
        let long_sequence = Sequence::gen_random_seq(black_box(100000));
    }));*/

    /*let long_sequence = Sequence::gen_random_seq(black_box(100000));

    c.bench_function("normal lorf", |b| b.iter(|| {
        long_sequence.find_lorf();
    }));
    c.bench_function("concurrent lorf", |b| b.iter(|| {
        long_sequence.concurrent_find_lorf();
    }));*/

    let sequences = [
        Sequence::gen_random_seq(black_box(100)), // 32.756 us for concurrent
        Sequence::gen_random_seq(black_box(1000)), //40.622 us for concurrent
        Sequence::gen_random_seq(black_box(10000)), //107.72 us for concurrent
        Sequence::gen_random_seq(black_box(100000)), //1.8189 ms for concurrent
        Sequence::gen_random_seq(black_box(100000)), //1.8598 ms for concurrent
    ];

    let mut norm_group = c.benchmark_group("concurrent lorf group");
    for seq in sequences.iter() {
        norm_group.bench_with_input(BenchmarkId::from_parameter(seq), seq, |b, seq| {
            b.iter(|| {
                seq.concurrent_find_lorf();
            });
        });
    }
    norm_group.finish();

    let mut con_group = c.benchmark_group("concurrent lorf group");
    for seq in sequences.iter() {
        con_group.bench_with_input(BenchmarkId::from_parameter(seq), seq, |b, seq| {
            b.iter(|| {
                seq.concurrent_find_lorf();
            });
        });
    }
    con_group.finish();
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
