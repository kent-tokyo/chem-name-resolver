use criterion::{black_box, criterion_group, criterion_main, Criterion};

fn bench_resolve(c: &mut Criterion) {
    c.bench_function("resolve_ethanol", |b| {
        b.iter(|| chem_name_resolver::resolve(black_box("ethanol")))
    });
    c.bench_function("resolve_2chlorobutane", |b| {
        b.iter(|| chem_name_resolver::resolve(black_box("2-chlorobutane")))
    });
    c.bench_function("normalize_greek", |b| {
        b.iter(|| chem_name_resolver::normalize(black_box("α-D-glucose")))
    });
}

criterion_group!(benches, bench_resolve);
criterion_main!(benches);
