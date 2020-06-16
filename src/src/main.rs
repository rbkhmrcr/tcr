#![allow(unused_variables)]
#![allow(dead_code)]
#![allow(non_snake_case)]
#![feature(test)]
extern crate bn;
extern crate rand;
extern crate test;

use bn::*;
use rand::thread_rng;

pub struct Witness {
    x0: Fr,
    x1: Fr,
    x2: Fr,
    b: Fr, // though has to be restricted to 0 / 1
}

pub fn genwitness() -> Witness {
    let rng = &mut thread_rng();
    let x0 = Fr::random(rng);
    let x1 = Fr::random(rng);
    let x2 = Fr::random(rng);
    let b = Fr::random(rng);
    Witness { x0, x1, x2, b }
}

pub struct Statement {
    g0: G1,
    g1: G1,
    h0: G1,
    h1: G1,
    y: G1,
}

pub fn genstatement() -> Statement {
    let rng = &mut thread_rng();
    let g0 = G1::random(rng);
    let g1 = G1::random(rng);
    let h0 = G1::random(rng);
    let h1 = G1::random(rng);
    let y = G1::random(rng);
    Statement { g0, g1, h0, h1, y }
}

pub struct Proof {
    R0: G1,
    R1: G1,
    R2: G1,
    R3: G1,
    S0: G1,
    S1: G1,
    T: G1,
    d0: Fr,
    d1: Fr,
    u: Fr,
    v: Fr,
    w: Fr,
}

pub struct Commitments {
    c0: G1,
    c1: G1,
    c2: G1,
    c3: G1,
    c4: G1,
}

pub struct Computed {
    A0: G1,
    A1: G1,
    A2: G1,
    A3: G1,
    A4: G1,
    A5: G1,
}

pub fn make_commitments(c: Commitments) {
    let A0 = c.c0;
    let A1 = c.c1;
    let A2 = c.c2;
    let A3 = c.c1 - c.c2;
    let A4 = c.c1 - c.c2 + c.c4;
    let A5 = c.c3;
}

pub fn deposit(rng: &mut rand::ThreadRng, amt: Fr, ps: &Statement) -> (G1, G1, Fr) {
    let r = Fr::random(rng);
    let C = ps.g1 * amt + ps.h0 * r;
    // A = h0^a, c, z = c amt + a.
    let a = Fr::random(rng);
    let A = ps.h0 * a;
    let c = Fr::random(rng);
    let z = c * amt + r;
    (C, A, z)
}

pub fn vote1(vote: Fr, wgt: Fr, ps: &Statement) -> (G1, G1, G1, Fr) {
    let rng = &mut thread_rng();
    let x = Fr::random(rng);
    let c0 = ps.g0 * x;
    let c1 = ps.g1 * vote + ps.h0 * x;
    // A = g^a, c, z = cx + a
    let a = Fr::random(rng);
    let A = ps.g0 * a;
    let c = Fr::random(rng);
    let z = c * x + a;
    (c0, c1, A, z)
}

pub fn vote2(vote: Fr, x: Fr, ps: &Statement) -> (G1, G1, G1) {
    let rng = &mut thread_rng();
    let witness = genwitness();
    let s = Fr::random(rng);
    let Y = G1::random(rng);
    let c2 = ps.g1 * vote + Y * x;
    let c3 = ps.g0 * s;
    let c4 = (-ps.h0 + Y) * x + ps.h1 * s;
    let proof = genproof(&witness, ps);
    (c2, c3, c4)
}

pub fn genproof(pw: &Witness, ps: &Statement) -> Proof {
    let rng = &mut thread_rng();
    let r0 = Fr::random(rng);
    let r1 = Fr::random(rng);
    let s0 = Fr::random(rng);
    let s1 = Fr::random(rng);
    let s2 = Fr::random(rng);
    let R0 = ps.g0 * r0;
    let R1 = (ps.h0 - ps.y) * r0;
    let R2 = ps.h1 * r1;
    let R3 = ps.g0 * r1;
    let S0 = ps.g1 * s0 + ps.h0 * s1;
    let S1 = ps.g1 * s0 + ps.y * s1;
    let T = ps.g1 * (s0 * pw.b) + ps.h0 * s2;

    let a0 = Fr::random(rng);
    let a1 = Fr::random(rng);
    let a2 = Fr::random(rng);

    let d0 = r0 + a0 * pw.x0;
    let d1 = r1 + a1 * pw.x2;
    let u = s0 + a2 * pw.b;
    let v = s1 + a2 * pw.x1;
    let w = s2 + pw.x1 * (a2 - u);

    Proof {
        R0,
        R1,
        R2,
        R3,
        S0,
        S1,
        T,
        d0,
        d1,
        u,
        v,
        w,
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    use test::Bencher;

    #[bench]
    fn bench_deposit(b: &mut Bencher) {
        let rng = &mut thread_rng();
        let dep = Fr::random(rng);
        let ps = genstatement();
        b.iter(|| deposit(rng, dep, &ps));
    }

    #[bench]
    fn bench_update(b: &mut Bencher) {
        let rng = &mut thread_rng();
        let update = Fr::random(rng);
        let ps = genstatement();
        b.iter(|| deposit(rng, update, &ps));
    }

    #[bench]
    fn bench_vote1(b: &mut Bencher) {
        let ps = genstatement();
        let vote = Fr::one();
        let rng = &mut thread_rng();
        let weight = Fr::random(rng);
        b.iter(|| vote1(vote, weight, &ps));
    }

    #[bench]
    fn bench_vote2(b: &mut Bencher) {
        let pw = genwitness();
        let ps = genstatement();
        b.iter(|| genproof(&pw, &ps));
    }
}
