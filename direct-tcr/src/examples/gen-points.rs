extern crate bn;
extern crate rand;

use bn::*;

fn main() {
    let rng = &mut rand::thread_rng();

    // Construct private keys
    let alice_sk = Fr::random(rng);

    print!("g0 {:?}\n", G1::random(rng));
    print!("g1 {:?}\n", G1::random(rng));
    print!("h0 {:?}\n", G1::random(rng));
    print!("h1 {:?}\n", G1::random(rng));
    print!("y {:?}\n", G1::random(rng));

}
