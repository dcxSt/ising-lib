use ising_lib::lattice2d::Lattice2d;
use std::{thread, time};

fn main() {
    let sixteen_millis = time::Duration::from_millis(7); // wait 16 millis for 60fps
    let mut lattice = Lattice2d::new_basic([35,155]);
    lattice.disp_terminal();
    for _ in 0..1500 {
        for _ in 0..2000 {
            lattice.update();
        }
        thread::sleep(sixteen_millis); // wait 16 millis for 60 fps
        lattice.disp_terminal();
    }
}
