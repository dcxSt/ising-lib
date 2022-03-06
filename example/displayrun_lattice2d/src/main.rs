use ising_lib::lattice2d::Lattice2d;
use std::{thread, time};

fn main() {
    let seven_millis = time::Duration::from_millis(7);  // wait 7 millis to slow simulation
    let mut lattice = Lattice2d::new_basic([45,145]);   // create a lattice
    lattice.disp_terminal();
    for _ in 0..1200 {                  // 1200 frames 
        for _ in 0..2000 {
            lattice.update();           // update lattice 2000 times between each frame
        }
        thread::sleep(seven_millis);    // wait 7 millis for smooth flip-card video display
        lattice.disp_terminal();        // display the lattice in terminal
    }
}
