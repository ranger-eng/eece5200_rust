use ndarray::*;
use ndarray_linalg::*;

fn main() {
    let a: Array2<f64> = array![
        [-1.0,  0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 12.0]
    ];
    let (eigs, vecs) = a.eig().unwrap();

    println!["a:{}\neigs:{}\nvecs:{}",a, eigs, vecs];

}
