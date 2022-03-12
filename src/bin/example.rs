use ndarray::*;
use ndarray_linalg::*;

fn main() {
let a: Array2<f64> = array![
    [-1.01,  0.86, -4.60,  3.31, -4.81],
    [ 3.98,  0.53, -7.04,  5.29,  3.55],
    [ 3.30,  8.26, -3.89,  8.20, -1.51],
    [ 4.43,  4.96, -7.66, -7.33,  6.18],
    [ 7.31, -6.43, -6.16,  2.47,  5.58],
];
let (eigs, vecs) = a.eig().unwrap();

let a = a.map(|v| v.as_c());
for (&e, vec) in eigs.iter().zip(vecs.axis_iter(Axis(1))) {
    let ev = vec.map(|v| v * e);
    let av = a.dot(&vec);
    assert_close_l2!(&av, &ev, 1e-5);
}
}