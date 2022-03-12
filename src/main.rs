pub mod port3;

fn main() {
    // Init
    let n_h: usize = 4;
    let n_x: usize = 10;

    let mut x_array: Vec<f32> = vec![];
    for i in 0..n_x {
        if i <= 10 {x_array.push(1.0)}
        else {x_array.push(0.0)};
    }
    let h_array: Vec<f32> = vec![1.0, 0.5, 0.5*0.5, 0.5*0.5*0.5];
    let w_array: Vec<f32> = vec![1.1, 0.4, 0.5*0.5, 0.5*0.5*0.5];

    // System response
    let y_array = convolution(n_h, n_x, x_array.as_slice(), h_array.as_slice());
    
    // 
    let mut y_hat_array: Vec<f32> = vec![0.0; y_array.len()];

    for i in 0..y_hat_array.len() {
        for j in 0..w_array.len() {
            if i+j<x_array.len() {
                y_hat_array[i] = y_hat_array[i] + w_array[j]*x_array[i + j];
            }
        }
    }

    // test gess
    let A: Vec<f32> = vec![1.0,0.0,1.0,0.0];
    let B: Vec<f32> = vec![24.0,-19.0];
    let a_n: i32 = 2;
    let b_m: i32 = 1;
    
    println!("{:?}",y_hat_array);
}

fn convolution(n_h: usize, n_x: usize, x: & [f32], h: & [f32]) -> Vec<f32> {
    let n_y: usize = n_h + n_x -1;
    let mut x_stack_array: Vec<f32> = vec![0.0; n_h];
    let mut y_array: Vec<f32> = vec![0.0; n_y];

    for i in 0..n_y {

        if i < n_x {x_stack_array[i%n_h] = x[i]}
        else {x_stack_array[i%n_h] = 0.0};
        
        for j in 0..n_h {
            y_array[i] = y_array[i] + h[(n_h-j+i)%n_h]*x_stack_array[j%n_h];
        }
        print!("y:{}\n", y_array[i]);
    }

    return y_array;
}

