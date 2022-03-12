use std::os::raw::c_int;
use std::os::raw::c_float;

extern {
    fn gess_(N: *mut c_int, A: *mut c_float, IA: *mut c_int, B: *mut c_float, IB: *mut c_int, NB: *mut c_int, COND: *mut c_float);
}

pub fn gess(a: Vec<f32>, b: Vec<f32>, mut a_n: i32, mut b_m: i32) -> Vec<f32> {
    let mut a_tmp = a.clone();
    let mut x = b.clone();
    let mut cond: f32 = 0.0;

    println!("a:{:?}\nb:{:?}\na_n:{}\nb_m{}",a,b,a_n,b_m);
    unsafe {
        gess_(&mut a_n as *mut c_int, &mut a_tmp[0] as *mut c_float, &mut a_n as *mut c_int, &mut x[0] as *mut c_float, &mut a_n as *mut c_int, &mut b_m as *mut c_int, &mut cond as *mut c_float);
    }
    return x; 
}

extern {
    pub fn GESS();
}
