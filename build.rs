use cmake;


fn main() {
    let dst = cmake::build("./modules");
    println!("cargo:rustc-link-search=native={}", dst.display());
    println!("cargo:rustc-link-lib=static=port3");
    println!("cargo:rustc-link-arg-bins=-L/usr/local/Cellar/gcc/11.2.0_3/lib/gcc/11");
    println!("cargo:rustc-link-arg-bins=-lgfortran");
    println!("cargo:rustc-link-arg-bins=-lm");
    println!("cargo:rerun-if-changed=build.rs");
}


