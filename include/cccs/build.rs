// extern crate cc;
// extern crate cmake;

// use cmake::Config;
use std::env;

fn main() {
    // let dst = Config::new("vendor/spoa")
    //        .define("CMAKE_BUILD_TYPE","Release")
    //        .build();

    // println!("cargo:rustc-link-search=native={}", dst.display());
    // println!("cargo:rustc-link-lib=static=spoa");

    let out_dir = env::var("OUT_DIR").unwrap();
    println!("cargo:rustc-flags=-L {}/lib64/ -L {}/lib/", &out_dir, &out_dir);
    
    // cc::Build::new()
    //     .cpp(true)
    //     .shared_flag(false)
    //     .static_flag(true)
    //     .flag_if_supported("-O3")
    //     .flag_if_supported("-D_GNU_SOURCE")
    //     .flag_if_supported("-Wall")
    //     .flag_if_supported("-fPIC")
    //     .flag_if_supported("-std=c++11")
    //     .flag_if_supported("-Ivendor/spoa/include")
    //     .flag_if_supported(&format!("-L{}/lib64 -L{}/lib", &out_dir, &out_dir))
    //     .flag_if_supported("-lspoa")
    //     .file("src/poa.cpp")
    //     .compile("poa");
}
