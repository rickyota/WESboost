use std::process::Command;
use walkdir::WalkDir;

fn main() {
    // makefile here

    let out_dir_c: String = std::env::var_os("OUT_DIR")
        .unwrap()
        .to_string_lossy()
        .to_string()
        + "/c/";
    //let out_dir_c: String = format!("{:?}/c/", std::env::var_os("OUT_DIR").unwrap());

    // compile luacall library: OS dependant
    if cfg!(target_os = "linux") {
        // makefile is using a special env variable
        // ex.  "${root}/projects_rust/target/release/build/pgenlib_rust-08f8f451ee789acf/out"
        //panic!("DESTDIR, {:?}", std::env::var_os("DEST_DIR").unwrap());
        // panic!("out_dir_c, {:?}", &out_dir_c);
        Command::new("sh")
            .arg("-c")
            .arg(&format!(
                "cd ./lib/pgenlib && make zstd RUSTDESTDIR={:?}",
                &out_dir_c
            ))
            .status()
            .expect("failed to make!");
        // Command::new("sh")
        //     .arg("-c")
        //     .arg("cd ./lib/pgenlib && make zstd")
        //     .arg(&format!("RUSTDESTDIR={:?}", &out_dir_c))
        //     .status()
        //     .expect("failed to make!");
        // &format!("RUSTDESTDIR={:?}", std::env::var_os("DEST_DIR").unwrap()),
        //.args(&["-f", "src/lua/luacall_linux.mak"])
        // .args(&[ "-c",
        //     &format!(
        //         "cd ./lib/pgenlib && make zstd RUSTDESTDIR={:?}",
        //         std::env::var_os("DEST_DIR").unwrap()
        //     ),
        // ])
    } else {
        panic!("Unsupported OS");
    }

    //let dirs = ["./src/lib/pgenlib/"];
    let dirs = ["./lib/pgenlib/"];
    // exclude src/unuse
    //let dirs = [
    //    "./lib/pgenlib/src/",
    //    "./lib/pgenlib/src/include/",
    //    "./lib/pgenlib/src/simde/",
    //    "./lib/pgenlib/src/simde/x86",
    //    "./lib/pgenlib/src/simde/x86/avx512",
    //];

    // all files with .cc or .cpp
    //let cpps: Vec<String> = vec!["./lib/pgenlib/pgenlibr_wrapc.cpp".to_string()];

    let cpps: Vec<String> = dirs
        .map(|dir| {
            WalkDir::new(dir)
                .into_iter()
                .map(|x| x.unwrap().path().display().to_string())
        })
        .into_iter()
        .flatten()
        .filter(|x| x.ends_with(".cpp") || x.ends_with(".cc"))
        .collect();

    // run `make zstd` in ./lib/pgenlib first
    let objs: Vec<String> = [&out_dir_c]
        .map(|dir| {
            WalkDir::new(dir)
                .into_iter()
                .map(|x| x.unwrap().path().display().to_string())
        })
        .into_iter()
        .flatten()
        .filter(|x| x.ends_with(".o"))
        .collect();
    // let objs: Vec<String> = dirs
    //     .map(|dir| {
    //         WalkDir::new(dir)
    //             .into_iter()
    //             .map(|x| x.unwrap().path().display().to_string())
    //     })
    //     .into_iter()
    //     .flatten()
    //     .filter(|x| x.ends_with(".o"))
    //     .collect();

    // cannot print in build.rs
    //println!("cpps: {:?}", cpps);

    // how to add openmp?
    // https://users.rust-lang.org/t/binding-openmp-c-function/40196/4

    // should be .rs with #[cxx::bridge]
    // see pgenlibr-src_ori/Makefile
    cxx_build::bridge("src/lib.rs")
        .files(&cpps)
        .objects(&objs)
        .include("lib/pgenlib/zstd/lib")
        .include("lib/pgenlib/libdeflate")
        .flag_if_supported("-fopenmp")
        .flag_if_supported("-g")
        .flag_if_supported("-DIGNORE_BUNDLED_ZSTD")
        .flag_if_supported("-DZSTD_MULTITHREAD")
        .compile("pgenlib-bridge");

    // not necessary here; add in Makefile
    // .flag_if_supported("-fPIE")
    // added -fPIE due to error below
    // = note: /usr/bin/ld: ${root}/projects_rust/target/release/deps/libpgenlib_rust-ca6021eb115fcc32.rlib(deflate_decompress.o): relocation R_X86_64_32 against `.text' can not be used when making a PIE object; recompile with -fPIC
    //       /usr/bin/ld: final link failed: Nonrepresentable section on output
    //       collect2: error: ld returned 1 exit status

    // .include("lib/pgenlib/include/")

    // ng: .flag_if_supported("-lz")
    println!("cargo:rustc-link-lib=z");
    println!("cargo:rustc-link-lib=zstd");
    println!("cargo:rustc-link-lib=pthread");
    // TODO
    // .flag_if_supported("-llibdeflate")
    // .flag_if_supported("-lm")
    // .flag_if_supported("-ldl")

    // .flag_if_supported("-fpermissive") // ignore type conversion error in zstd
    //.flag_if_supported("-fpthread")
    //.flag_if_supported("-fzstd")

    // This requires -lzstd ??
    // .flag_if_supported("-DIGNORE_BUNDLED_ZSTD")
    // .flag_if_supported("-DZSTD_MULTITHREAD")

    //.flag_if_supported("-static")

    // seems not necessary
    //.flag_if_supported("-std=c++")
    //.flag_if_supported("-std=c++11")

    // seems not necessary
    //.flag_if_supported("-lgomp")

    // arg for rustc compile
    // This is same as
    // $ export RUSTFLAGS='-C link-args=-fopenmp'
    // Better use build.rs to make build.sh simple
    //
    println!("cargo:rustc-link-arg=-fopenmp");
    // seems not necessary
    // [ref](https://github.com/rust-lang/cc-rs/issues/266)
    //println!("cargo:rustc-link-lib=gomp");

    println!("cargo:rerun-if-changed=src/lib.rs");
    println!("cargo:rerun-if-changed=lib/pgenlib/Makefile");
    println!("cargo:rerun-if-changed=lib/pgenlib/Makefile.src");
    // all cpp and h
    dirs.map(|dir| {
        WalkDir::new(dir)
            .into_iter()
            .map(|x| x.unwrap().path().display().to_string())
    })
    .into_iter()
    .flatten()
    .filter(|x| {
        x.ends_with(".cpp")
            || x.ends_with(".cc")
            || x.ends_with(".hpp")
            || x.ends_with(".h")
            || x.ends_with(".o")
    })
    .for_each(|x| println!("cargo:rerun-if-changed={:?}", x));
}
