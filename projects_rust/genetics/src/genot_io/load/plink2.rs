use crate::genot::prelude::*;
use crate::{alloc, FillMissing, GenotFile, GenotFiles};
use crate::{genot_io, vec};
use crate::{genotype, FillMissingGroup};
use pgenlib::pgenlib_ffi as pgenlib;
use rayon::prelude::*;
use std::collections::HashMap;
use std::ffi::CString;
use std::os::unix::ffi::OsStrExt;
use std::time::Instant;

pub fn generate_genot_snv_file_plink2(
    fin_genot: &GenotFile,
    // multi-allele
    m_in_i: usize,
    // TODO
    // m_snv_i: usize,
    //allelei: usize,
    n: usize,
    use_samples: Option<&[bool]>,
    fill_missing_mode: bool,
) -> GenotSnv {
    //let m_in = genot_io::compute_num_snv_file(fin_genot);

    let (snvs_in, file_snv_allele_idx_in) = genot_io::load_snvs_file(fin_genot, Some(0));
    let m_in = snvs_in.len();

    let n_in = match use_samples {
        None => n,
        Some(usev) => usev.len(),
    };

    let use_samples_idx = convert_use_samples_pgenlib(use_samples, n);
    assert_eq!(use_samples_idx.len(), n);

    let fin_genotype = fgenotype_pgenlib(fin_genot);
    let fin_snv = fsnv_pgenlib(fin_genot);

    let genot_v = load_genot_snv_buf(
        fin_genotype,
        fin_snv,
        m_in,
        m_in_i,
        &file_snv_allele_idx_in,
        n_in,
        &use_samples_idx,
        n,
    );
    //println!("genot_v {:?}", genot_v);

    let mut g_snv = GenotSnv::new_empty(n);
    assign_gsnv_from_genot_i8(&mut g_snv.as_genot_snv_mut_snv(), &genot_v);

    //if !use_missing {
    if fill_missing_mode {
        g_snv.as_genot_snv_mut_snv().fill_missing_mode();
    }

    g_snv
}

fn fgenotype_pgenlib(fin_genot: &GenotFile) -> CString {
    let fin_genotype = fin_genot.genotype_file();
    let fin_str = CString::new(fin_genotype.as_os_str().as_bytes()).unwrap();
    fin_str
}

fn fsnv_pgenlib(fin_genot: &GenotFile) -> CString {
    let fin_genotype = fin_genot.snv_file();
    let fin_str = CString::new(fin_genotype.as_os_str().as_bytes()).unwrap();
    fin_str
}

fn convert_use_samples_pgenlib(use_samples: Option<&[bool]>, n: usize) -> Vec<i32> {
    let use_samples_idx: Vec<i32> = match use_samples {
        None => Vec::from_iter(0..(n as i32)),
        Some(usev) => usev
            .iter()
            .enumerate()
            .filter(|(_, b)| **b)
            .map(|(i, _)| i.try_into().unwrap())
            .collect(),
    };

    use_samples_idx
}

fn load_genot_snv_buf(
    fin_genotype: CString,
    fin_snv: CString,
    m_in: usize,
    m_in_i: usize,
    file_snv_allele_idx_in: &[(usize, usize, usize)],
    n_in: usize,
    use_samples_idx: &[i32],
    //mut use_samples_idx: Vec<i32>,
    n: usize,
    //) -> Vec<f64> {
) -> Vec<i8> {
    //let m_start = mi;
    //let m_end = mi + 1;
    let mut use_snvs = vec![false; m_in];
    use_snvs[m_in_i] = true;

    let mut genot_v = vec![0i8; n];

    load_genot_snvs_extract_buf(
        fin_genotype,
        fin_snv,
        //m_start,
        //m_end,
        &use_snvs,
        &file_snv_allele_idx_in,
        n_in,
        use_samples_idx,
        n,
        &mut genot_v,
    );
    genot_v
}

/* // {0:0, 1:1, 2:2, -3:3}
pub fn assign_pred_from_genot(pred: &mut GenotSnvMut, buf_mi: &[f64]) {
    for (ni, dosage) in buf_mi.iter().enumerate() {
        if *dosage < 0.0 {
            // missing
            pred.set_unchecked(3, ni);
        } else {
            let d = *dosage as u8;
            pred.set_unchecked(d, ni);
        }
    }
} */

pub fn assign_gsnv_from_genot_i8(gsnv: &mut GenotSnvMut, buf_mi: &[i8]) {
    //buf_mi.iter().enumerate().for_each(|(ni, dosage)| {
    buf_mi.iter().enumerate().for_each(|(ni, dosage)| {
        let d = *dosage as u8;
        gsnv.set_unchecked(d, ni);
    });

    // now, missing is 3 in pgenlib
    //for (ni, dosage) in buf_mi.iter().enumerate() {
    //    if *dosage < 0i8 {
    //        // missing
    //        pred.set_unchecked(3, ni);
    //    } else {
    //        let d = *dosage as u8;
    //        pred.set_unchecked(d, ni);
    //    }
    //}
}

// max and min size allocated in pgenlib
// 64 GB
const BUF_SIZE_PGENLIB_MAX: usize = 64 * 1024 * 1024 * 1024;
// 1 GB
const BUF_SIZE_PGENLIB_MIN: usize = 1 * 1024 * 1024 * 1024;

fn mem_buf_limit(m: usize, n: usize, mem: Option<usize>) -> usize {
    // do not check available_memory if mem is provided
    let mem_avail = match mem {
        Some(x) => Some(x),
        None => alloc::get_available_memory(),
    };
    log::debug!("available mem: {:?} bytes", mem_avail);

    let genot_byte = Genot::byte(m, n);

    let mem_buf = match mem_avail {
        Some(x) => {
            log::debug!(
                "genot + min pgenlib vs available mem, {:.3} GB + {:.3} GB vs {:.3} GB",
                alloc::mem_gb(genot_byte),
                alloc::mem_gb(BUF_SIZE_PGENLIB_MIN),
                alloc::mem_gb(x),
            );
            if genot_byte + BUF_SIZE_PGENLIB_MIN > x {
                panic!("Memory insufficient on preparing Genot.")
            }
            x - genot_byte
        }
        None => {
            log::debug!(
                "Could not get available memory; assume there is {} GB available memory.",
                alloc::mem_gb(BUF_SIZE_PGENLIB_MAX)
            );
            BUF_SIZE_PGENLIB_MAX - genot_byte
        }
    };

    mem_buf
}

// load whole is fastest
// TODO: use class to calculate index etc.
pub fn generate_genot_plink2(
    fins_genot: &GenotFiles,
    m_snv: usize,
    // m: usize, // = m_snv + m_group
    n: usize,
    use_snvs: Option<&[bool]>,
    file_snv_allele_idx_in: &[(usize, usize, usize)],
    group_to_m_in: Option<Vec<Vec<usize>>>,
    // group_to_m_in: Option<HashMap<usize, Vec<usize>>>,
    use_samples: Option<&[bool]>,
    fill_missing: Option<FillMissing>,
    // fill_missing_mode: bool,
    fill_missing_group: Option<FillMissingGroup>,
    //use_missing: bool,
    mem: Option<usize>,
    buf_num_snv: Option<usize>, // for testing
) -> Genot {
    let start = Instant::now();

    let m_group = if group_to_m_in.is_some() {
        group_to_m_in.as_ref().unwrap().len()
    } else {
        0
    };

    let m = m_snv + m_group;

    log::debug!("to prepare Genot plink2 m, n: {}, {}", m, n);
    let mem_buf_limit = mem_buf_limit(m, n, mem);

    // TODO: better way
    let use_snvs_v = vec![true; m_snv];
    let use_snvs = match use_snvs {
        Some(x) => x,
        None => &use_snvs_v,
    };

    let n_in = match use_samples {
        None => n,
        Some(usev) => usev.len(),
    };

    let use_samples_idx = convert_use_samples_pgenlib(use_samples, n);
    assert_eq!(use_samples_idx.len(), n);

    let mut g = Genot::new_zeros(m, n);

    let (mut g_snvs, mut g_group) = if group_to_m_in.is_some() {
        let (g_snvs, g_group) = g.split_genot_mut(m_snv);
        assert_eq!(g_snvs.m(), m_snv);
        assert!(g_group.m() == m_group);
        (g_snvs, Some(g_group))
    } else {
        (g.as_genot_mut(), None)
    };

    //if is_split_chrom {
    load_genot_files(
        fins_genot,
        &mut g_snvs,
        &mut g_group,
        use_snvs,
        file_snv_allele_idx_in,
        group_to_m_in,
        use_samples_idx,
        fill_missing_group,
        n_in,
        n,
        Some(mem_buf_limit),
        buf_num_snv,
    );
    //} else {
    //    let fin_genotype = fgenotype_pgenlib(fin_genot, None);

    //    assign_genot(
    //        &mut g.as_genot_mut(),
    //        fin_genotype,
    //        m_in,
    //        n_in,
    //        &use_samples_idx,
    //        n,
    //        use_snvs,
    //        Some(mem_buf_limit),
    //    );
    //}

    // should be done just after loading in perspect of cache locality
    // missing
    super::fill_missing_g_snvs(&mut g_snvs, fill_missing);

    let end = start.elapsed();
    log::info!("It took {} seconds to generate genot.", end.as_secs());

    g
}

fn load_genot_files(
    fins_genot: &GenotFiles,
    g_snvs: &mut GenotMut,
    g_group: &mut Option<GenotMut>,
    use_snvs: &[bool],
    file_snv_allele_idx_in: &[(usize, usize, usize)],
    group_to_m_in: Option<Vec<Vec<usize>>>,
    // group_to_m_in: Option<HashMap<usize, Vec<usize>>>,
    use_samples_idx: Vec<i32>,
    fill_missing_group: Option<FillMissingGroup>,
    n_in: usize,
    n: usize,
    mem_buf_limit: Option<usize>,
    buf_num_snv: Option<usize>, // for testing
) {
    // min of mem_buf, MAX_BED
    // already done in generate_genot_plink2() but do it again for mem_buf_limit=None
    // avoid loading size of pgen file in plink2 since it is slow
    let buf_size_max = buf_size_limit(mem_buf_limit);
    log::debug!("buf_size_max {}", buf_size_max);
    let mut buf: Vec<i8> = vec![0; buf_size_max];

    let mut m_begin = 0;
    let mut m_in_begin = 0;
    for (file_i, fin_genot) in fins_genot.files().iter().enumerate() {
        log::debug!("Loading file {:?}", fin_genot);
        let m_in_file = genot_io::compute_num_snv_allele_file(file_snv_allele_idx_in, file_i);
        //println!("file_snv_allele_idx_in {:?}", &file_snv_allele_idx_in);
        //println!("m_in_file {}", m_in_file);
        //let m_in_chrom = genot_io::compute_num_snv_file(fin_genot);
        let m_in_end = m_in_begin + m_in_file;
        log::debug!("m_in_file {}", m_in_file);
        log::debug!("m_in_end {}", m_in_end);
        let m_file = vec::count_true(&use_snvs[m_in_begin..m_in_end]);
        let m_end = m_begin + m_file;

        let group_to_m_in_file =
            super::group_to_m_in_range(group_to_m_in.as_ref(), m_in_begin, m_in_end);

        let is_group_empty = group_to_m_in_file.as_ref().map_or(true, |x| x.is_empty());

        // skip if no snvs to load
        if m_file == 0 && is_group_empty {
            m_in_begin = m_in_end;
            continue;
        }

        let fin_genotype = fgenotype_pgenlib(fin_genot);
        let fin_snv = fsnv_pgenlib(fin_genot);

        buf = load_genot_file(
            &mut g_snvs.as_genot_snvs_mut(m_begin, m_end),
            fin_genotype,
            fin_snv,
            m_in_file,
            n_in,
            &use_samples_idx,
            n,
            &use_snvs[m_in_begin..m_in_end],
            &file_snv_allele_idx_in[m_in_begin..m_in_end],
            Some(buf),
            Some(buf_size_max),
            g_group.as_mut(),
            group_to_m_in_file.as_ref(),
            fill_missing_group,
            buf_num_snv,
        );

        m_begin = m_end;
        m_in_begin = m_in_end;
    }
    assert_eq!(m_in_begin, use_snvs.len(), "Sth wrong.");
}

// TOFIX: mem should be subtracted by pgenlibr::Load() memory
fn buf_size_limit(mem: Option<usize>) -> usize {
    let buf_size_limit: usize = match mem {
        Some(x) => x.min(BUF_SIZE_PGENLIB_MAX),
        None => {
            log::debug!(
                "Could not get available memory; assume there is {:.3} GB available memory.",
                alloc::mem_gb(BUF_SIZE_PGENLIB_MAX)
            );
            BUF_SIZE_PGENLIB_MAX
        }
    };
    log::debug!("buf_size_limit: {:.3} GB", alloc::mem_gb(buf_size_limit));

    buf_size_limit
}

// m_in_i, buf_i, m_i, group_m_i are all different
fn load_genot_file(
    g_snv_file: &mut GenotMut,
    fin_genotype: CString,
    fin_snv: CString,
    m_in_file: usize,
    n_in: usize,
    use_samples_idx: &[i32],
    n: usize,
    use_snvs: &[bool],
    file_snv_allele_idx_in: &[(usize, usize, usize)],
    buf: Option<Vec<i8>>,
    mem_buf: Option<usize>,
    mut g_group: Option<&mut GenotMut>,
    group_to_m_in_file: Option<&Vec<Vec<usize>>>,
    // group_to_m_in_file: Option<&HashMap<usize, Vec<usize>>>,
    fill_missing_group: Option<FillMissingGroup>,
    buf_num_snv: Option<usize>, // for testing
) -> Vec<i8> {
    let buf_num_snv: usize = match buf_num_snv {
        Some(x) => x,
        None => {
            let buf_size_limit = buf_size_limit(mem_buf);
            // 1 byte (i8) per count in pgenlib
            let byte_per_snv = n * 1;
            // f64: 8 byte per count in pgenlib
            // let byte_per_snv = n * 8;
            let buf_num_snv_limit: usize = buf_size_limit / byte_per_snv;
            buf_num_snv_limit.min(m_in_file)
        }
    };
    // 1 byte (i8) per count in pgenlib
    // let byte_per_snv = n * 1;
    // f64: 8 byte per count in pgenlib
    // let byte_per_snv = n * 8;
    // let buf_num_snv_limit: usize = buf_size_limit / byte_per_snv;
    // let buf_num_snv: usize = buf_num_snv_limit.min(m_in_file);
    let buf_size = buf_num_snv * n;

    // buf can be larger than buf_size
    let mut buf: Vec<i8> = match buf {
        Some(v) => v,
        None => vec![0i8; buf_size],
    };
    // buf will be extended later if buf_size < buf.len()
    //assert!(buf.len() >= buf_size);

    //let mut buf = vec![0i8; buf_num_snv * n];
    //let mut buf = vec![0.0f64; buf_num_snv * n];

    let mut m_in_begin_buf = 0;
    let mut m_begin_buf = 0;
    loop {
        log::debug!("m_in_begin_buf: {}", m_in_begin_buf);
        let m_in_buf = buf_num_snv.min(m_in_file - m_in_begin_buf);
        log::debug!("m_in_buf: {}", m_in_buf);
        let m_in_end_buf = m_in_begin_buf + m_in_buf;
        log::debug!("m_in_end_buf: {}", m_in_end_buf);

        let group_to_m_in_buf =
            super::group_to_m_in_range(group_to_m_in_file, m_in_begin_buf, m_in_end_buf);

        let use_snvs_snv = &use_snvs[m_in_begin_buf..m_in_end_buf];
        let use_snvs_group = super::convert_group_to_use(group_to_m_in_buf.as_ref(), m_in_buf);
        let use_snvs_buf = vec::or_bool_vec(use_snvs_snv, &use_snvs_group);

        // snv only
        let m_buf = vec::count_true(use_snvs_snv);
        log::debug!("m_buf: {}", m_buf);
        let m_end_buf = m_begin_buf + m_buf;

        // loaded snvs to buf
        let m_num_buf = vec::count_true(&use_snvs_buf);
        if m_num_buf != 0 {
            load_genot_snvs_extract_buf(
                fin_genotype.clone(),
                fin_snv.clone(),
                &use_snvs_buf,
                &file_snv_allele_idx_in[m_in_begin_buf..m_in_end_buf],
                n_in,
                use_samples_idx,
                n,
                &mut buf,
            );
            //println!("buf {:?}", &buf[..10]);

            assign_genot_buf(
                &use_snvs_buf,
                use_snvs_snv,
                n,
                &mut buf,
                &mut g_snv_file.as_genot_snvs_mut(m_begin_buf, m_end_buf),
                &mut g_group,
                group_to_m_in_buf.as_ref(),
                fill_missing_group,
            );
        }

        m_begin_buf = m_end_buf;
        m_in_begin_buf = m_in_end_buf;
        assert!(m_in_begin_buf <= m_in_file);
        if m_in_begin_buf == m_in_file {
            break;
        }
    }
    assert_eq!(m_in_begin_buf, m_in_file);
    buf
}

fn assign_genot_buf(
    use_snvs_buf: &[bool],
    use_snvs_snv: &[bool],
    n: usize,
    buf: &mut Vec<i8>,
    g_snv_buf: &mut GenotMut,
    g_group: &mut Option<&mut GenotMut>,
    group_to_m_in_buf: Option<&Vec<Vec<usize>>>,
    // group_to_m_in_buf: Option<&HashMap<usize, Vec<usize>>>,
    fill_missing_group: Option<FillMissingGroup>,
) {
    // for snv m_i -> buf_i
    let (m_to_buf, _) = genotype::create_m_to_buf(&use_snvs_buf, use_snvs_snv);
    assign_genot_extract_buf(g_snv_buf, buf, &m_to_buf);

    if let Some(g_group) = g_group {
        // for group m_in_i -> buf_i
        let (m_in_to_buf, _) = genotype::create_m_in_to_buf(&use_snvs_buf);

        let group_to_m_in = group_to_m_in_buf.unwrap();
        g_group
            .iter_snv_mut()
            .enumerate()
            .par_bridge()
            .for_each_with(
                GenotSnv::new_empty(n),
                |g_snv_tmp, (group_i, mut g_group)| {
                    let group_m_in_is = &group_to_m_in[group_i];
                    // let group_m_in_is = group_to_m_in.get(&mi_loaded);

                    // if let Some(group_m_in_is) = group_m_in_is {
                    group_m_in_is.iter().for_each(|&group_m_in_i| {
                        let buf_i = m_in_to_buf[&group_m_in_i];
                        // let group_mi = m_in_to_m[&group_m_in_i];
                        let buf_mi = &buf[buf_i * n..(buf_i + 1) * n];

                        // initialize g_snv_tmp
                        g_snv_tmp.fill_0();

                        assign_gsnv_from_genot_i8(&mut g_snv_tmp.as_genot_snv_mut_snv(), &buf_mi);
                        match fill_missing_group {
                            Some(FillMissingGroup::Ref) => g_snv_tmp.fill_missing_ref(),
                            None => {}
                        };

                        g_group.or_binary(&g_snv_tmp.as_genot_snv());
                    });
                    // }
                    // otherwise, do nothing
                },
            );
    }
}

// // m_in_i, buf_i, m_i, group_m_i are all different
// fn assign_genot(
//     g_file: &mut GenotMut,
//     fin_genotype: CString,
//     fin_snv: CString,
//     m_in_file: usize,
//     n_in: usize,
//     use_samples_idx: &[i32],
//     n: usize,
//     use_snvs: &[bool],
//     file_snv_allele_idx_in: &[(usize, usize, usize)],
//     buf: Option<Vec<i8>>,
//     mem_buf: Option<usize>,
//     mut g_group: Option<&mut GenotMut>,
//     group_to_m_in: Option<&HashMap<usize, Vec<usize>>>,
//     fill_missing_group: Option<FillMissingGroup>,
//     buf_num_snv: Option<usize>, // for testing
// ) -> Vec<i8> {
//     let buf_num_snv: usize = match buf_num_snv {
//         Some(x) => x,
//         None => {
//             let buf_size_limit = buf_size_limit(mem_buf);
//             // 1 byte (i8) per count in pgenlib
//             let byte_per_snv = n * 1;
//             // f64: 8 byte per count in pgenlib
//             // let byte_per_snv = n * 8;
//             let buf_num_snv_limit: usize = buf_size_limit / byte_per_snv;
//             buf_num_snv_limit.min(m_in_file)
//         }
//     };
//     // 1 byte (i8) per count in pgenlib
//     // let byte_per_snv = n * 1;
//     // f64: 8 byte per count in pgenlib
//     // let byte_per_snv = n * 8;
//     // let buf_num_snv_limit: usize = buf_size_limit / byte_per_snv;
//     // let buf_num_snv: usize = buf_num_snv_limit.min(m_in_file);
//     let buf_size = buf_num_snv * n;

//     // buf can be larger than buf_size
//     let mut buf: Vec<i8> = match buf {
//         Some(v) => v,
//         None => vec![0i8; buf_size],
//     };
//     // buf will be extended later if buf_size < buf.len()
//     //assert!(buf.len() >= buf_size);

//     //let mut buf = vec![0i8; buf_num_snv * n];
//     //let mut buf = vec![0.0f64; buf_num_snv * n];

//     let mut m_in_begin_loaded = 0;
//     let mut m_begin_loaded = 0;
//     loop {
//         log::debug!("m_in_begin_loaded: {}", m_in_begin_loaded);
//         let m_in_read = buf_num_snv.min(m_in_file - m_in_begin_loaded);
//         log::debug!("m_in_read: {}", m_in_read);

//         let m_in_end_loaded = m_in_begin_loaded + m_in_read;
//         let use_snvs_loaded = &use_snvs[m_in_begin_loaded..m_in_end_loaded];

//         let use_snvs_group = convert_group_to_use(group_to_m_in, m_in_read);

//         let use_snvs_buf = vec::or_bool_vec(use_snvs_loaded, &use_snvs_group);

//         // m_loaded_snv
//         let m_loaded = vec::count_true(use_snvs_loaded);
//         log::debug!("m_read: {}", m_loaded);
//         log::debug!("m_in_end_loded: {}", m_in_end_loaded);

//         let file_snv_allele_idx_in_loaded =
//             &file_snv_allele_idx_in[m_in_begin_loaded..m_in_end_loaded];

//         let m_end_loaded = m_begin_loaded + m_loaded;

//         let m_buf_loaded = vec::count_true(&use_snvs_buf);
//         if m_buf_loaded != 0 {
//             let (m_to_buf, _) = genotype::create_m_to_buf(&use_snvs_buf, use_snvs_loaded);

//             load_genot_snvs_extract_buf(
//                 fin_genotype.clone(),
//                 fin_snv.clone(),
//                 //m_in_begin_loaded,
//                 //m_in_end_loaded,
//                 &use_snvs_buf,
//                 // use_snvs_loaded,
//                 file_snv_allele_idx_in_loaded,
//                 n_in,
//                 use_samples_idx,
//                 n,
//                 &mut buf,
//             );
//             //println!("buf {:?}", &buf[..10]);

//             let mut g_file_part = g_file.as_genot_snvs_mut(m_begin_loaded, m_end_loaded);

//             // FIX: use m_loaded to m
//             assign_genot_extract_buf(&mut g_file_part, &buf, &m_to_buf);

//             if let Some(g_group) = &mut g_group {
//                 let group_to_m_in = group_to_m_in.unwrap();

//                 // FIX: wrong
//                 let (m_in_to_buf, _) = genotype::create_m_in_to_buf(&use_snvs_buf);

//                 g_group
//                     .iter_snv_mut()
//                     .enumerate()
//                     .par_bridge()
//                     .for_each_with(
//                         GenotSnv::new_empty(n),
//                         |g_snv_tmp, (mi_loaded, mut g_snv)| {
//                             let group_m_in_is = group_to_m_in.get(&mi_loaded);

//                             if let Some(group_m_in_is) = group_m_in_is {
//                                 group_m_in_is.iter().for_each(|&group_m_in_i| {
//                                     let buf_i = m_in_to_buf[&group_m_in_i];
//                                     // let group_mi = m_in_to_m[&group_m_in_i];
//                                     let buf_mi = &buf[buf_i * n..(buf_i + 1) * n];

//                                     // initialize g_snv_tmp
//                                     g_snv_tmp.fill_0();

//                                     assign_gsnv_from_genot_i8(
//                                         &mut g_snv_tmp.as_genot_snv_mut_snv(),
//                                         &buf_mi,
//                                     );
//                                     match fill_missing_group {
//                                         Some(FillMissingGroup::Ref) => g_snv_tmp.fill_missing_ref(),
//                                         None => {}
//                                     };

//                                     g_snv.or_binary(&g_snv_tmp.as_genot_snv());
//                                 });
//                             }
//                             // otherwise, do nothing
//                         },
//                     );
//             }
//         }

//         m_begin_loaded = m_end_loaded;
//         m_in_begin_loaded = m_in_end_loaded;
//         assert!(m_in_begin_loaded <= m_in_file);
//         if m_in_begin_loaded == m_in_file {
//             break;
//         }
//     }
//     assert_eq!(m_in_begin_loaded, m_in_file);
//     buf
// }

/* fn assign_genot(
    g_chrom: &mut GenotMut,
    fin_genot: CString,
    m_in_chrom: usize,
    n_in: usize,
    use_samples_idx: Vec<i32>,
    //use_samples_idx: Vec<i32>,
    n: usize,
    use_snvs: &[bool],
    // TODO: buf: Option<Vec<i32>> // for chrom
) {
    let buf_size_limit: usize = BUF_SIZE_PGENLIB_MAX;
    // 8 byte per count in pgenlib
    let byte_per_snv = n * 8;
    let buf_num_snv_limit: usize = buf_size_limit / byte_per_snv;
    let buf_num_snv: usize = buf_num_snv_limit.min(m_in_chrom);
    //let buf_size: usize = buf_num_snv * byte_per_snv;
    //assert_eq!(buf_size % byte_per_snv, 0);
    //assert!(buf_size <= buf_size_limit);

    // TMP
    //let buf_num_snv = 10;

    let mut buf = vec![0i8; buf_num_snv * n];
    //let mut buf = vec![0.0f64; buf_num_snv * n];

    let mut m_in_begin_loaded = 0;
    let mut m_begin_loaded = 0;
    loop {
        log::debug!("m_in_begin_loaded: {}", m_in_begin_loaded);
        let m_in_read = buf_num_snv.min(m_in_chrom - m_in_begin_loaded);
        log::debug!("m_in_read: {}", m_in_read);

        let m_in_end_loaded = m_in_begin_loaded + m_in_read;
        let use_snvs_loaded = &use_snvs[m_in_begin_loaded..m_in_end_loaded];
        let (_, m_read) = genotype::create_m_to_m_in(use_snvs_loaded);
        log::debug!("m_read: {}", m_read);
        log::debug!("m_in_end_loded: {}", m_in_end_loaded);

        let m_end_loaded = m_begin_loaded + m_read;

        if m_read != 0 {
            //let buf = load_genot_whole_buf(fin_genot, m_in, n_in, use_samples_idx, n);
            //load_genot_snvs_extract_buf(
            load_genot_snvs_buf(
                fin_genot.clone(),
                m_in_begin_loaded,
                m_in_end_loaded,
                n_in,
                //use_samples_idx.clone(),
                //&mut use_samples_idx,
                &use_samples_idx,
                n,
                &mut buf,
            );
            //println!("buf {:?}", &buf[..10]);

            let mut g_chrom_part = g_chrom.as_genot_snvs_mut(m_begin_loaded, m_end_loaded);
            assign_genot_buf(&mut g_chrom_part, &buf, use_snvs_loaded);
        }

        m_begin_loaded = m_end_loaded;
        m_in_begin_loaded = m_in_end_loaded;
        assert!(m_in_begin_loaded <= m_in_chrom);
        if m_in_begin_loaded == m_in_chrom {
            break;
        }
    }
    assert_eq!(m_in_begin_loaded, m_in_chrom);
} */

fn load_genot_snvs_extract_buf(
    fin_genotype: CString,
    fin_snv: CString,
    use_snvs: &[bool],
    file_snv_allele_idx_in: &[(usize, usize, usize)],
    n_in: usize,
    use_samples_idx: &[i32],
    n: usize,
    buf: &mut Vec<i8>,
    //buf: &mut Vec<f64>,
) {
    let snv_idx: Vec<i32> = file_snv_allele_idx_in
        .iter()
        .map(|x| x.1.try_into().unwrap())
        .collect();
    let allele_idx: Vec<i32> = file_snv_allele_idx_in
        .iter()
        .map(|x| x.2.try_into().unwrap())
        .collect();

    // too large nthr might leads mem error.
    let max_threads_pgenlib = 32;
    let nthr = rayon::current_num_threads().min(max_threads_pgenlib);
    log::debug!("nthr in rust {}", nthr);

    let m_read = vec::count_true(use_snvs);
    buf.resize(m_read * n, 0i8);

    let m_in = use_snvs.len();

    unsafe {
        let _ = pgenlib::pgenreader_load_snvs_extract(
            buf.as_mut_ptr(),
            fin_genotype.as_ptr(),
            fin_snv.as_ptr(),
            //m_in_start.try_into().unwrap(),
            //m_in_end.try_into().unwrap(),
            m_in.try_into().unwrap(),
            //m_read.try_into().unwrap(),
            use_snvs.as_ptr(),
            snv_idx.as_ptr(),
            allele_idx.as_ptr(),
            n_in.try_into().unwrap(),
            use_samples_idx.as_ptr(),
            n.try_into().unwrap(),
            nthr.try_into().unwrap(),
        );
    }
    // println!("buf {:?}", &buf);
    //println!("buf_tmp {:?}", &buf_tmp[..10]);
    //println!("genot_v {:?}", genot_v);
    //bufv
}

// #[allow(dead_code)]
// fn load_genot_snvs_buf(
//     fin_genotype: CString,
//     fin_snv: CString,
//     m_start: usize,
//     m_end: usize,
//     n_in: usize,
//     use_samples_idx: &Vec<i32>,
//     n: usize,
//     buf: &mut Vec<i8>,
// ) {
//     let m_in = m_end - m_start;
//     let use_snvs = vec![true; m_in];

//     load_genot_snvs_extract_buf(
//         fin_genotype,
//         fin_snv,
//         m_start,
//         m_end,
//         &use_snvs,
//         n_in,
//         use_samples_idx,
//         n,
//         //&mut buf,
//         buf,
//     );
//     //genot_v

//     /*
//     let nthr = rayon::current_num_threads();
//     //let nthr = 4;
//     log::debug!("nthr in rust {}", nthr);

//     //log::debug!("buf size {}",(m_end - m_start) * n);
//     buf.resize((m_end - m_start) * n, 0.0f64);

//     // TODO: convert f64->i8 in pgenlib
//     unsafe {
//         let _ = pgenlib::pgenreader_load_snvs(
//             buf.as_mut_ptr(),
//             fin_genot.as_ptr(),
//             m_start.try_into().unwrap(),
//             m_end.try_into().unwrap(),
//             n_in.try_into().unwrap(),
//             use_samples_idx.as_ptr(),
//             //use_samples_idx.as_mut_ptr(),
//             n.try_into().unwrap(),
//             nthr.try_into().unwrap(),
//         );
//     }
//     //println!("genot_v {:?}", genot_v);
//     //bufv */
// }

// /// won't work for 1M SNVs x 300K samples for f64
// /// not tried for i8
// #[allow(dead_code)]
// fn load_genot_whole_buf(
//     fin_genotype: CString,
//     fin_snv: CString,
//     m_in: usize,
//     n_in: usize,
//     use_samples_idx: &[i32],
//     n: usize,
// ) -> Vec<i8> {
//     let m_start = 0;
//     let m_end = m_in;
//     let use_snvs = vec![true; m_in];

//     let mut genot_v = vec![0i8; m_in * n];

//     load_genot_snvs_extract_buf(
//         fin_genotype,
//         fin_snv,
//         m_start,
//         m_end,
//         &use_snvs,
//         n_in,
//         use_samples_idx,
//         n,
//         &mut genot_v,
//     );
//     genot_v

//     /*     //let mut genot_v = vec![0.0f64; m_in * n];
//     unsafe {
//         let _ = pgenlib::pgenreader_load_whole(
//             genot_v.as_mut_ptr(),
//             fin_genot.as_ptr(),
//             m_in.try_into().unwrap(),
//             n_in.try_into().unwrap(),
//             use_samples_idx.as_mut_ptr(),
//             n.try_into().unwrap(),
//             nthr.try_into().unwrap(),
//         );
//     }
//     //println!("genot_v {:?}", genot_v);
//     genot_v */
// }

fn assign_genot_extract_buf(g: &mut GenotMut, buf: &[i8], m_to_buf: &HashMap<usize, usize>) {
    let n = g.n();

    g.iter_snv_mut()
        .enumerate()
        .par_bridge()
        .for_each(|(mi, mut g_snv)| {
            let buf_i = m_to_buf[&mi];
            let buf_mi = &buf[buf_i * n..(buf_i + 1) * n];
            // let buf_mi = &buf[mi * n..(mi + 1) * n];

            //println!("buf_mi {:?}", &buf_mi[..10]);
            assign_gsnv_from_genot_i8(&mut g_snv, &buf_mi);
        });
}

// fn assign_genot_extract_buf(g: &mut GenotMut, buf: &[i8]) {
//     //assert
//     //assert_eq!(g.m(), m_read);
//     let n = g.n();

//     g.iter_snv_mut()
//         .enumerate()
//         .par_bridge()
//         .for_each(|(mi, mut g_snv)| {
//             let buf_mi = &buf[mi * n..(mi + 1) * n];

//             //println!("buf_mi {:?}", &buf_mi[..10]);
//             assign_gsng_from_genot_i8(&mut g_snv, &buf_mi);
//         });
// }

//fn assign_genot_buf(g: &mut GenotMut, buf: &[f64], use_snvs: &[bool]) {
#[allow(dead_code)]
fn assign_gsnvs_buf(g: &mut GenotMut, buf: &[i8], use_snvs: &[bool]) {
    let (m_to_m_in, m_read) = genotype::create_m_to_m_in(use_snvs);

    //assert
    assert_eq!(g.m(), m_read);
    let n = g.n();

    g.iter_snv_mut()
        .enumerate()
        .par_bridge()
        .for_each(|(mi, mut g_snv)| {
            let m_in_i = m_to_m_in[&mi];
            let buf_mi = &buf[m_in_i * n..(m_in_i + 1) * n];

            //println!("buf_mi {:?}", &buf_mi[..10]);
            assign_gsnv_from_genot_i8(&mut g_snv, &buf_mi);
        });
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{sample, snv};
    use std::path::PathBuf;

    fn setup_test3() -> (
        GenotFiles,
        Vec<bool>,
        usize,
        usize,
        Vec<bool>,
        Vec<(usize, usize, usize)>,
        Vec<bool>,
    ) {
        let fin = PathBuf::from("../../test/data/toy3/genot");
        let fins_genot = GenotFiles::new_plink2vzs(vec![fin.clone()]);

        let m_in: usize = genot_io::compute_num_snv(&fins_genot);
        log::debug!("{}", m_in);
        let n_in: usize = genot_io::compute_num_sample(&fins_genot);
        log::debug!("{}", n_in);
        // load snvs
        let (snvs_in, file_snv_allele_idx_in) = genot_io::load_snvs(&fins_genot);
        let (use_snvs, m) = snv::make_use_snvs_buf(None, &snvs_in);
        let (use_samples, n) = sample::make_use_samples_buf(None, &fins_genot);
        let ys = vec![
            true, true, true, true, true, false, false, false, false, false,
        ];

        (
            fins_genot,
            ys,
            m,
            n,
            use_snvs,
            file_snv_allele_idx_in,
            use_samples,
        )
    }

    fn setup_test3_part() -> (
        GenotFiles,
        Vec<bool>,
        usize,
        usize,
        Vec<bool>,
        Vec<(usize, usize, usize)>,
        Vec<bool>,
    ) {
        let fin = PathBuf::from("../../test/data/toy3/genot");
        let fins_genot = GenotFiles::new_plink2vzs(vec![fin.clone()]);

        let m_in: usize = genot_io::compute_num_snv(&fins_genot);
        log::debug!("{}", m_in);
        let n_in: usize = genot_io::compute_num_sample(&fins_genot);
        log::debug!("{}", n_in);
        let m = 2;
        let use_snvs = vec![true, false, true];
        let n = 5;
        let use_samples = vec![
            false, true, false, true, false, true, false, true, false, true,
        ];
        let file_snv_allele_idx_in = vec![(0, 0, 0), (0, 1, 0), (0, 2, 0)];
        let ys = vec![
            true, true, true, true, true, false, false, false, false, false,
        ];

        (
            fins_genot,
            ys,
            m,
            n,
            use_snvs,
            file_snv_allele_idx_in,
            use_samples,
        )
    }

    fn setup_test3ref() -> (
        GenotFiles,
        Vec<bool>,
        usize,
        usize,
        Vec<bool>,
        Vec<(usize, usize, usize)>,
        Vec<bool>,
    ) {
        let fin = PathBuf::from("../../test/data/toy3/genot.ref");
        let fins_genot = GenotFiles::new_plink2vzs(vec![fin.clone()]);

        let m_in: usize = genot_io::compute_num_snv(&fins_genot);
        log::debug!("{}", m_in);
        let n_in: usize = genot_io::compute_num_sample(&fins_genot);
        log::debug!("{}", n_in);
        let (snvs_in, file_snv_allele_idx_in) = genot_io::load_snvs(&fins_genot);
        let (use_snvs, m) = snv::make_use_snvs_buf(None, &snvs_in);
        let (use_samples, n) = sample::make_use_samples_buf(None, &fins_genot);
        let ys = vec![
            true, true, true, true, true, false, false, false, false, false,
        ];

        (
            fins_genot,
            ys,
            m,
            n,
            use_snvs,
            file_snv_allele_idx_in,
            use_samples,
        )
    }

    fn setup_test3_plink2() -> (
        GenotFiles,
        Vec<bool>,
        usize,
        usize,
        Vec<bool>,
        Vec<(usize, usize, usize)>,
        Vec<bool>,
    ) {
        let fin = PathBuf::from("../../test/data/toy3/genot");
        let fins_genot = GenotFiles::new_plink2vzs(vec![fin.clone()]);

        let m_in: usize = genot_io::compute_num_snv(&fins_genot);
        log::debug!("{}", m_in);
        let n_in: usize = genot_io::compute_num_sample(&fins_genot);
        log::debug!("{}", n_in);
        // load snvs
        let (snvs_in, file_snv_allele_idx_in) = genot_io::load_snvs(&fins_genot);
        let (use_snvs, m) = snv::make_use_snvs_buf(None, &snvs_in);
        let (use_samples, n) = sample::make_use_samples_buf(None, &fins_genot);
        let ys = vec![
            true, true, true, true, true, false, false, false, false, false,
        ];

        (
            fins_genot,
            ys,
            m,
            n,
            use_snvs,
            file_snv_allele_idx_in,
            use_samples,
        )
    }

    // multi-files
    fn setup_test7_plink2() -> (
        GenotFiles,
        Vec<bool>,
        usize,
        usize,
        Vec<bool>,
        Vec<(usize, usize, usize)>,
        Vec<bool>,
    ) {
        let fins = vec![
            PathBuf::from("../../test/data/toy7/genot_1.ref"),
            PathBuf::from("../../test/data/toy7/genot_2.ref"),
        ];
        let fins_genot = GenotFiles::new_plink2vzs(fins.clone());

        let m_in: usize = genot_io::compute_num_snv(&fins_genot);
        log::debug!("{}", m_in);
        let n_in: usize = genot_io::compute_num_sample(&fins_genot);
        log::debug!("{}", n_in);
        // load snvs
        let (_, file_snv_allele_idx_in) = genot_io::load_snvs(&fins_genot);
        let use_snvs = vec![true, true, false, true, false, true, true, true, true, true];
        let m = vec::count_true(&use_snvs);
        let (use_samples, n) = sample::make_use_samples_buf(None, &fins_genot);
        let ys = vec![true, true, true, true, true, true, true, true, true, true];

        (
            fins_genot,
            ys,
            m,
            n,
            use_snvs,
            file_snv_allele_idx_in,
            use_samples,
        )
    }

    // multi-allelic
    fn setup_test8_plink2() -> (
        GenotFiles,
        Vec<bool>,
        usize,
        usize,
        Vec<bool>,
        Vec<(usize, usize, usize)>,
        Vec<bool>,
    ) {
        let fin = PathBuf::from("../../test/data/toy8/genot.ref");
        let fins_genot = GenotFiles::new_plink2vzs(vec![fin.clone()]);

        let m_in: usize = genot_io::compute_num_snv(&fins_genot);
        log::debug!("{}", m_in);
        let n_in: usize = genot_io::compute_num_sample(&fins_genot);
        log::debug!("{}", n_in);
        // load snvs
        let (snvs_in, file_snv_allele_idx_in) = genot_io::load_snvs(&fins_genot);
        let (use_snvs, m) = snv::make_use_snvs_buf(None, &snvs_in);
        let (use_samples, n) = sample::make_use_samples_buf(None, &fins_genot);
        let ys = vec![true, true, true, true, true, true, true, true, true, true];

        (
            fins_genot,
            ys,
            m,
            n,
            use_snvs,
            file_snv_allele_idx_in,
            use_samples,
        )
    }

    #[test]
    fn test_assign_pred_from_bed() {
        let mut g = GenotSnv::new_empty(6);
        // [2, 0, 3, 0, 1, 0]
        //let pbuf = vec![2.0f64, 0.0, 3.0, 0.0, 1.0, 0.0];
        let pbuf = vec![2i8, 0, 3, 0, 1, 0];

        assign_gsnv_from_genot_i8(&mut g.as_genot_snv_mut_snv(), &pbuf);
        assert_eq!(g.vals(), vec![2u8, 0, 3, 0, 1, 0]);
    }

    #[test]
    fn test_generate_genot_snv_plink2() {
        let (fin_genot, _, _, n, _, _, use_samples) = setup_test3_plink2();
        //let use_snvs = vec![true; use_snvs.len()];
        //let use_samples = vec![true; use_samples.len()];
        //let g = generate_genot_plink(&fin, gfmt, m, n, &use_snvs, Some(&use_samples), true);
        let m_in_i = 2;
        let g = generate_genot_snv_file_plink2(
            &fin_genot.files()[0],
            m_in_i,
            n,
            Some(&use_samples),
            false,
        );
        assert_eq!(g.vals(), vec![2, 0, 1, 0, 1, 2, 0, 1, 0, 3]);
    }

    #[test]
    fn test_generate_genot_snv_plink2vzs() {
        let (fin_genot, _, _, n, _, _, use_samples) = setup_test3();
        //let use_snvs = vec![true; use_snvs.len()];
        //let use_samples = vec![true; use_samples.len()];
        //let g = generate_genot_plink(&fin, gfmt, m, n, &use_snvs, Some(&use_samples), true);
        let m_in_i = 2;
        let g = generate_genot_snv_file_plink2(
            &fin_genot.files()[0],
            m_in_i,
            n,
            Some(&use_samples),
            false,
        );
        assert_eq!(g.vals(), vec![2, 0, 1, 0, 1, 2, 0, 1, 0, 3]);
    }

    #[test]
    fn test_generate_genot_plink2vzs_ref() {
        let (fin_genot, _, _, n, _, _, use_samples) = setup_test3ref();
        let m_in_i = 2;
        let g = generate_genot_snv_file_plink2(
            &fin_genot.files()[0],
            m_in_i,
            n,
            Some(&use_samples),
            false,
        );
        assert_eq!(g.vals(), vec![0, 2, 1, 2, 1, 0, 2, 1, 2, 3]);
    }

    #[test]
    fn test_generate_genot_snv_plink2vzs_part() {
        let (fin_genot, _, _, n, _, _, use_samples) = setup_test3_part();
        let m_in_i = 2;
        let g = generate_genot_snv_file_plink2(
            &fin_genot.files()[0],
            m_in_i,
            n,
            Some(&use_samples),
            false,
        );
        assert_eq!(g.vals(), vec![0, 0, 2, 1, 3]);
    }

    #[test]
    fn test_generate_genot_plink2vzs_part() {
        let (fin_genot, _, m, n, use_snvs, file_snv_allele_idx_in, use_samples) =
            setup_test3_part();
        let g = generate_genot_plink2(
            &fin_genot,
            m,
            n,
            Some(&use_snvs),
            &file_snv_allele_idx_in,
            None,
            Some(&use_samples),
            None,
            None,
            None,
            Some(1),
        );
        let mut giter = g.iter_snv();
        // [2,1,0,0,0,2,1,0,0,3]
        assert_eq!(giter.next().unwrap().vals(), vec![1, 0, 2, 0, 3]);
        assert_eq!(giter.next().unwrap().vals(), vec![0, 0, 2, 1, 3]);
        assert_eq!(giter.next(), None);
    }

    #[test]
    fn test_generate_genot_plink2vzs_aa() {
        let (fin_genot, _, m, n, use_snvs, file_snv_allele_idx_in, use_samples) = setup_test3();
        let g = generate_genot_plink2(
            &fin_genot,
            m,
            n,
            Some(&use_snvs),
            &file_snv_allele_idx_in,
            None,
            Some(&use_samples),
            None,
            None,
            None,
            Some(1),
        );
        let mut giter = g.iter_snv();
        // [2,1,0,0,0,2,1,0,0,3]
        assert_eq!(
            giter.next().unwrap().vals(),
            vec![2, 1, 0, 0, 0, 2, 1, 0, 0, 3]
        );
        assert_eq!(
            giter.next().unwrap().vals(),
            vec![2, 0, 1, 1, 0, 1, 0, 2, 0, 3]
        );
        assert_eq!(
            giter.next().unwrap().vals(),
            vec![2, 0, 1, 0, 1, 2, 0, 1, 0, 3]
        );
        assert_eq!(giter.next(), None);
    }

    #[test]
    fn test_generate_genot_plink2vzs_toy7() {
        // TODO
        let (fin_genot, _, m, n, use_snvs, file_snv_allele_idx_in, use_samples) =
            setup_test7_plink2();
        let g = generate_genot_plink2(
            &fin_genot,
            m,
            n,
            Some(&use_snvs),
            &file_snv_allele_idx_in,
            None,
            Some(&use_samples),
            None,
            None,
            None,
            Some(3),
        );
        let mut giter = g.iter_snv();
        // chr1:100:A>C
        assert_eq!(giter.next().unwrap().vals(), vec![0, 0, 0, 0, 0, 0, 0, 1]);
        // chr1:80:T>A
        assert_eq!(giter.next().unwrap().vals(), vec![0, 1, 0, 0, 0, 1, 2, 3]);
        // chr1:90:A>C
        assert_eq!(giter.next().unwrap().vals(), vec![3, 0, 1, 2, 0, 0, 1, 3]);
        // chr2:150:A>C
        assert_eq!(giter.next().unwrap().vals(), vec![0, 1, 0, 0, 2, 0, 0, 1]);
        // chr2:150:A>G
        assert_eq!(giter.next().unwrap().vals(), vec![0, 0, 1, 0, 0, 2, 0, 1]);
        // chr2:150:A>T
        assert_eq!(giter.next().unwrap().vals(), vec![0, 0, 0, 1, 0, 0, 2, 0]);
        // chr2:10:ACC>A
        assert_eq!(giter.next().unwrap().vals(), vec![3, 0, 0, 1, 1, 0, 2, 3]);
        // chr2:10:ACC>AC
        assert_eq!(giter.next().unwrap().vals(), vec![3, 1, 0, 0, 0, 2, 0, 3]);
        assert_eq!(giter.next(), None);
    }

    #[test]
    fn test_load_genot_snvs_extract_buf_toy8() {
        // TODO
        let (fin_genot, _, m, n, use_snvs, file_snv_allele_idx_in, use_samples) =
            setup_test8_plink2();

        let g = generate_genot_plink2(
            &fin_genot,
            m,
            n,
            Some(&use_snvs),
            &file_snv_allele_idx_in,
            None,
            Some(&use_samples),
            None,
            None,
            None,
            Some(3),
        );
        let mut giter = g.iter_snv();
        // chr1:1:A>C
        assert_eq!(
            giter.next().unwrap().vals(),
            vec![3, 0, 1, 2, 0, 0, 0, 0, 0, 0]
        );
        // chr1:2:A>AC
        assert_eq!(
            giter.next().unwrap().vals(),
            vec![3, 0, 1, 2, 0, 0, 1, 0, 0, 0]
        );
        // chr1:2:A>ACC
        assert_eq!(
            giter.next().unwrap().vals(),
            vec![3, 0, 0, 0, 1, 2, 1, 0, 0, 0]
        );
        // chr1:3:A>ACC
        assert_eq!(
            giter.next().unwrap().vals(),
            vec![3, 0, 0, 1, 0, 0, 2, 0, 1, 1]
        );
        // chr1:3:A>AC
        assert_eq!(
            giter.next().unwrap().vals(),
            vec![3, 0, 1, 0, 0, 2, 0, 0, 1, 0]
        );
        // chr1:3:A>ACCC
        assert_eq!(
            giter.next().unwrap().vals(),
            vec![3, 0, 0, 0, 1, 0, 0, 2, 0, 1]
        );
        // chr1:4:A>ACC
        assert_eq!(
            giter.next().unwrap().vals(),
            vec![3, 0, 0, 1, 0, 0, 2, 0, 1, 1]
        );
        // chr1:4:A>AC
        assert_eq!(
            giter.next().unwrap().vals(),
            vec![3, 0, 1, 0, 0, 2, 0, 0, 1, 0]
        );
        // chr1:4:A>ACCC
        assert_eq!(
            giter.next().unwrap().vals(),
            vec![3, 0, 0, 0, 1, 0, 0, 2, 0, 1]
        );
        // chr1:5:A>ACCCC
        assert_eq!(
            giter.next().unwrap().vals(),
            vec![3, 0, 0, 0, 0, 1, 1, 2, 0, 0]
        );
        // chr1:5:A>AC
        assert_eq!(
            giter.next().unwrap().vals(),
            vec![3, 0, 1, 0, 0, 0, 1, 0, 0, 0]
        );
        // chr1:5:A>ACC
        assert_eq!(
            giter.next().unwrap().vals(),
            vec![3, 0, 0, 1, 0, 0, 0, 0, 1, 0]
        );
        // chr1:5:A>ACCC
        assert_eq!(
            giter.next().unwrap().vals(),
            vec![3, 0, 0, 0, 1, 0, 0, 0, 1, 0]
        );
        // chr1:6:AC>A
        assert_eq!(
            giter.next().unwrap().vals(),
            vec![3, 2, 1, 0, 0, 0, 0, 0, 0, 0]
        );
        // chr1:6:AC>ACC
        assert_eq!(
            giter.next().unwrap().vals(),
            vec![3, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        );
        assert_eq!(giter.next(), None);
    }
}
