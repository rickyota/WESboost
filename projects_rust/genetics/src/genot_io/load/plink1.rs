use crate::genot::prelude::*;
use crate::{alloc, FillMissing, FillMissingGroup, GenotFile, GenotFiles};
use crate::{genot_io, vec};
use crate::{genotype, textfile};

use rayon::prelude::*;
use std::collections::HashMap;
use std::error::Error;
use std::fs::File;
use std::io::{prelude::*, BufRead, BufReader, Read, SeekFrom};
use std::path::Path;
use std::time::Instant;
//use std::io::{prelude::*, BufReader, SeekFrom};
//use std::io::{prelude::*, Read, SeekFrom}; // for seek

/*
// save
// load single snv
pub fn generate_genot_snv_plink(
    fin: &Path,
    gfmt: GenotFormat,
    mi: usize,
    n: usize,
    use_samples: Option<&[bool]>,
    fill_missing: bool,
) -> GenotSnv {
    let reader = BufReader::new(File::open(genot_io::fname_plinks_genot(fin, gfmt, None)).unwrap());

    let mut g_snv = GenotSnv::new_empty(n);

    load_snv_read(
        &mut g_snv.as_genot_snv_mut_snv(),
        reader,
        mi,
        use_samples,
        n,
    );

    //if !use_missing {
    if fill_missing {
        //super::fill_missing_snv(&mut g_snv.as_genot_snv_mut_snv());
        g_snv.as_genot_snv_mut_snv().fill_missing_mode()
    }

    g_snv
}

fn load_snv_read<R: BufRead + Seek>(
    g_snv: &mut GenotSnvMut,
    mut reader: R,
    mi: usize,
    use_samples: Option<&[bool]>,
    n: usize,
) {
    let n_in = if let Some(v) = use_samples {
        v.len()
    } else {
        n
    };
    assert!(n <= n_in);

    // load to buf
    let byte_per_snv = genot_io::bed_per_snv_size(n_in);
    let buf_size = byte_per_snv;
    let mut buf: Vec<B8_2> = vec![0; buf_size];
    let buf_begin = 3 + mi * byte_per_snv;
    reader.seek(SeekFrom::Start(buf_begin as u64)).unwrap();
    let loaded_byte = reader.read(&mut buf).unwrap();
    assert_eq!(loaded_byte, byte_per_snv);

    assign_pred_from_bed(g_snv, &buf, use_samples);
}
*/

const BUF_SIZE_BED_MAX: usize = 64 * 1024 * 1024 * 1024;
const BUF_SIZE_BED_MIN: usize = 1 * 1024 * 1024 * 1024;

// Maximum file size of GenotFiles
//fn load_buf_size_max(fin: &Path, gfmt: GenotFormat) -> usize {
fn load_file_size_max(fins_genot: &GenotFiles) -> usize {
    let mut buf_size_max = 0usize;
    for fin_genot in fins_genot.files().iter() {
        let file_size = textfile::file_size(&fin_genot.genotype_file());
        if let Some(x) = file_size {
            buf_size_max = buf_size_max.max(x);
        }
    }
    buf_size_max
}

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
                alloc::mem_gb(BUF_SIZE_BED_MIN),
                alloc::mem_gb(x),
            );
            if genot_byte + BUF_SIZE_BED_MIN > x {
                panic!("Memory insufficient on preparing Genot.")
            }
            x - genot_byte
        }
        None => {
            log::debug!(
                "Could not get available memory; assume there is {} GB available memory.",
                alloc::mem_gb(BUF_SIZE_BED_MAX)
            );
            BUF_SIZE_BED_MAX - genot_byte
        }
    };

    mem_buf
}

fn buf_read(fin_genot: &GenotFile) -> BufReader<File> {
    // using  File.read() did not speed up
    //let fin_chrom = plink::fname_chrom(fin, Some(chrom_i));
    // 512 KB
    //let reader_cap = 512usize * 1024;
    // 1 MB
    let reader_cap = 1usize * 1024 * 1024;
    let reader =
        BufReader::with_capacity(reader_cap, File::open(fin_genot.genotype_file()).unwrap());
    reader
}

fn buf_size_limit(mem: Option<usize>) -> usize {
    let buf_size_limit = match mem {
        Some(x) => x.min(BUF_SIZE_BED_MAX),
        None => {
            log::debug!(
                "Could not get available memory; assume there is {:.3} GB available memory.",
                alloc::mem_gb(BUF_SIZE_BED_MAX)
            );
            BUF_SIZE_BED_MAX
        }
    };
    log::debug!("buf_size_limit: {} GB", alloc::mem_gb(buf_size_limit));
    buf_size_limit
}

fn buf_size_max(mem_buf: Option<usize>, fins_genot: &GenotFiles) -> usize {
    let buf_size_limit = buf_size_limit(mem_buf);

    log::debug!("buf_size_limit: {:.3} GB", alloc::mem_gb(buf_size_limit));

    let buf_size_max = load_file_size_max(fins_genot).min(buf_size_limit);

    buf_size_max
}

// pub fn group_to_m_in_file(
//     group_to_m_in: Option<&HashMap<usize, Vec<usize>>>,
//     m_in_begin: usize,
//     m_in_end: usize,
// ) -> Option<HashMap<usize, Vec<usize>>> {
//     if group_to_m_in.is_none() {
//         return None;
//     }

//     let mut group_to_m_in_chrom = HashMap::new();
//     let group_to_m_in = group_to_m_in.unwrap();

//     for (group_i, m_in_set) in group_to_m_in.iter() {
//         let m_in_set_chrom: Vec<usize> = m_in_set
//             .iter()
//             .filter(|&mi| m_in_begin <= *mi && *mi < m_in_end)
//             .map(|&mi| mi - m_in_begin)
//             .collect();
//         if !m_in_set_chrom.is_empty() {
//             group_to_m_in_chrom.insert(*group_i, m_in_set_chrom);
//         }
//     }

//     Some(group_to_m_in_chrom)
// }

/// Sequentially load part of bed and convert into predictions since bed could be too large to load on mem.
/// predictions: m*2*n
// load whole is fastest
pub fn generate_genot_plink1(
    fins_genot: &GenotFiles,
    m_snv: usize,
    // m: usize, // = m_snv + m_group
    n: usize,
    use_snvs: Option<&[bool]>,
    group_to_m_in: Option<Vec<Vec<usize>>>,
    // group_to_m_in: Option<HashMap<usize, Vec<usize>>>,
    use_samples: Option<&[bool]>,
    fill_missing: Option<FillMissing>,
    fill_missing_group: Option<FillMissingGroup>,
    //use_missing: bool,
    //fill_missing_mode: bool,
    mem: Option<usize>,
    //m_set: Option<usize>,
) -> Genot {
    let start = Instant::now();

    let m_group = if group_to_m_in.is_some() {
        group_to_m_in.as_ref().unwrap().len()
    } else {
        0
    };

    let m = m_snv + m_group;

    log::debug!("to prepare Genot m, n: {}, {}", m, n);
    let mem_buf_limit = mem_buf_limit(m, n, mem);

    // TODO: better way
    let use_snvs_v = vec![true; m_snv];
    let use_snvs = match use_snvs {
        Some(x) => x,
        None => &use_snvs_v,
    };

    let mut g = Genot::new_zeros(m, n);

    let (mut g_snvs, mut g_group) = if group_to_m_in.is_some() {
        let (g_snvs, g_group) = g.split_genot_mut(m_snv);
        assert_eq!(g_snvs.m(), m_snv);
        assert_eq!(g_group.m(), m_group);
        (g_snvs, Some(g_group))
    } else {
        (g.as_genot_mut(), None)
    };

    // Two patterns to assign predictions
    // bed file is split into chrom or not
    //let is_split_chrom = fins_genot.judge_split_chrom();

    //if is_split_chrom {
    load_genot_files(
        //mem,
        fins_genot,
        &mut g_snvs,
        &mut g_group,
        //g_snvs,
        use_snvs,
        group_to_m_in,
        use_samples,
        fill_missing_group,
        Some(mem_buf_limit),
    );
    //} else {
    // If input is one file
    //let reader = buf_read(fins_genot, None);
    //assign_genot_toggle(
    //    &mut g_snvs.as_genot_mut(),
    //    reader,
    //    use_snvs,
    //    use_samples,
    //    None,
    //    Some(mem_buf),
    //    g_group.as_mut(),
    //    group_to_m_in.as_ref(),
    //);
    //}

    // missing
    super::fill_missing_g_snvs(&mut g_snvs, fill_missing);

    let end = start.elapsed();
    log::info!("It took {} seconds to generate genot.", end.as_secs());

    g
}

// TODO: clean
/// If input is one file
///let reader = buf_read(fins_genot, None);
///assign_genot_toggle(
///    &mut g_snvs.as_genot_mut(),
///    reader,
///    use_snvs,
///    use_samples,
///    None,
///    Some(mem_buf),
///    g_group.as_mut(),
///    group_to_m_in.as_ref(),
///);
fn load_genot_files(
    fins_genot: &GenotFiles,
    g_snvs: &mut GenotMut,
    g_group: &mut Option<GenotMut>,
    use_snvs: &[bool],
    group_to_m_in: Option<Vec<Vec<usize>>>,
    // group_to_m_in: Option<HashMap<usize, Vec<usize>>>,
    use_samples: Option<&[bool]>,
    fill_missing_group: Option<FillMissingGroup>,
    mem_buf_limit: Option<usize>,
    //mem: Option<usize>,
) {
    let mut m_begin = 0;
    let mut m_in_begin = 0;

    // min of mem_buf, MAX_BED, and file size
    let buf_size_max = buf_size_max(mem_buf_limit, fins_genot);
    log::debug!("buf_size_max {}", buf_size_max);
    let mut buf: Vec<B8_2> = vec![0; buf_size_max];

    for (_file_i, fin_genot) in fins_genot.files().iter().enumerate() {
        log::debug!("Loading file {:?}", fin_genot);
        let m_in_file = genot_io::compute_num_snv_file_bi_allelic(fin_genot);
        let m_in_end = m_in_begin + m_in_file;
        log::debug!("m_in_file {}", m_in_file);
        log::debug!("m_in_end {}", m_in_end);
        let m_file = vec::count_true(&use_snvs[m_in_begin..m_in_end]);
        let m_end = m_begin + m_file;

        let group_to_m_in_file =
            super::group_to_m_in_range(group_to_m_in.as_ref(), m_in_begin, m_in_end);

        // None -> true
        // Some(x) -> x.is_empty()
        let is_group_empty = group_to_m_in_file.as_ref().map_or(true, |x| x.is_empty());

        // skip if no snvs to load
        if m_file == 0 && is_group_empty {
            m_in_begin = m_in_end;
            continue;
        }

        // TODO
        //plink::check_valid_bed(fin_chrom, None, m_in_chrom, n_in).unwrap();

        let reader = buf_read(fin_genot);

        // reuse buf
        buf = load_genot_file_toggle(
            &mut g_snvs.as_genot_snvs_mut(m_begin, m_end),
            reader,
            &use_snvs[m_in_begin..m_in_end],
            use_samples,
            Some(buf),
            Some(buf_size_max),
            g_group.as_mut(),
            group_to_m_in_file.as_ref(),
            fill_missing_group,
        );

        m_begin = m_end;
        m_in_begin = m_in_end;
    }
    assert_eq!(m_in_begin, use_snvs.len(), "Sth wrong.");
}

fn load_genot_file_toggle<R: BufRead + Seek>(
    g_snv_file: &mut GenotMut,
    reader: R,
    //mut reader: R,
    use_snvs: &[bool],
    use_samples: Option<&[bool]>,
    buf: Option<Vec<B8_2>>,
    mem: Option<usize>,
    // mut is necessary
    // https://stackoverflow.com/questions/74763962/when-is-mut-required-when-passing-an-optionmut-t
    //mut g_group: Option<&mut GenotMut>,
    g_group: Option<&mut GenotMut>,
    group_to_m_in_file: Option<&Vec<Vec<usize>>>,
    // group_to_m_in_file: Option<&HashMap<usize, Vec<usize>>>,
    fill_missing_group: Option<FillMissingGroup>,
) -> Vec<B8_2> {
    //let _m = vec::count_true(&use_snvs);
    //if (m as f64)/(use_snvs.len() as f64)<0.1{

    // error on ukbrap
    //    log::debug!("use load part()");
    //    assign_genot_loadpart(
    //        g_chrom,
    //        reader,
    //        use_snvs,
    //        use_samples,
    //        buf,
    //        mem,
    //        g_group,
    //        group_to_m_in,
    //        )
    //}else{
    //    log::debug!("use load whole()");
    load_genot_file(
        g_snv_file,
        reader,
        use_snvs,
        use_samples,
        buf,
        mem,
        g_group,
        group_to_m_in_file,
        fill_missing_group,
        None,
    )
    //}
}

/*
// directly use assign_predictions() applicable to both
fn assign_predictions_toggle<R: BufRead + Seek>(
    g_chrom: &mut GenotMut,
    mut reader: R,
    use_snvs: &[bool],
    use_samples: Option<&[bool]>,
) {
    let file_size = reader.seek(SeekFrom::End(0)).unwrap() as usize;
    log::debug!("seek {}", file_size);

    log::debug!("Use whole strategy to load.");
    assign_predictions_whole(g_chrom, reader, use_snvs, use_samples);
    //log::debug!("Use chunk strategy to load.");
    //assign_predictions_chunk(g_chrom, reader, use_snvs, use_samples);

    //// 64 GB
    //let thres_size = 64 * 1024 * 1024 * 1024;
    //// whole have to allocate memory same size as file
    //if file_size > thres_size {
    //    log::debug!("Use chunk strategy to load.");
    //    assign_predictions_chunk(g_chrom, reader, use_snvs, use_samples);
    //} else {
    //    log::debug!("Use whole strategy to load.");
    //    assign_predictions_whole(g_chrom, reader, use_snvs, use_samples);
    //}
} */

// TODO: use SIMD
// can I implement buf: Option<&mut Vec<u8>> ?
// now unknown way to declare Vec<u8> with longer lifetime when buf=None
// What happens if no var is loaded in each loop? can I shorten time?
/// available for both whole and chunk
fn load_genot_file<R: BufRead + Seek>(
    g_snvs_file: &mut GenotMut,
    mut reader: R,
    use_snvs: &[bool],
    use_samples: Option<&[bool]>,
    // cannot use &mut Vec<B8_2> since it cannot resize??
    buf: Option<Vec<B8_2>>,
    buf_size_max: Option<usize>,
    // mut is necessary
    // https://stackoverflow.com/questions/74763962/when-is-mut-required-when-passing-an-optionmut-t
    mut g_group: Option<&mut GenotMut>,
    group_to_m_in_file: Option<&Vec<Vec<usize>>>,
    // group_to_m_in_file: Option<&HashMap<usize, Vec<usize>>>,
    fill_missing_group: Option<FillMissingGroup>,
    buf_num_snv: Option<usize>, // for testing
) -> Vec<B8_2> {
    let m_in_chrom = use_snvs.len();

    let n_in = match use_samples {
        Some(v) => v.len(),
        None => g_snvs_file.n(),
    };

    let n = match use_samples {
        Some(v) => vec::count_true(v),
        None => n_in,
    };

    // FIXME: available mem?
    // check if 32 GB? 16GB? remains or not
    // reading smaller than buf_size is not error but reason is unknown
    // [here](https://doc.rust-lang.org/std/io/trait.Read.html#tymethod.read)
    // -> to deal with this, every time use Seek to change position

    let byte_per_snv = bed_per_snv_size(n_in);
    let buf_num_snv = match buf_num_snv {
        Some(x) => x,
        None => {
            // already done in load_genot_files, but do again for mem_buf=None
            let buf_size_max = buf_size_limit(buf_size_max);
            // let byte_per_snv = bed_per_snv_size(n_in);
            let buf_num_snv_limit: usize = buf_size_max / byte_per_snv;
            buf_num_snv_limit.min(m_in_chrom)
        }
    };

    // already done in load_genot_files, but do again for mem_buf=None
    // let buf_size_max = buf_size_limit(buf_size_max);
    // let byte_per_snv = bed_per_snv_size(n_in);
    // let buf_num_snv_limit: usize = buf_size_max / byte_per_snv;
    // let buf_num_snv: usize = buf_num_snv_limit.min(m_in_chrom);
    let buf_size: usize = buf_num_snv * byte_per_snv;
    assert_eq!(buf_size % byte_per_snv, 0);
    // assert!(buf_size <= buf_size_max);

    //let mut reader = BufReader::new(File::open(plink::fname_bed(fin_chrom, None)).unwrap());
    // TODO: check buf_size < available

    // first, create large buf here and resize later in loop
    let mut buf: Vec<B8_2> = match buf {
        Some(v) => v,
        None => vec![0; buf_size],
    };
    // even if buf size < buf_size, resize later

    let mut m_in_begin_buf = 0;
    let mut m_begin_buf = 0;
    // read buf length in one loop
    loop {
        //log::debug!(
        //    "m_in_begin, m_begin {},{}",
        //    m_in_begin_loaded, m_begin_loaded
        //);

        // next start is m_in_begin_load
        let buf_next_begin = 3 + m_in_begin_buf * byte_per_snv;
        //log::debug!("buf_next_begin {}", buf_next_begin);
        reader.seek(SeekFrom::Start(buf_next_begin as u64)).unwrap();
        // https://stackoverflow.com/questions/37079342/what-is-the-most-efficient-way-to-read-a-large-file-in-chunks-without-loading-th

        let m_in_buf = buf_num_snv.min(m_in_chrom - m_in_begin_buf);
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
            // use fill() might be fast?
            // This might be smaller than buf
            let buf_size_ = m_in_buf * byte_per_snv;
            buf.resize(buf_size_, 0);
            log::debug!("buf_size_: {}", buf_size_);
            assert_eq!(buf_size_ % byte_per_snv, 0);
            assert_eq!(buf.len(), buf_size_);

            // read_exact??
            reader.read_exact(&mut buf).unwrap();
            assert_eq!(buf.len(), buf_size_);
            log::debug!("loaded: {}", buf.len());

            assign_genot_buf(
                byte_per_snv,
                &mut buf,
                &mut g_snvs_file.as_genot_snvs_mut(m_begin_buf, m_end_buf),
                &use_snvs_snv,
                // m_to_m_in,
                use_samples,
                &mut g_group,
                group_to_m_in_buf.as_ref(),
                n,
                fill_missing_group,
            );
        }

        m_begin_buf = m_end_buf;
        m_in_begin_buf = m_in_end_buf;
        assert!(m_in_begin_buf <= m_in_chrom);
        if m_in_begin_buf == m_in_chrom {
            break;
        }
    }
    assert_eq!(m_in_begin_buf, m_in_chrom);
    buf
}

fn assign_genot_buf(
    byte_per_snv: usize,
    buf: &mut Vec<u8>,
    g_snv_buf: &mut GenotMut,
    // m_to_m_in: HashMap<usize, usize>,
    use_snvs_snv: &[bool],
    use_samples: Option<&[bool]>,
    g_group: &mut Option<&mut GenotMut>,
    group_to_m_in: Option<&Vec<Vec<usize>>>,
    // group_to_m_in: Option<&HashMap<usize, Vec<usize>>>,
    n: usize,
    fill_missing_group: Option<FillMissingGroup>,
) {
    // x plan 1 create vec of pred_use and  pred_use.iter()...
    //  -> unsafe when pred_use is duplicated
    // x plan 2 predictions.par_chunks_mut().filter(use)
    // o plan 3 predictions_chrom[m_buf..m_buf_end].par_chunks_mut()

    //log::debug!("m_in_end, m_end {},{}", m_in_end_loaded, m_end_loaded);

    // for snv m_i -> buf_i
    let (m_to_buf, _) = genotype::create_m_to_m_in(use_snvs_snv);

    g_snv_buf
        .iter_snv_mut()
        .enumerate()
        .par_bridge()
        .for_each(|(mi_loaded, mut g_snv)| {
            let buf_i = m_to_buf[&mi_loaded];
            // let m_in_read_i = m_to_m_in[&mi_loaded];
            let buf_mi = &buf[buf_i * byte_per_snv..(buf_i + 1) * byte_per_snv];

            assign_gsnv_from_bed(&mut g_snv, &buf_mi, use_samples);
        });

    // https://stackoverflow.com/questions/74763962/when-is-mut-required-when-passing-an-optionmut-t
    if let Some(g_group) = g_group {
        // m_in_i = buf_i in plink1

        let group_to_m_in = group_to_m_in.unwrap();
        g_group
            .iter_snv_mut()
            .enumerate()
            .par_bridge()
            .for_each_with(
                GenotSnv::new_empty(n),
                |g_snv_tmp, (group_i, mut g_group_snv)| {
                    let group_m_in_is = &group_to_m_in[group_i];
                    // let group_m_in_is = group_to_m_in.get(&group_i);

                    // if let Some(group_m_in_is) = group_m_in_is {
                    group_m_in_is.iter().for_each(|&group_m_in_i| {
                        let buf_i = group_m_in_i;
                        let buf_mi = &buf[buf_i * byte_per_snv..(buf_i + 1) * byte_per_snv];

                        // initialize g_snv_tmp
                        g_snv_tmp.fill_0();

                        assign_gsnv_from_bed(
                            &mut g_snv_tmp.as_genot_snv_mut_snv(),
                            &buf_mi,
                            use_samples,
                        );
                        match fill_missing_group {
                            Some(FillMissingGroup::Ref) => g_snv_tmp.fill_missing_ref(),
                            None => {}
                        };

                        g_group_snv.or_binary(&g_snv_tmp.as_genot_snv());
                    });
                    // }
                    // otherwise, do nothing
                },
            );
    }
}

// TODO: use SIMD
// can I implement buf: Option<&mut Vec<u8>> ?
// now unknown way to declare Vec<u8> with longer lifetime when buf=None
// What happens if no var is loaded in each loop? can I shorten time?
/// available for both whole and chunk
// fn assign_genot<R: BufRead + Seek>(
//     g_snvs_file: &mut GenotMut,
//     mut reader: R,
//     use_snvs: &[bool],
//     use_samples: Option<&[bool]>,
//     // cannot use &mut Vec<B8_2> since it cannot resize??
//     buf: Option<Vec<B8_2>>,
//     buf_size_max: Option<usize>,
//     // mut is necessary
//     // https://stackoverflow.com/questions/74763962/when-is-mut-required-when-passing-an-optionmut-t
//     mut g_group: Option<&mut GenotMut>,
//     group_to_m_in: Option<&HashMap<usize, Vec<usize>>>,
//     fill_missing_group: Option<FillMissingGroup>,
//     buf_num_snv: Option<usize>, // for testing
// ) -> Vec<B8_2> {
//     let m_in_chrom = use_snvs.len();

//     let n_in = match use_samples {
//         Some(v) => v.len(),
//         None => g_snvs_file.n(),
//     };

//     // FIXME: available mem?
//     // check if 32 GB? 16GB? remains or not
//     // reading smaller than buf_size is not error but reason is unknown
//     // [here](https://doc.rust-lang.org/std/io/trait.Read.html#tymethod.read)
//     // -> to deal with this, every time use Seek to change position

//     let byte_per_snv = bed_per_snv_size(n_in);
//     let buf_num_snv = match buf_num_snv {
//         Some(x) => x,
//         None => {
//             // already done in load_genot_files, but do again for mem_buf=None
//             let buf_size_max = buf_size_limit(buf_size_max);
//             // let byte_per_snv = bed_per_snv_size(n_in);
//             let buf_num_snv_limit: usize = buf_size_max / byte_per_snv;
//             buf_num_snv_limit.min(m_in_chrom)
//         }
//     };

//     // already done in load_genot_files, but do again for mem_buf=None
//     // let buf_size_max = buf_size_limit(buf_size_max);
//     // let byte_per_snv = bed_per_snv_size(n_in);
//     // let buf_num_snv_limit: usize = buf_size_max / byte_per_snv;
//     // let buf_num_snv: usize = buf_num_snv_limit.min(m_in_chrom);
//     let buf_size: usize = buf_num_snv * byte_per_snv;
//     assert_eq!(buf_size % byte_per_snv, 0);
//     // assert!(buf_size <= buf_size_max);

//     //let mut reader = BufReader::new(File::open(plink::fname_bed(fin_chrom, None)).unwrap());
//     // TODO: check buf_size < available

//     // first, create large buf here and resize later in loop
//     let mut buf: Vec<B8_2> = match buf {
//         Some(v) => v,
//         None => vec![0; buf_size],
//     };
//     // even if buf size < buf_size, resize later

//     let mut m_in_begin_loaded = 0;
//     let mut m_begin_loaded = 0;
//     // read buf length in one loop
//     loop {
//         //log::debug!(
//         //    "m_in_begin, m_begin {},{}",
//         //    m_in_begin_loaded, m_begin_loaded
//         //);

//         // next start is m_in_begin_load
//         let buf_next_begin = 3 + m_in_begin_loaded * byte_per_snv;
//         //log::debug!("buf_next_begin {}", buf_next_begin);
//         reader.seek(SeekFrom::Start(buf_next_begin as u64)).unwrap();
//         // https://stackoverflow.com/questions/37079342/what-is-the-most-efficient-way-to-read-a-large-file-in-chunks-without-loading-th

//         let m_in_read = buf_num_snv.min(m_in_chrom - m_in_begin_loaded);
//         log::debug!("m_in_read: {}", m_in_read);

//         let m_in_end_loaded = m_in_begin_loaded + m_in_read;
//         let use_snvs_loaded = &use_snvs[m_in_begin_loaded..m_in_end_loaded];

//         // also check group

//         let (m_to_m_in, m_read) = genotype::create_m_to_m_in(use_snvs_loaded_with_group);

//         let m_end_loaded = m_begin_loaded + m_read;

//         if m_read != 0 {
//             // use fill() might be fast?
//             // This might be smaller than buf
//             let buf_size_ = m_in_read * byte_per_snv;
//             buf.resize(buf_size_, 0);
//             log::debug!("buf_size_: {}", buf_size_);
//             assert_eq!(buf_size_ % byte_per_snv, 0);
//             assert_eq!(buf.len(), buf_size_);

//             // read_exact??
//             reader.read_exact(&mut buf).unwrap();
//             assert_eq!(buf.len(), buf_size_);
//             log::debug!("loaded: {}", buf.len());

//             // x plan 1 create vec of pred_use and  pred_use.iter()...
//             //  -> unsafe when pred_use is duplicated
//             // x plan 2 predictions.par_chunks_mut().filter(use)
//             // o plan 3 predictions_chrom[m_buf..m_buf_end].par_chunks_mut()

//             //log::debug!("m_in_end, m_end {},{}", m_in_end_loaded, m_end_loaded);

//             let mut g_file_part = g_snvs_file.as_genot_snvs_mut(m_begin_loaded, m_end_loaded);

//             g_file_part
//                 .iter_snv_mut()
//                 .enumerate()
//                 .par_bridge()
//                 .for_each(|(mi_loaded, mut g_snv)| {
//                     let m_in_read_i = m_to_m_in[&mi_loaded];
//                     let buf_mi = &buf[m_in_read_i * byte_per_snv..(m_in_read_i + 1) * byte_per_snv];

//                     assign_gsnv_from_bed(&mut g_snv, &buf_mi, use_samples);
//                 });

//             // https://stackoverflow.com/questions/74763962/when-is-mut-required-when-passing-an-optionmut-t
//             if let Some(g_group) = &mut g_group {
//                 // FIX: should be group_to_m_loaded
//                 let group_to_m_in = group_to_m_in.unwrap();

//                 let n = match use_samples {
//                     Some(v) => vec::count_true(v),
//                     None => n_in,
//                 };

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
//                                     let m_in_read_i = group_m_in_i;
//                                     // group_m_in_i is already m_in index
//                                     // let m_in_read_i = m_to_m_in[&group_m_in_i];
//                                     let buf_mi = &buf[m_in_read_i * byte_per_snv
//                                         ..(m_in_read_i + 1) * byte_per_snv];

//                                     // initialize g_snv_tmp
//                                     g_snv_tmp.fill_0();

//                                     assign_gsnv_from_bed(
//                                         &mut g_snv_tmp.as_genot_snv_mut_snv(),
//                                         &buf_mi,
//                                         use_samples,
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

//             // Error
//             //if let Some(x) = g_group {
//         }

//         m_begin_loaded = m_end_loaded;
//         m_in_begin_loaded = m_in_end_loaded;
//         assert!(m_in_begin_loaded <= m_in_chrom);
//         if m_in_begin_loaded == m_in_chrom {
//             break;
//         }
//     }
//     assert_eq!(m_in_begin_loaded, m_in_chrom);
//     buf
// }

/* // TODO: use SIMD
fn assign_predictions_whole<R: BufRead + Seek>(
    g_chrom: &mut GenotMut,
    mut reader: R,
    use_snvs: &[bool],
    use_samples: Option<&[bool]>,
) {
    //log::debug!("g_chrom len: {}", g_chrom.len());

    let n_in = if let Some(v) = use_samples {
        v.len()
    } else {
        g_chrom.n()
    };

    //plink::check_valid_bed(fin_chrom, None, m_in_chrom, n_in).unwrap();

    // seek is fast
    let file_size = reader.seek(SeekFrom::End(0)).unwrap() as usize;
    log::debug!("seek {}", file_size);
    let buf_size = file_size - 3;
    log::debug!("buf_size {}", buf_size);
    if buf_size == 0 {
        panic!("Buffer size is 0.");
    }

    let byte_per_snv = plink::bed_per_snv_size(n_in);
    assert_eq!(buf_size % byte_per_snv, 0);

    /*     log::debug!(
        "Loading buffer size vs available mem, {} bytes vs {} bytes",
        buf_size,
        alloc::get_available_memory()
    );
    if buf_size > alloc::get_available_memory() {
        panic!("Memory insufficient on loading bed file.");
    } */

    let mut buf: Vec<B8_2> = vec![0; buf_size];

    let current_pos = reader.seek(SeekFrom::Start(3)).unwrap();
    assert_eq!(current_pos, 3);
    reader.read_exact(&mut buf).unwrap();
    log::debug!("loaded: {}", buf.len());
    // do not use; somehow elongate buf size
    //let loaded_byte = reader.read_to_end(&mut buf).unwrap();
    //log::debug!("loaded: {} {}", loaded_byte, buf_size);
    //println!("loaded: {} {}", loaded_byte, buf_size);
    //println!("buflen {:?}", buf.len());
    //println!("buf {:?}", buf);

    //assert_eq!(loaded_byte, buf.len(), "loaded byte and buf len");
    //assert_eq!(loaded_byte, buf_size, "reader did not read whole buf size.");
    //assert_eq!(loaded_byte % byte_per_snv, 0);

    let (m_to_m_in, _) = genot_index::create_m_to_m_in(use_snvs);

    g_chrom
        .iter_snv_mut()
        .enumerate()
        .par_bridge()
        .for_each(|(mi_loaded, mut g_snv)| {
            let m_in_read_i = m_to_m_in[&mi_loaded];
            let buf_mi = &buf[m_in_read_i * byte_per_snv..(m_in_read_i + 1) * byte_per_snv];

            assign_pred_from_bed(&mut g_snv, &buf_mi, use_samples);
        });
} */

// // TODO: use SIMD
//// slower using .read() than .read_exact()
//// use when file size is too large
//fn assign_predictions_chunk_old2<R: BufRead + Seek>(
//    g_chrom: &mut GenotMut,
//    mut reader: R,
//    use_snvs: &[bool],
//    use_samples: Option<&[bool]>,
//    //buf: Option<Vec<B8_2>>,
//    mem: Option<usize>,
//    // mut is necessary
//    // https://stackoverflow.com/questions/74763962/when-is-mut-required-when-passing-an-optionmut-t
//    mut g_group: Option<&mut GenotMut>,
//    group_to_m_in: Option<&HashMap<usize, Vec<usize>>>,
//) {
//    let m_in_chrom = use_snvs.len();
//
//    let n_in = match use_samples {
//        Some(v) => v.len(),
//        None => g_chrom.n(),
//    };
//
//    //plink::check_valid_bed(fin_chrom, None, m_in_chrom, n_in).unwrap();
//
//    // check if 32 GB? 16GB? remains or not
//    // reading smaller than buf_size is not error but reason is unknown
//    // [here](https://doc.rust-lang.org/std/io/trait.Read.html#tymethod.read)
//    // -> to deal with this, every time use Seek to change position
//    //let buf_size_limit: usize = 64 * 1024 * 1024 * 1024;
//    let buf_size_limit: usize = 4 * 1024 * 1024 * 1024;
//
//    let byte_per_snv = bed_per_snv_size(n_in);
//    let buf_num_snv_limit: usize = buf_size_limit / byte_per_snv;
//    let buf_num_snv: usize = buf_num_snv_limit.min(m_in_chrom);
//    let buf_size: usize = buf_num_snv * byte_per_snv;
//    assert_eq!(buf_size % byte_per_snv, 0);
//    assert!(buf_size <= buf_size_limit);
//
//    //let mut reader = BufReader::new(File::open(plink::fname_bed(fin_chrom, None)).unwrap());
//    // TODO: check buf_size < available
//    let mut buf: Vec<B8_2> = vec![0; buf_size];
//
//    let mut m_in_begin_loaded = 0;
//    let mut m_begin_loaded = 0;
//    // read buf length in one loop
//    loop {
//        //log::debug!(
//        //    "m_in_begin, m_begin {},{}",
//        //    m_in_begin_loaded, m_begin_loaded
//        //);
//
//        // next start is m_in_begin_load
//        let buf_next_begin = 3 + m_in_begin_loaded * byte_per_snv;
//        //log::debug!("buf_next_begin {}", buf_next_begin);
//        reader.seek(SeekFrom::Start(buf_next_begin as u64)).unwrap();
//        // https://stackoverflow.com/questions/37079342/what-is-the-most-efficient-way-to-read-a-large-file-in-chunks-without-loading-th
//        // use fill() might be fast?
//        // This might be smaller than buf
//        // read_exact??
//        //let loaded_byte = reader.read_exact(&mut buf).unwrap();
//        let loaded_byte = reader.read(&mut buf).unwrap();
//        //if loaded_byte != buf_size, last buffer or could not load whole buffer
//        log::debug!("loaded: {} {}", loaded_byte, buf_size);
//        //}
//        if loaded_byte == 0 {
//            // all buf read
//            break;
//        } else {
//            //log::debug!("buf {:?}", buf.len());
//            // might fail when only part was loaded
//            //assert_eq!(loaded_byte % byte_per_snv, 0);
//            let m_in_read = loaded_byte / byte_per_snv;
//            // if m_in_read==0, loaded_byte < byte_per_snv, which means loaded_byte is less than one SNV.
//            assert!(m_in_read > 0, "loaded_byte is less than byte of one SNV.");
//
//            // x plan 1 create vec of pred_use and  pred_use.iter()...
//            //  -> unsafe when pred_use is duplicated
//            // x plan 2 predictions.par_chunks_mut().filter(use)
//            // o plan 3 predictions_chrom[m_buf..m_buf_end].par_chunks_mut()
//
//            let m_in_end_loaded = m_in_begin_loaded + m_in_read;
//            let use_snvs_loaded = &use_snvs[m_in_begin_loaded..m_in_end_loaded];
//
//            let (m_to_m_in, m_read) = genot_index::create_m_to_m_in(use_snvs_loaded);
//
//            // here for when m_read==0
//            // or move below and do not continue for m_read==0 since unnecessary
//            //m_in_begin_loaded = m_in_end_loaded;
//
//            //if m_read == 0 {
//            //    continue;
//            //}
//
//            let m_end_loaded = m_begin_loaded + m_read;
//
//            //log::debug!("m_in_end, m_end {},{}", m_in_end_loaded, m_end_loaded);
//
//            let mut g_chrom_part = g_chrom.as_genot_snvs_mut(m_begin_loaded, m_end_loaded);
//
//            g_chrom_part
//                .iter_snv_mut()
//                .enumerate()
//                .par_bridge()
//                .for_each(|(mi_loaded, mut g_snv)| {
//                    let m_in_read_i = m_to_m_in[&mi_loaded];
//                    let buf_mi = &buf[m_in_read_i * byte_per_snv..(m_in_read_i + 1) * byte_per_snv];
//
//                    assign_pred_from_bed(&mut g_snv, &buf_mi, use_samples);
//                });
//
//            m_begin_loaded = m_end_loaded;
//            // bug! when m_read==0, will not renewed
//            m_in_begin_loaded = m_in_end_loaded;
//        }
//    }
//    assert_eq!(m_in_chrom, m_in_begin_loaded);
//}

// error on rap
// rayon impossible since cannot read inside
fn _assign_genot_loadpart<R: BufRead + Seek>(
    g_chrom: &mut GenotMut,
    mut reader: R,
    use_snvs: &[bool],
    use_samples: Option<&[bool]>,
    // no use inside
    buf_in: Option<Vec<B8_2>>,
    _mem: Option<usize>,
    g_group: Option<&mut GenotMut>,
    //mut g_group: Option<&mut GenotMut>,
    group_to_m_in: Option<&HashMap<usize, Vec<usize>>>,
) -> Vec<B8_2> {
    if g_group.is_some() {
        unimplemented!("ny");
    }
    if group_to_m_in.is_some() {
        unimplemented!("ny");
    }

    //let m_in_chrom = use_snvs.len();

    // assume m and m_in in m_to_m_in are sorted
    // FIXME: use_snvs should be only chrom
    //let m_to_m_in = create_m_chrom_to_m_in_chrom(use_snvs);
    //let m_in_to_m = create_m_in_to_m_chrom(use_snvs, m_begin);

    let n_in = match use_samples {
        Some(v) => v.len(),
        None => g_chrom.n(),
    };

    // for debug
    //log::debug!("ys[0]:{:#010b}", ys[0]);

    //plink::check_valid_bed(fin_chrom, None, m_in_chrom, n_in).unwrap();

    let (m_to_m_in, m_read) = genotype::create_m_to_m_in(use_snvs);
    //let (m_to_m_in, m_for_read) = genot_index::create_m_to_m_in(use_snvs);

    let byte_per_snv = bed_per_snv_size(n_in);
    //let byte_per_snv = plink::bed_per_snv_size(n_in);
    let buf_size = byte_per_snv;

    //let mut reader = BufReader::new(File::open(plink::fname_bed(fin_chrom, None)).unwrap());
    let mut buf: Vec<B8_2> = vec![0; buf_size];
    assert_eq!(buf.len(), buf_size);

    assert_eq!(g_chrom.m(), m_read, "Sth wrong");

    // rayon impossible since cannot seek
    g_chrom
        .iter_snv_mut()
        .enumerate()
        .for_each(|(m_read_i, mut g_snv)| {
            log::debug!("m_read_i: {}", m_read_i);

            // for debug
            log::debug!("m_read_in: {}", m_to_m_in[&m_read_i]);
            let buf_begin = 3 + m_to_m_in[&m_read_i] * byte_per_snv;
            reader.seek(SeekFrom::Start(buf_begin as u64)).unwrap();
            reader.read_exact(&mut buf).unwrap();
            assert_eq!(buf.len(), buf_size);
            //let loaded_byte = reader.read(&mut buf).unwrap();
            //assert_eq!(loaded_byte, byte_per_snv);

            //let g_snv = g_chrom.to_genot_twin_snv_mut(m_read_i);
            //let pred = predictions_snv_s_mut(predictions_chrom, m_read_i, n);
            assign_gsnv_from_bed(&mut g_snv, &buf, use_samples);
        });

    // return buf
    // TODO: all right when calling assign_genot() after loadpart
    match buf_in {
        Some(v) => v,
        None => buf,
    }
}

//// rayon impossible since cannot read inside
//pub fn assign_genot_file_loadpart(
//    g_chrom: &mut GenotMut,
//    fin_chrom: &str,
//    use_snvs: &[bool],
//    use_samples: &[bool],
//) {
//    let m_in_chrom = use_snvs.len();
//    // log::debug!("use_samples: {:?}", use_samples);
//
//    // assume m and m_in in m_to_m_in are sorted
//    // FIXME: use_snvs should be only chrom
//    //let m_to_m_in = create_m_chrom_to_m_in_chrom(use_snvs);
//    //let m_in_to_m = create_m_in_to_m_chrom(use_snvs, m_begin);
//
//    let n_in = use_samples.len();
//
//    // for debug
//    //log::debug!("ys[0]:{:#010b}", ys[0]);
//
//    plink::check_valid_bed(fin_chrom, None, m_in_chrom, n_in).unwrap();
//
//    let (m_to_m_in, m_for_read) = genot_index::create_m_to_m_in(use_snvs);
//
//    let byte_per_snv = plink::bed_per_snv_size(n_in);
//    let buf_size = byte_per_snv;
//
//    let mut reader = BufReader::new(File::open(plink::fname_bed(fin_chrom, None)).unwrap());
//    let mut buf: Vec<B8_2> = vec![0; buf_size];
//
//    assert_eq!(g_chrom.m(), m_for_read, "Sth wrong");
//
//    // rayon impossible since cannot seek
//    g_chrom
//        .iter_snv_mut()
//        .enumerate()
//        .for_each(|(m_read_i, mut g_snv)| {
//            let buf_begin = 3 + m_to_m_in[&m_read_i] * byte_per_snv;
//            reader.seek(SeekFrom::Start(buf_begin as u64)).unwrap();
//            let loaded_byte = reader.read(&mut buf).unwrap();
//            assert_eq!(loaded_byte, byte_per_snv);
//
//            //let g_snv = g_chrom.to_genot_twin_snv_mut(m_read_i);
//            //let pred = predictions_snv_s_mut(predictions_chrom, m_read_i, n);
//            assign_pred_from_bed(&mut g_snv, &buf, use_samples);
//        });
//}
//
pub fn bed_per_snv_size(n: usize) -> usize {
    (n + 3) / 4
}

pub fn calculate_bed_size_genotype(m: usize, n: usize) -> usize {
    m * bed_per_snv_size(n)
    //m * ((n + 3) / 4)
}

pub fn calculate_bed_size(m: usize, n: usize) -> usize {
    3 + calculate_bed_size_genotype(m, n)
    //3 + m * ((n + 3) / 4)
}

//  TODO: should unwrap here?
/// return bed_size if valid, error otherwise
//pub fn check_valid_bed(fin: &str, n: usize, m: usize) -> Result<usize, String> {
pub fn check_valid_bed_file(
    //fin: &Path,
    //gfmt: GenotFormat,
    fin_genot: &GenotFile,
    //chrom: Option<&Chrom>,
    m: usize,
    n: usize,
) -> Result<usize, Box<dyn Error>> {
    let fin_bed = fin_genot.genotype_file();
    //let fin_bed = fname_plinks_genot(fin, gfmt, chrom);
    // check if open
    let mut reader = File::open(fin_bed)?;

    // check if size is correct
    let f_end: usize = reader.seek(SeekFrom::End(0)).unwrap() as usize;
    log::debug!("file end {}", f_end);
    let bed_size = calculate_bed_size(m, n);
    if f_end != bed_size {
        return Err(format!(
            "File size of .bed is wrong: {} vs correct {}.",
            f_end, bed_size
        )
        .into());
    }

    // check if the first 3 bytes are correct.
    reader.seek(SeekFrom::Start(0)).unwrap();
    let mut buf: Vec<u8> = Vec::with_capacity(n);
    unsafe {
        buf.set_len(3);
    }
    reader.read_exact(&mut buf).unwrap();
    //log::debug!("{:?}", buf);
    if buf != vec![0x6cu8, 0x1b, 0x01] {
        return Err("Magic number of .bed file is wrong.".into());
    }
    Ok(bed_size)
}

/// g_snv must be all false
/// use g_snv.fill_0() to make sure
pub fn assign_gsnv_from_bed(
    g_snv: &mut GenotSnvMut,
    buf_mi: &[B8_2],
    use_samples: Option<&[bool]>,
) {
    if let Some(use_sample) = use_samples {
        let mut ni = 0;
        for (n_in_i, v) in use_sample.iter().enumerate() {
            if *v {
                // bedcode
                let bcode = buf_to_ped_code(buf_mi, n_in_i);
                g_snv.set_bed_code_init_unchecked(bcode, ni);
                ni += 1;
            }
        }
    } else {
        for ni in 0..g_snv.n() {
            // bedcode
            let bcode = buf_to_ped_code(buf_mi, ni);
            g_snv.set_bed_code_init_unchecked(bcode, ni);
        }
    }
}

#[inline]
pub fn buf_to_ped_code(buf: &[B8_2], ni: usize) -> u8 {
    //buf[ni // 4], ni % 3
    byte_to_ped_code(buf[ni >> 2], ni & 3)
}

#[inline]
pub fn buf_to_count(buf: &[B8_2], ni: usize) -> u8 {
    //buf[ni // 4], ni % 3
    byte_to_count(buf[ni >> 2], ni & 3)
}

// {00: 2, 01: 3, 10: 1, 11: 0}
const CODE_TO_COUNT_AR: [u8; 4] = [2, 3, 1, 0];

// plink BED code to minor allele counts
// {00: 2, 01: 3, 10: 1, 11: 0}
#[inline]
fn code_to_count(v: u8) -> u8 {
    CODE_TO_COUNT_AR[v as usize]
}

#[inline]
fn byte_to_ped_code(v: B8_2, i: usize) -> u8 {
    (v & (0x03 << (i << 1))) >> (i << 1)
}

// extract ith code (= 2i~2i+1 bits from lower)
// ex. byte_to_count(0b11011000, 2) = code_to_count(0b10) = 1
#[inline]
fn byte_to_count(v: B8_2, i: usize) -> u8 {
    //#define count_pl(c, i) (code2count((((c) & (0x03 << (i << 1))) >> (i << 1))))
    code_to_count(byte_to_ped_code(v, i))
    //code_to_count((v & (0x03 << (i << 1))) >> (i << 1))
}

// TODO: this might be clear if use .flat_map()
// no rayon here
//fn load_x(x: &mut [u8], buf_mi: &[B8_2], use_samples: &[bool]) {
//    let mut ni = 0;
//    for (n_in_i, v) in use_samples.iter().enumerate() {
//        if *v {
//            x[ni] = buf_to_count(buf_mi, n_in_i);
//            ni += 1;
//        }
//    }
//    // ng: x could be larger thant n
//    //assert_eq!(ni, x.len());
//
//    super::missing_to_mode(x);
//}

// as fast as load_byte2 but complecated
#[allow(dead_code)]
fn load_byte1(fin: &Path, n: usize, m: usize) -> Vec<u8> {
    //let bed_size = calculate_bed_size(n, m);
    let v_size = calculate_bed_size_genotype(m, n);
    log::debug!("{}", v_size);
    let mut v: Vec<u8> = Vec::with_capacity(v_size);
    unsafe {
        v.set_len(v_size);
    }
    let mut reader = BufReader::new(File::open(fin).unwrap());
    let n: usize = 100_000_000;
    //let n: usize = 1_000_000;
    let mut buf: Vec<u8> = Vec::with_capacity(n);

    // first for 3
    unsafe {
        buf.set_len(3);
    }
    reader.read_exact(&mut buf).unwrap();
    log::debug!("{:?}", buf);

    unsafe {
        buf.set_len(n);
    }
    let mut loop_times = 0;
    // if we know the size is 9, then
    for i in 0..((v_size - 1) / n) {
        reader.read_exact(&mut buf).unwrap();
        // https://stackoverflow.com/questions/28219231/how-to-idiomatically-copy-a-slice
        v[n * i..n * (i + 1)].copy_from_slice(&buf);
        loop_times += 1;
    }
    log::debug!("loop times{}", loop_times);
    // for remaining
    let n_remain: usize = v_size - (v_size - 1) / n * n;
    unsafe {
        buf.set_len(n_remain);
    }
    reader.read_exact(&mut buf).unwrap();
    v[(v_size - 1) / n * n..v_size].copy_from_slice(&buf);
    v
}

// fast and simple
/// Use this function to load bytes.
#[allow(dead_code)]
fn load_byte2(fin: &Path, n: usize, m: usize) -> Vec<u8> {
    // you can use "seek"
    // use this!!
    //let bed_size = calculate_bed_size(n, m);
    let v_size = calculate_bed_size_genotype(m, n);
    let mut v: Vec<u8> = Vec::with_capacity(v_size);
    unsafe {
        v.set_len(v_size);
    }
    log::debug!("v_size: {}", v_size);

    let mut reader = BufReader::new(File::open(fin).unwrap());
    let n: usize = 100_000_000;
    //let n: usize = 1_000_000;
    let mut buf: Vec<u8> = Vec::with_capacity(n);

    // first for 3
    unsafe {
        buf.set_len(3);
    }
    reader.read_exact(&mut buf).unwrap();
    log::debug!("{:?}", buf);

    unsafe {
        buf.set_len(n);
    }

    let mut loop_times = 0;
    let mut i_ptr: usize = 0;
    loop {
        match reader.read(&mut buf).unwrap() {
            0 => break,
            n => {
                log::debug!("i_ptr,n:{},{}", i_ptr, n);
                let buf = &buf[..n];
                v[i_ptr..i_ptr + n].copy_from_slice(&buf);
                i_ptr = i_ptr + n;
            }
        }
        loop_times += 1;
    }
    log::debug!("loop times{}", loop_times);
    v
}

#[allow(dead_code)]
fn load_byte3(fin: &str, n: usize, m: usize) -> Vec<u8> {
    let mut file = File::open(fin).unwrap();
    // this is ok but not sure fast or slow
    let mut buf: Vec<u8> = Vec::new();
    // meaningless, lastly cap is 32
    //let n = 9;
    //let mut buf: Vec<u8> = Vec::with_capacity(n + 2);
    // this wont work
    //unsafe {
    //    buf.set_len(n);
    //}
    log::debug!("len: {}", buf.len());
    log::debug!("cap: {}", buf.capacity());
    let _ = file.read_to_end(&mut buf).unwrap();
    log::debug!("len: {}", buf.len());
    // larger than cap
    log::debug!("cap: {}", buf.capacity());
    log::debug!("buf[0]: {}", buf[0]);
    //log::debug!("{:?}", buf);

    //let v_size = calculate_bed_genotype_size(n, m);
    //tmp
    let v_size = calculate_bed_size(n, m);
    let mut v: Vec<u8> = Vec::with_capacity(v_size);
    unsafe {
        v.set_len(v_size);
    }

    (&mut v).copy_from_slice(&buf);
    v
}

pub fn run_byte(fin: &Path, n: usize, m: usize) -> Vec<u8> {
    // fast but complecated
    log::debug!("way 1");
    let istart = Instant::now();
    let v = load_byte1(fin, n, m);
    log::debug!("load_byte1: {:?}", Instant::now().duration_since(istart));
    log::debug!("{}", v[0]);

    // as fast as way 1
    // also easy to write
    log::debug!("way 2");
    let istart = Instant::now();
    let v = load_byte2(fin, n, m);
    log::debug!("load_byte2: {:?}", Instant::now().duration_since(istart));
    log::debug!("{}", v[0]);

    /*
    // seems slow -> could be because file size was larger than mem size?
    log::debug!("way 3");
    let istart = Instant::now();
    let v = load_byte3(fin, n, m);
    log::debug!("load_byte3: {:?}", Instant::now().duration_since(istart));
    log::debug!("{}", v[0]);
    */

    v
}

#[cfg(test)]
mod tests {

    use super::*;
    use crate::samples;
    use crate::{sample, snv};
    use std::{io::Cursor, path::PathBuf};

    fn setup_test() -> (GenotFiles, Vec<bool>, usize, usize, Vec<bool>, Vec<bool>) {
        let fin = PathBuf::from("../../test/data/toy1/genot");
        //let gfmt = GenotFormat::Plink1;
        let fin_genot = GenotFiles::new_plink1(vec![fin]);
        //let fin_snv = None;
        //let fin_sample = None;

        let m_in: usize = genot_io::compute_num_snv(&fin_genot);
        log::debug!("{}", m_in);
        //let n_in: usize = genot_io::compute_num_sample(&fin_genot).unwrap();
        //println!("{}", n_in);
        // load snvs
        let (snvs_in, _) = genot_io::load_snvs(&fin_genot);
        let (use_snvs, m) = snv::make_use_snvs_buf(None, &snvs_in);
        //let (m, use_snvs) = snv::make_use_snvs_buf(None, &snvs_in);
        //let (m, use_snvs) = snv::make_use_snvs(fin_snv, &snvs_in);
        //let (m,use_snvs: Vec<bool>) = plink::make_use_snv(fin, snvs_in);
        let (use_samples, n) = sample::make_use_samples_buf(None, &fin_genot);
        //let (n, use_samples) = sample::make_use_samples(fin_sample, &fin_genot);
        //let ys: Vec<bool> = io_genot::load_ys(&fin, gfmt, None, None, &use_samples);
        let sample_id_to_n = samples::create_sample_id_to_n(&fin_genot, Some(&use_samples));
        let ys: Vec<bool> = genot_io::load_ys_buf(&fin_genot, None, None, &sample_id_to_n).unwrap();

        (fin_genot, ys, m, n, use_snvs, use_samples)
    }

    // #[test]
    // fn test_group_to_m_in_chrom() {
    //     let group_to_m_in: HashMap<usize, Vec<usize>> = HashMap::from_iter(vec![
    //         (0, vec![0, 1]),
    //         (1, vec![2, 3]),
    //         (2, vec![4, 5]),
    //         (3, vec![6, 7]),
    //         (4, vec![8, 9]),
    //         (5, vec![10, 5]),
    //     ]);

    //     let m_in_begin = 3;
    //     let m_in_end = 7;

    //     let group_to_m_in_chrom = group_to_m_in_file(Some(&group_to_m_in), m_in_begin, m_in_end);

    //     // extracted in the range
    //     // (1, [3])
    //     // (2, [4, 5])
    //     // (3, [6])
    //     // (5, [5])
    //     //
    //     // after adjusting for m_in_begin
    //     // (1, [0])
    //     // (2, [1, 2])
    //     // (3, [3])
    //     // (5, [2])
    //     //
    //     let group_to_m_in_chrom_ans: HashMap<usize, Vec<usize>> = HashMap::from_iter(vec![
    //         (1, vec![0]),
    //         (2, vec![1, 2]),
    //         (3, vec![3]),
    //         (5, vec![2]),
    //     ]);

    //     assert_eq!(group_to_m_in_chrom.unwrap(), group_to_m_in_chrom_ans);
    // }

    #[test]
    fn test_generate_genot_plink_whole() {
        let (fin_genot, _, m, n, use_snvs, use_samples) = setup_test();
        //let gfmt = GenotFormat::Plink1;
        //let use_snvs = vec![true; use_snvs.len()];
        //let use_samples = vec![true; use_samples.len()];
        let g = generate_genot_plink1(
            &fin_genot,
            //gfmt,
            m,
            n,
            Some(&use_snvs),
            None,
            Some(&use_samples),
            None,
            // false,
            None,
            None,
        );
        //let g = generate_genot_plink(&fin, gfmt, m, n, Some(&use_snvs), Some(&use_samples), true);
        let mut giter = g.iter_snv();
        assert_eq!(giter.next().unwrap().vals(), vec![2, 0, 1, 0, 0]);
        assert_eq!(giter.next().unwrap().vals(), vec![1, 0, 2, 1, 0]);
        assert_eq!(giter.next().unwrap().vals(), vec![0, 2, 0, 1, 1]);
    }

    /*     #[test]
    fn test_generate_predictions_whole2() {
        let (fin, _, m, n, use_snvs, use_samples) = setup_test2();
        //let use_snvs = vec![true; use_snvs.len()];
        //let use_samples = vec![true; use_samples.len()];
        generate_predictions(&fin, m, n, &use_snvs, Some(&use_samples), true);
    } */

    #[test]
    fn test_assign_gsnv_from_bed() {
        let mut g = GenotSnv::new_empty(4);
        // [2, 0, 3, 0, 1, 0]
        let pbuf = vec![0b11_01_11_00, 0b00_00_11_10];

        let use_samples = vec![true, false, true, true, true, false];
        assign_gsnv_from_bed(&mut g.as_genot_snv_mut_snv(), &pbuf, Some(&use_samples));
        assert_eq!(g.vals(), vec![2, 3, 0, 1]);
    }

    #[test]
    fn test_assign_genot() {
        let mut g = Genot::new_zeros(2, 3);

        let reader: Vec<u8> = vec![0x6c, 0x1b, 0x01, 0b00_01_11_00, 0b00_10_00_11];
        let mut cur = Cursor::new(reader.as_slice());
        load_genot_file(
            &mut g.as_genot_mut(),
            &mut cur,
            //reader.as_slice(),
            &[true, true],
            None,
            None,
            None,
            None,
            None,
            None,
            Some(1),
        );

        assert_eq!(g.vals_snv(0), vec![2, 0, 3]);
        assert_eq!(g.vals_snv(1), vec![0, 2, 1]);
    }

    #[test]
    fn test_assign_genot_group() {
        let mut g_snvs = Genot::new_zeros(3, 5);
        let mut g_group = Genot::new_zeros(2, 5);

        let group_m_to_m_in: Vec<Vec<usize>> = vec![vec![0, 1], vec![0, 1, 3]];
        //HashMap::from_iter(vec![(1, vec![0, 1, 2])]);

        // snv0: [0, 0, 0, 1, 3]
        // snv1: [0, 0, 1, 0, 0]
        // snv2: [2, 2, 2, 2, 2] // no use
        // snv3: [0, 1, 2, 3, 1]
        //
        // group0: snv0,1
        // group1:  snv0,1,3

        // after fill missing to ref
        // snv0: [0, 0, 0, 1, 0]
        // snv1: [0, 0, 1, 0, 0]
        // snv2: [2, 2, 2, 2, 2] // no use
        // snv3: [0, 1, 2, 0, 1]
        //
        // group would be
        // group0: [0, 0, 1, 1, 0]
        // group1:  [0, 1, 1, 1, 1]
        let reader: Vec<u8> = vec![
            0x6c,
            0x1b,
            0x01,
            0b10_11_11_11,
            0b00_00_00_01,
            0b11_10_11_11,
            0b00_00_00_11,
            0b11_11_11_11,
            0b11_11_11_11,
            0b01_00_10_11,
            0b00_00_00_10,
        ];
        let mut cur = Cursor::new(reader.as_slice());
        load_genot_file(
            &mut g_snvs.as_genot_mut(),
            &mut cur,
            &[true, true, false, true],
            None,
            None,
            None,
            Some(&mut g_group.as_genot_mut()),
            Some(&group_m_to_m_in),
            Some(FillMissingGroup::Ref),
            Some(2),
            // None,
        );

        assert_eq!(g_snvs.vals_snv(0), vec![0, 0, 0, 1, 3]);
        assert_eq!(g_snvs.vals_snv(1), vec![0, 0, 1, 0, 0]);
        assert_eq!(g_snvs.vals_snv(2), vec![0, 1, 2, 3, 1]);
        assert_eq!(g_group.vals_snv(0), vec![0, 0, 1, 1, 0]);
        assert_eq!(g_group.vals_snv(1), vec![0, 1, 1, 1, 1]);
    }

    /*
    #[test]
    fn test_load_snv_read() {
        let mut g = GenotSnv::new_empty(3);

        let reader: Vec<u8> = vec![0x6c, 0x1b, 0x01, 0b00_01_11_00, 0b00_10_00_11];
        let mut cur = Cursor::new(reader.as_slice());
        load_snv_read(
            &mut g.as_genot_snv_mut_snv(),
            &mut cur,
            //reader.as_slice(),
            1,
            Some(&[true, true, true]),
            3,
        );

        assert_eq!(g.vals(), vec![0, 2, 1]);
    }
    */

    #[test]
    fn test_get_bed_size() {
        let m: usize = 3;
        let n: usize = 5;
        let bed_size: usize = calculate_bed_size(m, n);
        assert_eq!(bed_size, 9);

        let m: usize = 3;
        let n: usize = 4;
        let bed_size: usize = calculate_bed_size(m, n);
        assert_eq!(bed_size, 6);
    }

    #[test]
    fn test_check_valid_bed() {
        let fin = PathBuf::from("../../test/data/toy1/genot");
        //let gfmt = GenotFormat::Plink1;
        //let fins_genot = GenotFiles::new_plink1(vec![fin]);
        let fin_genot = GenotFile::Plink1(fin);

        let m: usize = genot_io::compute_num_snv_file_bi_allelic(&fin_genot);
        //let m: usize = genot_io::compute_num_snv_file(&fin_genot);
        let n: usize = genot_io::compute_num_sample_file(&fin_genot);
        log::debug!("m,n: {},{}", m, n);
        let bed_size = check_valid_bed_file(&fin_genot, m, n).unwrap();
        assert_eq!(bed_size, calculate_bed_size(m, n));
    }

    #[test]
    #[should_panic]
    fn test_check_valid_bed_panic_nofile() {
        // file does not exist
        let fin = PathBuf::from("./toy");
        let fin_genot = GenotFile::Plink1(fin);

        check_valid_bed_file(&fin_genot, 3, 2).unwrap();
    }

    #[test]
    #[should_panic]
    fn test_check_valid_bed_panic() {
        let fin = PathBuf::from("../../test/data/toy1/genot");
        let fin_genot = GenotFile::Plink1(fin);

        let m: usize = genot_io::compute_num_snv_file_bi_allelic(&fin_genot);
        //let m: usize = genot_io::compute_num_snv_file(&fin_genot);
        let n: usize = genot_io::compute_num_sample_file(&fin_genot);
        // size is wrong
        check_valid_bed_file(&fin_genot, m, n - 1).unwrap();
    }

    #[test]
    fn test_code_to_count() {
        let vs_input = vec![0b00, 0b01, 0b10, 0b11];
        let expects = vec![2, 3, 1, 0];
        for (v, exp) in vs_input.iter().zip(expects.iter()) {
            assert_eq!(code_to_count(*v), *exp);
        }
    }

    #[test]
    #[should_panic]
    fn test_code_to_count_panic() {
        let v = 4;
        code_to_count(v);
    }

    #[test]
    fn test_byte_to_count() {
        let v = 0b00_01_10_11;
        let expects = vec![0, 1, 3, 2];

        for (i, exp) in expects.iter().enumerate() {
            assert_eq!(byte_to_count(v, i), *exp);
        }
    }

    #[test]
    fn test_buf_to_count() {
        let buf = vec![0b00_01_10_11, 0b00_00_11_01];
        // read the lowest of the first byte
        assert_eq!(buf_to_count(&buf, 0), 0);
        // read the higest of the first byte
        assert_eq!(buf_to_count(&buf, 3), 2);
        // read the lowest of the second byte
        assert_eq!(buf_to_count(&buf, 4), 3);
        // read the highest of the second byte
        assert_eq!(buf_to_count(&buf, 7), 2);
    }

    #[test]
    fn test_load_byte2() {
        let fin = PathBuf::from("../../test/data/toy1/genot");
        //let gfmt = GenotFormat::Plink1;
        let fin_genot = GenotFile::Plink1(fin);

        let fin_fam = fin_genot.sample_file();
        //let fin_fam = genot_io::fname_plinks_sample(&fin, gfmt, None);
        //let fin_fam = fin.clone() + ".fam";
        let n: usize = textfile::compute_num_line_text(&fin_fam, None).unwrap();

        let fin_bim = fin_genot.snv_file();
        //let fin_bim = genot_io::fname_plinks_snv(&fin, gfmt, None);
        //let fin_bim = fin.clone() + ".bim";
        let m: usize = textfile::compute_num_line_text(&fin_bim, None).unwrap();

        let fin_bed = fin_genot.genotype_file();
        //let fin_bed = genot_io::fname_plinks_genot(&fin, gfmt, None);
        //let fin_bed = fin.clone() + ".bed";
        load_byte2(&fin_bed, n, m);
    }
}
