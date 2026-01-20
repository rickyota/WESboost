pub mod plink1;
#[cfg(feature = "plink2")]
pub mod plink2;

use rayon::prelude::*;
use std::collections::HashMap;

use crate::genot::prelude::*;
use crate::{genot_io, vec, FillMissing, FillMissingGroup, GenotFile, GenotFiles};

/*
// This is not fast for plink2
//  unused now: for small number of snvs called from python
//  mi for multiple files
/// mi is in fin file
pub fn generate_genot_snv(
    fin: &Path,
    gfmt: GenotFormat,
    mi: usize,
    n: usize,
    use_samples: Option<&[bool]>,
    //use_missing: bool,
    fill_missing: bool,
) -> GenotSnv {
    match gfmt {
        GenotFormat::Plink1 => {
            plink::generate_genot_snv_plink(fin, gfmt, mi, n, use_samples, fill_missing)
        }
        GenotFormat::Plink2 | GenotFormat::Plink2Vzs => {
            call_generate_genot_snv_plink2(fin, gfmt, mi, n, use_samples, fill_missing)
            //if cfg!(feature = "plink2") {
            //    plink2::load_snv_plink2(fin, gfmt, mi, n, use_samples, use_missing)
            //} else {
            //    panic!("Cannot use plink2 in this program feature. Use --feature plink2");
            //}
        }
    }
}

#[cfg(feature = "plink2")]
fn call_generate_genot_snv_plink2(
    fin: &Path,
    gfmt: GenotFormat,
    mi: usize,
    n: usize,
    use_samples: Option<&[bool]>,
    fill_missing: bool,
) -> GenotSnv {
    plink2::generate_genot_snv_plink2(fin, gfmt, mi, n, use_samples, fill_missing)
}

#[cfg(not(feature = "plink2"))]
fn call_generate_genot_snv_plink2(
    _: &Path,
    _: GenotFormat,
    _: usize,
    _: usize,
    _: Option<&[bool]>,
    _: bool,
) -> GenotSnv {
    panic!("Cannot use plink2 in this program feature. Use --feature plink2");
}
*/

/// Generate Vector of predictions.
pub fn generate_genot_simple(
    fins_genot: &GenotFiles,
    use_snvs: Option<&[bool]>,
    use_samples: Option<&[bool]>,
    fill_missing: Option<FillMissing>,
    //fill_missing_mode: bool,
    m: Option<usize>,
    n: Option<usize>,
    mem: Option<usize>,
) -> Genot {
    let m = match m {
        Some(x) => x,
        None => match use_snvs {
            Some(x) => vec::count_true(x),
            None => genot_io::compute_num_snv(&fins_genot),
        },
    };

    let n = match n {
        Some(x) => x,
        None => match use_samples {
            Some(x) => vec::count_true(x),
            None => genot_io::compute_num_sample(fins_genot),
        },
    };

    let (_, file_snv_allele_idx_in) = genot_io::load_snvs(&fins_genot);

    generate_genot(
        fins_genot,
        m,
        n,
        use_snvs,
        &file_snv_allele_idx_in,
        None,
        use_samples,
        fill_missing,
        //fill_missing_mode,
        None,
        mem,
    )
}

pub fn generate_genot(
    fins_genot: &GenotFiles,
    m_snv: usize,
    n: usize,
    // snv only. not for group
    use_snvs: Option<&[bool]>,
    file_snv_allele_idx_in: &[(usize, usize, usize)],
    group_to_m_in: Option<Vec<Vec<usize>>>,
    // group_to_m_in: Option<HashMap<usize, Vec<usize>>>,
    use_samples: Option<&[bool]>,
    fill_missing: Option<FillMissing>,
    fill_missing_group: Option<FillMissingGroup>,
    //fill_missing_mode: bool,
    mem: Option<usize>,
) -> Genot {
    // assume all files are in the same format, which is checked in GenotFiles::check()
    match fins_genot.files()[0] {
        GenotFile::Plink1(_) => plink1::generate_genot_plink1(
            fins_genot,
            m_snv,
            n,
            use_snvs,
            group_to_m_in,
            use_samples,
            fill_missing,
            fill_missing_group,
            //fill_missing_mode,
            mem,
        ),
        GenotFile::Plink2(_) | GenotFile::Plink2Vzs(_) => {
            call_generate_genot_plink2(
                fins_genot,
                m_snv,
                n,
                use_snvs,
                file_snv_allele_idx_in,
                group_to_m_in,
                use_samples,
                fill_missing,
                fill_missing_group,
                //fill_missing_mode,
                mem,
            )
        }
    }
}

#[cfg(feature = "plink2")]
fn call_generate_genot_plink2(
    fins_genot: &GenotFiles,
    m_snv: usize,
    n: usize,
    use_snvs: Option<&[bool]>,
    file_snv_allele_idx_in: &[(usize, usize, usize)],
    group_to_m_in: Option<Vec<Vec<usize>>>,
    // group_to_m_in: Option<HashMap<usize, Vec<usize>>>,
    use_samples: Option<&[bool]>,
    fill_missing: Option<FillMissing>,
    //fill_missing_mode: bool,
    fill_missing_group: Option<FillMissingGroup>,
    mem: Option<usize>,
) -> Genot {
    plink2::generate_genot_plink2(
        fins_genot,
        m_snv,
        n,
        use_snvs,
        file_snv_allele_idx_in,
        group_to_m_in,
        use_samples,
        fill_missing,
        //fill_missing_mode,
        fill_missing_group,
        mem,
        None,
    )
}

#[cfg(not(feature = "plink2"))]
fn call_generate_genot_plink2(
    _: &GenotFiles,
    _: usize,
    _: usize,
    _: Option<&[bool]>,
    _: &[(usize, usize, usize)],
    _: Option<Vec<Vec<usize>>>,
    // _: Option<HashMap<usize, Vec<usize>>>,
    _: Option<&[bool]>,
    _: Option<FillMissing>,
    _: Option<FillMissingGroup>,
    // _: bool,
    _: Option<usize>,
) -> Genot {
    panic!("Cannot use plink2 in this program feature. Use --feature plink2");
}

/* // error: plink2 is not loaded without feature=plink2
fn call_generate_genot_plink2(
    fin: &Path,
    gfmt: GenotFormat,
    m: usize,
    n: usize,
    use_snvs: Option<&[bool]>,
    use_samples: Option<&[bool]>,
    use_missing: bool,
) -> Genot {
    if cfg!(feature = "plink2") {
        plink2::generate_genot_plink2(fin, gfmt, m, n, use_snvs, use_samples, use_missing)
    } else {
        panic!("Cannot use plink2 in this program feature. Use --feature plink2");
    }
} */

pub fn group_to_m_in_range(
    group_to_m_in: Option<&Vec<Vec<usize>>>,
    // group_to_m_in: Option<&HashMap<usize, Vec<usize>>>,
    m_in_begin: usize,
    m_in_end: usize,
) -> Option<Vec<Vec<usize>>> {
    // ) -> Option<HashMap<usize, Vec<usize>>> {
    if group_to_m_in.is_none() {
        return None;
    }

    let group_to_m_in = group_to_m_in.unwrap();

    let mut group_to_m_in_range = Vec::with_capacity(group_to_m_in.len());
    for m_in_group in group_to_m_in.iter() {
        let m_in_group_range: Vec<usize> = m_in_group
            .iter()
            .filter(|&mi| m_in_begin <= *mi && *mi < m_in_end)
            .map(|&mi| mi - m_in_begin)
            .collect();
        // add even if m_in_group is empty
        group_to_m_in_range.push(m_in_group_range);
    }

    // for (group_i, m_in_set) in group_to_m_in.iter() {
    //     let m_in_set_chrom: Vec<usize> = m_in_set
    //         .iter()
    //         .filter(|&mi| m_in_begin <= *mi && *mi < m_in_end)
    //         .map(|&mi| mi - m_in_begin)
    //         .collect();
    //     if !m_in_set_chrom.is_empty() {
    //         group_to_m_in_file.insert(*group_i, m_in_set_chrom);
    //     }
    // }

    Some(group_to_m_in_range)
}

/// fill missing for g_snv only
pub fn fill_missing_g_snvs(g_snvs: &mut GenotMut, fill_missing: Option<FillMissing>) {
    if fill_missing.is_some() {
        match fill_missing.unwrap() {
            FillMissing::Mode => {
                g_snvs
                    .iter_snv_mut()
                    .par_bridge()
                    .for_each(|mut g_snv| g_snv.fill_missing_mode());
            }
            FillMissing::Ref => {
                g_snvs
                    .iter_snv_mut()
                    .par_bridge()
                    .for_each(|mut g_snv| g_snv.fill_missing_ref());
            }
        }
    }

    // no use
    // print missing rate for group
    // if let Some(g_group) = &mut g_group {
    // let group_missing_count: Vec<usize> = g_group
    //    .iter_snv()
    //    .par_bridge()
    //    .map(|g_snv| g_snv.count_missing())
    //    .collect();

    // let (hist_missing, hist_bins) = vec::histogram(&group_missing_count, 10);
    // log::info!("Histgram of missing count for group");
    // hist_bins
    //    .iter()
    //    .zip(hist_missing.iter())
    //    .for_each(|(bin, count)| {
    //        log::info!("<{}: {}", bin, count);
    //    });
    // log::info!(
    //    ">={}: {}",
    //    hist_bins.last().unwrap(),
    //    hist_missing.last().unwrap()
    // );
    // }
}

fn convert_group_to_use(
    group_to_m_in: Option<&Vec<Vec<usize>>>,
    // group_to_m_in: Option<&HashMap<usize, Vec<usize>>>,
    m_in_buf: usize,
) -> Vec<bool> {
    let use_snvs_group: Vec<bool> = match &group_to_m_in {
        Some(group_to_m_in) => {
            let mut use_snvs_group = vec![false; m_in_buf];
            group_to_m_in.iter().for_each(|m_in_is| {
                m_in_is.iter().for_each(|&m_in_i| {
                    use_snvs_group[m_in_i] = true;
                });
            });
            use_snvs_group
        }
        None => vec![false; m_in_buf],
    };
    use_snvs_group
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_group_to_m_in_chrom() {
        let group_to_m_in: Vec<Vec<usize>> = vec![
            vec![0, 1],
            vec![2, 3],
            vec![4, 5],
            vec![6, 7],
            vec![8, 9],
            vec![10, 5],
        ];
        // let group_to_m_in: HashMap<usize, Vec<usize>> = HashMap::from_iter(vec![
        //     (0, vec![0, 1]),
        //     (1, vec![2, 3]),
        //     (2, vec![4, 5]),
        //     (3, vec![6, 7]),
        //     (4, vec![8, 9]),
        //     (5, vec![10, 5]),
        // ]);

        let m_in_begin = 3;
        let m_in_end = 7;

        let group_to_m_in_chrom = group_to_m_in_range(Some(&group_to_m_in), m_in_begin, m_in_end);

        // extracted in the range
        // (1, [3])
        // (2, [4, 5])
        // (3, [6])
        // (5, [5])
        //
        // after adjusting for m_in_begin
        // (1, [0])
        // (2, [1, 2])
        // (3, [3])
        // (5, [2])
        //
        // let group_to_m_in_chrom_ans: HashMap<usize, Vec<usize>> = HashMap::from_iter(vec![
        let group_to_m_in_chrom_ans: Vec<Vec<usize>> =
            vec![vec![], vec![0], vec![1, 2], vec![3], vec![], vec![2]];

        assert_eq!(group_to_m_in_chrom.unwrap(), group_to_m_in_chrom_ans);
    }
}
