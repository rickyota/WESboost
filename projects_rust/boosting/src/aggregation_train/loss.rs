// use core::num;
use std::fs::File;
use std::io::{BufWriter, Write};

//use cmatrix::dense::{BaseCBits, BaseCMatrixMut, CBitsMut};
//use genetics::genot::prelude::*;
use cmatrix::prelude::*;
use genetics::vec;
use genetics::{regression, BaseGenot, Genot, Snvs};
use genetics::{samples::prelude::*, AggId};

pub fn error_rate(gbits: &CBitsRef, phe: &Phe, case_n: usize) -> f64 {
    let (_d1, n1, d0, _n0) = gbits.stat_contingency_table(phe.inner(), case_n);

    let error_rate = (d0 + n1) as f64 / phe.n() as f64;

    // should be slow
    // let mut tp = 0usize;
    // let mut fp_ = 0usize;

    // gbits.iter().zip(phe.iter()).for_each(|(gbit, y)| {
    //     if gbit {
    //         if y {
    //             tp += 1;
    //         } else {
    //             fp_ += 1;
    //         }
    //     }
    // });

    // let tp = tp as f64;
    // let fp_ = fp_ as f64;
    // let fn_ = case_n as f64 - tp;
    // let error_rate = (fp_ + fn_) / phe.n() as f64;
    return error_rate;
}

// return f64::NaN if any of table is 0
pub fn logistic_loss(gbits: &CBitsRef, phe: &Phe, case_n: usize) -> f64 {
    let (d1, n1, d0, n0) = gbits.stat_contingency_table(phe.inner(), case_n);

    let (d1, n1, d0, n0) = (d1 as f64, n1 as f64, d0 as f64, n0 as f64);

    if d1 == 0.0 || n1 == 0.0 || d0 == 0.0 || n0 == 0.0 {
        return f64::NAN;
    }

    // loss = - \sum_k n_k ln n_k + num_0 ln num_0 + num_1 ln num_1
    let loss = -d1 * d1.ln() - n1 * n1.ln() - d0 * d0.ln() - n0 * n0.ln()
        + (d1 + n1) * (d1 + n1).ln()
        + (d0 + n0) * (d0 + n0).ln();

    loss
}

// ng: always NaN
// pub fn logistic_loss_0(phe: &Phe, case_n: usize) -> f64 {
//     let gbits_0 = CBits::new_zeros(phe.n());
//     logistic_loss(&gbits_0.as_cbits_ref_b(), phe, case_n)
// }

//pub fn exp_loss() {}

/// log10(P) not -log10(P)
/// the smaller, the better
/// TODO: f64::NAN not -f64::NAN
pub fn logistic_log10_pval(gbits: &CBitsRef, phe: &Phe, case_n: usize) -> f64 {
    -logistic_mlog10_pval(gbits, phe, case_n)
}

/// -log10(P)
pub fn logistic_mlog10_pval(gbits: &CBitsRef, phe: &Phe, case_n: usize) -> f64 {
    match logistic_mlog10_pval_z_score_option(gbits, phe, case_n) {
        Some(x) => x.0,
        None => f64::NAN,
    }
}

pub fn logistic_mlog10_pval_z_score(gbits: &CBitsRef, phe: &Phe, case_n: usize) -> (f64, f64) {
    match logistic_mlog10_pval_z_score_option(gbits, phe, case_n) {
        Some(x) => x,
        None => (f64::NAN, f64::NAN),
    }
}

pub fn logistic_mlog10_pval_z_score_option(
    gbits: &CBitsRef,
    phe: &Phe,
    case_n: usize,
) -> Option<(f64, f64)> {
    let (d1, n1, d0, n0) = gbits.stat_contingency_table(phe.inner(), case_n);

    let observed = vec![vec![d1 as f64, d0 as f64], vec![n1 as f64, n0 as f64]];

    regression::logistic_regression_cont_table_2_2_mlog10_pval_zscore(&observed)
}

// error: digit loss
/// log10(P) not -log10(P)
/// the smaller, the better
// pub fn logistic_log10_pval(gbits: &CBitsRef, phe: &Phe, case_n: usize) -> f64 {
//     match logistic_pval(gbits, phe, case_n) {
//         Some(x) => x.log10(),
//         None => f64::NAN,
//     }
// }

// error: digit loss
// pub fn logistic_pval(gbits: &CBitsRef, phe: &Phe, case_n: usize) -> Option<f64> {
//     let (d1, n1, d0, n0) = gbits.stat_contingency_table(phe.inner(), case_n);

//     let observed = vec![vec![d1 as f64, d0 as f64], vec![n1 as f64, n0 as f64]];

//     regression::logistic_regression_cont_table_2_2_p_value(&observed)
// }

// pub fn chi_square_pval(gbits: &CBitsRef, phe: &Phe, case_n: usize) -> f64 {
//     let (d1, n1, d0, n0) = gbits.stat_contingency_table(phe.inner(), case_n);

//     let observed = vec![vec![d1 as f64, d0 as f64], vec![n1 as f64, n0 as f64]];

//     regression::chi_square_cont_table_p_value(&observed)
// }

// wrong code
//// assume no missing
//// TODO: confirm by confirm_no_missing();
//pub fn error_rate_tp_fn_add(
//    //tp_fn: (usize, usize),
//    prev: &CBitsRef,
//    gbits: &CBitsRef,
//    phe: &Phe,
//) -> (usize, usize) {
//    // TODO: simd ver.
//
//    let mut new = prev.clone_to_cbits();
//    let mut tp_fn_add = (0usize, 0);
//    //let mut tp_fn_add = tp_fn.clone();
//
//    prev.iter()
//        .zip(gbits.iter())
//        .zip(phe.iter())
//        .enumerate()
//        .for_each(|(i, ((prev, gbit), y))| {
//            // prev=false, gbit=true is only pattern that we need to update
//            //if !prev && gbit {
//            if gbit && !prev {
//                new.set_bool_unchecked_b(true, i);
//                if y {
//                    tp_fn_add.0 += 1;
//                } else {
//                    tp_fn_add.1 += 1;
//                }
//            }
//        });
//    return tp_fn_add;
//}

#[derive(Debug, Clone)]
pub enum AggLossStruct {
    // vec for each region
    // min loss of a region, (pos_start, pos_end)
    Loss(Vec<(f64, (usize, usize))>),
    // Loss(Vec<f64>, Vec<(usize, usize)>),
}

impl AggLossStruct {
    pub fn inner_mut(&mut self) -> &mut Vec<(f64, (usize, usize))> {
        match self {
            AggLossStruct::Loss(x) => x,
        }
    }

    pub fn new(m: usize) -> Self {
        AggLossStruct::Loss(vec![(f64::NAN, (usize::MAX, usize::MAX)); m])
    }

    pub fn search_topprop_n(&self, m_top: usize) -> (Vec<bool>, f64) {
        match self {
            AggLossStruct::Loss(loss) => {
                let m = loss.len();

                //log::debug!("loss {:?}", loss);
                if loss.iter().any(|(x, _)| x.is_nan()) {
                    log::debug!("loss any nan {:?}", loss.iter().any(|(x, _)| x.is_nan()));
                    log::debug!(
                        "loss nan {:?}",
                        loss.iter().filter(|(x, _)| x.is_nan()).map(|x| *x).count()
                    );
                }

                let mut loss_sort = loss.clone();

                loss_sort.iter_mut().for_each(|x| {
                    if x.0.is_nan() {
                        x.0 = f64::MAX;
                    }
                });

                vec::sort_float(&mut loss_sort);

                let top_n = m_top.min(m);
                if top_n >= loss_sort.len() {
                    return (vec![true; m], f64::MAX);
                }

                // the smaller, the better
                let loss_topprop = loss_sort[top_n - 1].0;
                let loss_topprop_top_outside = loss_sort[top_n].0;

                let mut use_snvs = vec![false; m];

                for (mi, loss) in loss.iter().enumerate() {
                    if loss.0 <= loss_topprop {
                        //if (*loss <= loss_topprop) && !(skip_snv.contains(&mi)) {
                        use_snvs[mi] = true;
                    }
                }

                log::debug!(
                    "#SNVs in use_snvs by loss {} in {}",
                    use_snvs.iter().filter(|b| **b).count(),
                    m
                );

                return (use_snvs, loss_topprop_top_outside);
            }
        }
    }

    pub fn string_write(
        &self,
        agg_ids: &Vec<AggId>,
        agg_to_m: &Vec<Vec<usize>>,
        snvs: &Snvs,
        genot: &Genot,
    ) -> String {
        match self {
            AggLossStruct::Loss(x) => {
                let snvs_ids = snvs.snv_ids();

                // TODO: add filter, maf
                let str_header =
                    "index\taggid\tchrom\tpos_start\tpos_end\tgroup_pos_start\tgroup_pos_end\tloss\tfreq\tnum_snvs\tsnvs"
                        .to_owned();
                let strings = x
                    .iter()
                    .zip(agg_to_m.iter())
                    .enumerate()
                    .map(|(i, ((loss, (group_start, group_end)), agg_to_m_i))| {
                        let agg_id = &agg_ids[i];

                        let (group_start, group_end, gfreq, num_group, snvs_group) =
                            match (*group_start, *group_end) {
                                (usize::MAX, usize::MAX) => (
                                    "NaN".to_string(),
                                    "NaN".to_string(),
                                    f64::NAN,
                                    0,
                                    "NaN".to_string(),
                                ),
                                (group_start_idx, group_end_idx) => {
                                    // group_start_idx, group_end_idx; index in agg_to_m_i

                                    let digit = 0;

                                    let mut gbits_or = CBits::new_zeros(genot.n());
                                    for j in group_start_idx..group_end_idx + 1 {
                                        gbits_or.or_bitwise(
                                            genot.genot_inner().as_cbits_ref(agg_to_m_i[j], digit),
                                        );
                                    }

                                    let gfreq = gbits_or.maf_group();

                                    let num_group = group_end_idx - group_start_idx + 1;
                                    let snvs_group = (group_start_idx..group_end_idx + 1)
                                        .map(|i| snvs_ids[agg_to_m_i[i]].idma().to_string())
                                        .collect::<Vec<String>>()
                                        .join(",")
                                        .to_string();

                                    let group_start_mi = agg_to_m_i[group_start_idx];
                                    let group_end_mi = agg_to_m_i[group_end_idx];

                                    (
                                        snvs_ids[group_start_mi].pos().to_string(),
                                        snvs_ids[group_end_mi].pos().to_string(),
                                        gfreq,
                                        num_group,
                                        snvs_group,
                                    )
                                }
                            };

                        format!(
                            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                            i,
                            agg_id.id(),
                            agg_id.chrom(),
                            agg_id.pos_start(),
                            agg_id.pos_end(),
                            group_start,
                            group_end,
                            loss,
                            gfreq,
                            num_group,
                            snvs_group,
                        )
                    })
                    .collect::<Vec<String>>();
                format!("{}\n{}", str_header, strings.join("\n"))
            }
        }
    }

    pub fn write_writer(
        &self,
        agg_ids: &Vec<AggId>,
        agg_to_m: &Vec<Vec<usize>>,
        snvs: &Snvs,
        genot: &Genot,
        writer: &mut std::io::BufWriter<File>,
    ) {
        let str_out = self.string_write(agg_ids, agg_to_m, snvs, genot);
        writer.write_all(&str_out.as_bytes()).unwrap();
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    //#[test]
    //fn test_error_rate_tp_fn() {
    //    let prev = CBits::new(&vec![false, false, false, true, true]);
    //    let gbits = CBits::new(&vec![true, true, false, false, true]);
    //    let phe = Phe::new(&vec![true, false, true, true, true]);

    //    //let tp_fn = (3, 5);
    //    let tp_fn_new = error_rate_tp_fn_add(&prev.as_cbits_ref_b(), &gbits.as_cbits_ref_b(), &phe);
    //    assert_eq!(tp_fn_new, (1, 1));
    //}
}
