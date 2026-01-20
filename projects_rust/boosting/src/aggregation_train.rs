pub mod loss;

use rayon::prelude::*;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::time::Instant;

use crate::aggregation_param::AggParamCommonTrait;
use crate::dout_file::{DoutFile, DoutParaFile};
use crate::SampleWeight;
use crate::{AggLossFunc, AggParam};
use cmatrix::prelude::*;
use genetics::genot::prelude::*;
use genetics::samples::prelude::*;
use genetics::{sample::BasePhe, BaseGenot, Dataset, Snvs};
use genetics::{vec, SampleScore};
use loss::AggLossStruct;

pub fn loss_whole_region(
    agg_params: &AggParam,
    genot: &Genot,
    phe: &Phe,
    region_snv_indexs: Vec<usize>,
    // snvs: &Snvs,
    // mut writer_agg_loss: Option<&mut BufWriter<File>>,
) -> Option<(f64, (usize, usize))> {
    if region_snv_indexs.len() == 0 {
        return None;
        // return (f64::NAN, (usize::MAX, usize::MAX));
    }

    let genot_cmatrix = genot.genot_inner();
    let case_n = phe.count();

    // digit_i=0 : dominant
    // missing is assumed to be ref
    let digit = 0usize;

    // >1.0
    // let mut error_rate_min = 2.0;
    // let mut loss_min = f64::MAX;
    // let mut loss_min_region: Option<(usize, usize)> = None;

    // let snv_ids = snvs.snv_ids();

    // index: 0, 1, 2, 3, 4
    // pos  :  1, 1, 2, 3, 3
    //            0-1
    //            0  -  2
    //            0      -      4
    //                    2
    //                    2  -  4
    //                         3-4
    //
    // for start,
    //      skip if pos[i]==pos[i-1]
    // for end,
    //      skip if pos[j]==pos[j+1]

    // the smaller, the better
    let fn_loss = match agg_params.loss_func() {
        AggLossFunc::ErrorRate => loss::error_rate,
        AggLossFunc::Logistic => loss::logistic_loss,
        AggLossFunc::Pval => loss::logistic_log10_pval,
    };

    // TODO: loss_0 = NaN for logistic loss
    // let loss_0 = match agg_params.loss_func() {
    //     AggLossFunc::ErrorRate => 0.0f64,
    //     AggLossFunc::Logistic => loss::logistic_loss_0(&phe, case_n),
    // };

    let i = 0;

    // log::debug!("i: {}", i);
    let mut gbits_or = genot_cmatrix
        .as_cbits_ref(region_snv_indexs[i], digit)
        .clone_to_cbits();

    for j in i..region_snv_indexs.len() {
        let gbits_j = genot_cmatrix.as_cbits_ref(region_snv_indexs[j], digit);

        gbits_or.or_bitwise(gbits_j.as_cbits_ref_b());
    }

    let loss_ij = fn_loss(&gbits_or.as_cbits_ref_b(), &phe, case_n);
    let loss_min = loss_ij;
    let loss_min_region = Some((0, region_snv_indexs.len() - 1));
    // bug
    // let loss_min_region = Some((0, region_snv_indexs.len()));

    match loss_min_region {
        Some(x) => Some((loss_min, x)),
        None => None,
    }
}

// TODO: simd
// best region is [i, j]
pub fn loss_min_region(
    agg_params: &AggParam,
    genot: &Genot,
    phe: &Phe,
    region_snv_indexs: Vec<usize>,
    snvs: &Snvs,
    mut writer_agg_loss: Option<&mut BufWriter<File>>,
) -> Option<(f64, (usize, usize))> {
    if region_snv_indexs.len() == 0 {
        return None;
        // return (f64::NAN, (usize::MAX, usize::MAX));
    }

    let genot_cmatrix = genot.genot_inner();
    let case_n = phe.count();

    // digit_i=0 : dominant
    // missing is assumed to be ref
    let digit = 0usize;

    // >1.0
    // let mut error_rate_min = 2.0;
    let mut loss_min = f64::MAX;
    let mut loss_min_region: Option<(usize, usize)> = None;

    let snv_ids = snvs.snv_ids();

    // index: 0, 1, 2, 3, 4
    // pos  :  1, 1, 2, 3, 3
    //            0-1
    //            0  -  2
    //            0      -      4
    //                    2
    //                    2  -  4
    //                         3-4
    //
    // for start,
    //      skip if pos[i]==pos[i-1]
    // for end,
    //      skip if pos[j]==pos[j+1]

    // the smaller, the better
    let fn_loss = match agg_params.loss_func() {
        AggLossFunc::ErrorRate => loss::error_rate,
        AggLossFunc::Logistic => loss::logistic_loss,
        AggLossFunc::Pval => loss::logistic_log10_pval,
    };

    // TODO: loss_0 = NaN for logistic loss
    // let loss_0 = match agg_params.loss_func() {
    //     AggLossFunc::ErrorRate => 0.0f64,
    //     AggLossFunc::Logistic => loss::logistic_loss_0(&phe, case_n),
    // };

    for i in 0..region_snv_indexs.len() {
        // assume all chroms are same
        // for duplicated pos, start from the first one.
        if i != 0 && snv_ids[region_snv_indexs[i]].pos() == snv_ids[region_snv_indexs[i] - 1].pos()
        {
            continue;
        }

        // log::debug!("i: {}", i);
        let mut gbits_or = genot_cmatrix
            .as_cbits_ref(region_snv_indexs[i], digit)
            .clone_to_cbits();

        // TODO: i -> i+1
        for j in i..region_snv_indexs.len() {
            let gbits_j = genot_cmatrix.as_cbits_ref(region_snv_indexs[j], digit);

            gbits_or.or_bitwise(gbits_j.as_cbits_ref_b());

            // should be after updating cbits_or
            // for duplicated pos, output loss for the last snv.
            if j != region_snv_indexs.len() - 1
                && snv_ids[region_snv_indexs[j]].pos() == snv_ids[region_snv_indexs[j] + 1].pos()
            {
                continue;
            }
            // log::debug!("j: {}", j);

            // TOFIX: log10p might be inf
            let loss_ij = fn_loss(&gbits_or.as_cbits_ref_b(), &phe, case_n);
            // log::debug!("loss_ij: {}", loss_ij);
            // let error_rate = loss::error_rate(&gbits_or.as_cbits_ref_b(), &phe, case_n);

            if let Some(writer) = writer_agg_loss.as_mut() {
                let (mlog10_pval, z_score) =
                    loss::logistic_mlog10_pval_z_score(&gbits_or.as_cbits_ref_b(), &phe, case_n);
                // let pval = loss::logistic_pval(&gbits_or.as_cbits_ref_b(), &phe, case_n);
                // let pval = pval.unwrap_or(f64::NAN);

                let snv_id_i = &snv_ids[region_snv_indexs[i]];
                let snv_id_j = &snv_ids[region_snv_indexs[j]];
                let str_out = format!(
                    "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.10e}\t{:.5e}\t{:.5e}\n",
                    i,
                    j,
                    snv_id_i.pos(),
                    snv_id_j.pos(),
                    snv_id_i.ida(),
                    snv_id_j.ida(),
                    j - i + 1,
                    gbits_or.maf_group(),
                    loss_ij, // loss
                    // loss_0 - loss_ij, // relative loss
                    mlog10_pval, // -log10_pval
                    // -pval.log10() // -log10_pval
                    z_score
                );
                writer.write_all(str_out.as_bytes()).unwrap();
            }

            if loss_ij.is_nan() {
                continue;
            }
            if loss_ij < loss_min
                && gbits_or.maf_group() > agg_params.group_freq_thre().unwrap_or(0.0)
            {
                loss_min = loss_ij;
                // loss_min_region = Some((region_snv_indexs[i], region_snv_indexs[j]));
                loss_min_region = Some((i, j));
            }
        }
    }

    // if loss_min_region.is_none() {
    //     panic!("loss_min_region is None");
    // }

    // return (loss_min, loss_min_region.unwrap());
    match loss_min_region {
        Some(x) => Some((loss_min, x)),
        None => None,
    }
}

//pub fn error_rate_tps_fns(dataset: &Dataset, indexs: Vec<usize>) -> (Vec<usize>, Vec<usize>) {
//    let mut tps = vec![0usize; indexs.len()];
//    let mut fns = vec![0usize; indexs.len()];
//
//    let genot = dataset.genot();
//    let genot_cmatrix = genot.genot_inner();
//    let phe = dataset.samples().phe_unwrap();
//
//    // digit_i=0 : dominant
//    let digit_i = 0usize;
//
//    let prev_v = CBits::new(&vec![false; genot_cmatrix.col_n()]);
//    let prev = prev_v.as_cbits_ref_b();
//
//    for (ind, snv_index) in indexs.iter().enumerate() {
//        let tp_prev = tps[ind - 1];
//        let fn_prev = fns[ind - 1];
//        let gbits = genot_cmatrix.as_cbits_ref(*snv_index, digit_i);
//        let tp_fn_add = loss::error_rate_tp_fn_add(&prev, &gbits, &phe);
//        let tp = tp_prev + tp_fn_add.0;
//        let fn_ = fn_prev + tp_fn_add.1;
//        tps[ind] = tp;
//        fns[ind] = fn_;
//    }
//
//    (tps, fns)
//}

// boosting_train::loss::calc::calc_loss_logit()
// pub fn calculate_loss_logit() {
//     let (coef, _, _) = coefficient::calculate_coef_logit_eps(
//         &genot.to_genot_snv(mi),
//         sample_weight.wzs_pad().unwrap(),
//         sample_weight.wls_pad().unwrap(),
//         phe,
//         epsilons_wzs,
//         epsilons_wls,
//         boost_param.eps(),
//         //boost_param.learning_rate(),
//         boost_param.eff_eps(),
//         boost_param.boost_type(),
//         true,
//         false,
//     );
//     //*loss = calculate_loss_gt_logit_simd_sm(
//     *loss = calc_loss_logit_mi()(
//         &genot.to_genot_snv(mi),
//         &coef,
//         sample_weight.wls_pad().unwrap(),
//         sample_weight.zs_pad().unwrap(),
//         //use_adjloss,
//         loss_max_theory,
//     );
// }

// pub fn aggregation(agg_params: &AggParam, dataset: &Dataset, dout: &DoutFile) {
pub fn aggregation(
    agg_params: &AggParam,
    dataset: &Dataset,
    dout: &DoutParaFile,
    agg_min_m: Option<usize>,
    is_write_agg_loss: bool,
) {
    let start_time = Instant::now();

    // see boosting_train::create_initial_interaction()

    // snvs index in thre region
    let agg_ids = dataset.snvs().agg_ids().unwrap();
    let agg_to_m = dataset.snvs().agg_to_m().unwrap();

    // log::debug!("agg_to_m: {:?}", agg_to_m);

    let use_snvs_top_n =
        match agg_min_m {
            Some(m_top) => {
                // loss whole region
                let mut loss_whole = AggLossStruct::new(agg_to_m.len());
                loss_whole.inner_mut().par_iter_mut().enumerate().for_each(
                    |(i, error_rate_region)| {
                        let ret = loss_whole_region(
                            agg_params,
                            dataset.genot(),
                            dataset.samples().phe_unwrap(),
                            agg_to_m[i].clone(),
                        );
                        if ret.is_some() {
                            *error_rate_region = ret.as_ref().unwrap().clone();
                        }
                    },
                );

                let mut writer_whole_loss = dout.bufwriter_fgroup_whole_loss();
                loss_whole.write_writer(
                    agg_ids,
                    agg_to_m,
                    dataset.snvs(),
                    dataset.genot(),
                    &mut writer_whole_loss,
                );

                log::debug!("top_m {}", m_top);
                let (use_snvs_top_n, loss_top_n) = loss_whole.search_topprop_n(m_top);
                log::debug!("loss top_m {}", loss_top_n);

                use_snvs_top_n
            }
            None => {
                vec![true; agg_to_m.len()]
            }
        };

    log::info!("Done whole region, {} sec", start_time.elapsed().as_secs());

    // loss min range
    // rayon
    let mut loss = AggLossStruct::new(agg_to_m.len());

    loss.inner_mut().iter_mut()
        .enumerate()
        .filter(|(i, _)| use_snvs_top_n[*i])
        .collect::<Vec<_>>()
        .par_iter_mut()
        .for_each(|(i,  error_rate_region)| {
            let i=*i;
            let mut writer_agg_loss = if is_write_agg_loss {
                let mut writer = dout.bufwriter_fagg_loss_min(agg_ids[i].id());
                let str_header =
                    "index_start\tindex_end\tpos_start\tpos_end\tida_start\tida_end\tnum_snvs\tfreq\tloss\tlog10_p\tz\n"
                        .to_string();
                writer.write_all(str_header.as_bytes()).unwrap();
                Some(writer)
            } else {
                None
            };

            let ret = loss_min_region(
                agg_params,
                dataset.genot(),
                dataset.samples().phe_unwrap(),
                agg_to_m[i].clone(),
                dataset.snvs(),
                writer_agg_loss.as_mut(),
            );
            if ret.is_some() {
                **error_rate_region = ret.unwrap();
            }
        });

    // TODO: add freq for group
    // TODO: num_snvs
    // TODO: list all vars in the group

    let mut writer_loss = dout.bufwriter_fgroup_loss();
    loss.write_writer(
        agg_ids,
        agg_to_m,
        dataset.snvs(),
        dataset.genot(),
        &mut writer_loss,
    );

    log::info!("Done min region, {} sec", start_time.elapsed().as_secs());
}

#[cfg(test)]
mod tests {
    use genetics::SnvId;

    use super::*;

    #[test]
    fn test_error_rate_tp_fn() {
        let genot = Genot::new(
            5,
            8,
            &vec![
                1, 1, 1, 1, 1, 1, 1, 1, // 0
                0, 0, 1, 2, 1, 0, 0, 0, // 1
                1, 2, 0, 0, 2, 0, 0, 0, // 2
                1, 2, 0, 0, 0, 0, 1, 0, // 3
                1, 1, 1, 1, 1, 1, 1, 1, // 4
            ],
        );

        // no duplicated pos
        let snv_ids = vec![
            SnvId::new(
                "rs1".to_string(),
                "chr1".to_string(),
                "1",
                "A".to_string(),
                "C".to_string(),
            ),
            SnvId::new(
                "rs2".to_string(),
                "chr1".to_string(),
                "2",
                "A".to_string(),
                "C".to_string(),
            ),
            SnvId::new(
                "rs3".to_string(),
                "chr1".to_string(),
                "3",
                "A".to_string(),
                "C".to_string(),
            ),
            SnvId::new(
                "rs4".to_string(),
                "chr1".to_string(),
                "4",
                "A".to_string(),
                "C".to_string(),
            ),
            SnvId::new(
                "rs5".to_string(),
                "chr1".to_string(),
                "5",
                "A".to_string(),
                "C".to_string(),
            ),
        ];
        let snvs = Snvs::new_from_snv_ids(snv_ids);

        let phe = Phe::new(&vec![true, true, true, true, false, false, false, false]);

        let region_snv_indexs = vec![1, 2, 3];

        let (error_rate, error_rate_region) = loss_min_region(
            &AggParam::default().set_loss_func("error"),
            &genot,
            &phe,
            region_snv_indexs,
            &snvs,
            None,
        )
        .unwrap();

        assert_eq!(error_rate, 0.125);
        assert_eq!(error_rate_region, (0, 1));
        // assert_eq!(error_rate_region, (1, 2));
    }

    #[test]
    fn test_error_rate_tp_fn_dup_pos() {
        // index: 0, 1, 2, 3, 4
        // pos  :  1, 1, 2, 3, 3
        //            0-1
        //            0  -  2
        //            0      -      4
        //                    2
        //                    2  -  4
        //                         3-4

        let genot = Genot::new(
            5,
            8,
            &vec![
                1, 1, 1, 1, 1, 1, 1, 1, // 0
                0, 0, 1, 2, 1, 0, 0, 0, // 1
                1, 2, 0, 0, 2, 0, 0, 0, // 2
                1, 2, 0, 0, 0, 0, 1, 0, // 3
                1, 1, 1, 1, 1, 1, 1, 1, // 4
            ],
        );

        // no duplicated pos
        let snv_ids = vec![
            SnvId::new(
                "rs1".to_string(),
                "chr1".to_string(),
                "1",
                "A".to_string(),
                "C".to_string(),
            ),
            SnvId::new(
                "rs1".to_string(),
                "chr1".to_string(),
                "1",
                "A".to_string(),
                "T".to_string(),
            ),
            SnvId::new(
                "rs2".to_string(),
                "chr1".to_string(),
                "2",
                "A".to_string(),
                "C".to_string(),
            ),
            SnvId::new(
                "rs3".to_string(),
                "chr1".to_string(),
                "3",
                "A".to_string(),
                "C".to_string(),
            ),
            SnvId::new(
                "rs3".to_string(),
                "chr1".to_string(),
                "3",
                "A".to_string(),
                "T".to_string(),
            ),
        ];
        let snvs = Snvs::new_from_snv_ids(snv_ids);

        let phe = Phe::new(&vec![true, true, true, true, false, false, false, false]);

        let region_snv_indexs = vec![0, 1, 2, 3, 4];

        let (error_rate, error_rate_region) = loss_min_region(
            &AggParam::default().set_loss_func("error"),
            &genot,
            &phe,
            region_snv_indexs,
            &snvs,
            None,
        )
        .unwrap();

        // TODO: correct answer
        // assert_eq!(error_rate, 0.125);
        // assert_eq!(error_rate_region, (1, 6));
    }

    #[test]
    fn test_gfreq_thre() {
        //TODO
    }
}
