//! Application of **Genoboost**.
//! Input plink file to run Genoboost.
//

mod aggregation_param;
mod aggregation_train;
mod boosting_param;
mod boosting_score;
pub mod boosting_train; // pub for bench
mod dout_file;
mod load_dataset;
mod score;
mod wgt_boost;
mod wgt_boosts;
mod wgt_io;

//use mysmartcore::linalg::basic::arrays::Array2;
//// mysmartcore v0.3.1
//// https://smartcorelib.org/user_guide/model_selection.html
//use mysmartcore::linalg::basic::matrix::DenseMatrix;
//use mysmartcore::model_selection::KFold;
//use mysmartcore::model_selection::*;

use std::time::Instant;

use aggregation_param::AggParamCommonTrait;
pub use aggregation_param::{AggLossFunc, AggParam};
pub use boosting_param::{
    BoostParam, BoostParamCommon, BoostParamCommonTrait, BoostParamLrs, BoostParamsTypes,
    BoostType, EffEps, Eps, IterationNumber,
};
use boosting_train::sample_weight::SampleWeight;
use boosting_train::ContingencyTable;
use dout_file::DoutParaFile;
pub use dout_file::{DoutFile, DoutScoreFile, WgtDoutOrFile};
use genetics::dataset::FillMissingScore;
use genetics::dataset_file::DatasetFile;
use genetics::genot_io;
use genetics::sample;
use genetics::textfile;
use genetics::FillMissing;
use genetics::FillMissingGroup;
use wgt_boost::WgtBoost;

// for prune
//use crate::boosting_train::sample_weight::SampleWeight;
//use boosting_train::loss::LossStruct;
//use genetics::snv;
//use genetics::{vec, Snvs};

#[macro_use]
extern crate assert_float_eq;

/// Run GenoBoost. Input plink file.
pub fn run_boosting(
    dout: &DoutFile,
    dfile: &DatasetFile,
    boost_params: &BoostParamLrs,
    is_resume: bool,
    is_write_loss: bool,
    is_initial_only: bool,
    prune_snv: Option<f64>,
    loss_for_prune: bool,
    make_major_a2_train: bool,
    mem: Option<usize>,
) {
    // check fwgt does not exist.
    if !is_resume {
        for learning_rate in boost_params.learning_rates().iter() {
            let dout_para = dout.dout_para_lr(*learning_rate);
            dout_para.check_file_wgt_not_exist();
        }
    }

    // TODO: if all wgt exceeds #SNVs, then exit here before loading genotype.

    dfile.check_valid_fin();

    let start_time = Instant::now();

    let (dataset, dataset_val) = load_dataset::load_dataset_boosting(
        dout,
        dfile,
        boost_params,
        prune_snv,
        make_major_a2_train,
        mem,
    );

    log::info!("Created dataset: {} sec", start_time.elapsed().as_secs());

    dout.create_dout();

    if loss_for_prune {
        log::info!("Loss for prune");
        let dout_para = dout.dout_loss_for_prune();
        // no use
        let mut writer = dout_para.bufwriter_fwgt_append();
        // lr=1.0
        let boost_param_lr = boost_params.param_lr_none();
        let _ = boosting_train::boosting_batch(
            &mut writer,
            &boost_param_lr,
            &dataset,
            dataset_val.as_ref(),
            is_resume,
            is_write_loss,
            Some(&dout_para),
            //fscore_start,
            dfile.score_start_buf(),
            true,
        );
        return;
    }

    let mut nsvs_acc_use: Vec<Option<(usize, f64)>> = vec![];
    for boost_param_lr in boost_params.clone().into_iter() {
        let learning_rate = boost_param_lr.learning_rate();
        log::debug!("learning rate: {:?}", learning_rate);
        let dout_para = dout.dout_para_lr(learning_rate);
        dout_para.create_dout_para();
        // TODO: writer should be inside of boosting()??
        let mut writer = dout_para.bufwriter_fwgt_append();
        let nsnv_acc_use = if boost_param_lr.batch_way().is_none() {
            log::info!("Run boosting");
            boosting_train::boosting(
                &mut writer,
                &boost_param_lr,
                &dataset,
                is_resume,
                is_write_loss,
                Some(&dout_para),
            )
        } else {
            let finitial_snv_interaction = dfile.fin_initial_snv();
            if boost_param_lr.boost_type().is_interaction() {
                log::info!("Run boosting");
                boosting_train::boosting_batch_interaction(
                    &mut writer,
                    &boost_param_lr,
                    &dataset,
                    dataset_val.as_ref(),
                    is_resume,
                    is_write_loss,
                    is_initial_only,
                    Some(&dout_para),
                    finitial_snv_interaction,
                )
            } else {
                //let fscore_start = dfile.fin_score_start();
                log::info!("Run boosting");
                boosting_train::boosting_batch(
                    &mut writer,
                    &boost_param_lr,
                    &dataset,
                    dataset_val.as_ref(),
                    is_resume,
                    is_write_loss,
                    Some(&dout_para),
                    //fscore_start,
                    dfile.score_start_buf(),
                    false,
                )
            }
        };
        nsvs_acc_use.push(nsnv_acc_use);
    }

    if nsvs_acc_use.iter().any(|&x| x.is_none()) {
        log::info!("Do not create the most accurate weight file across learning rates since validation file is not provided or training was suspended in some learning rates.")
    } else {
        let learning_rates = boost_params.learning_rates();
        let (nsnv_acc_max, lr) = nsvs_acc_use
            .iter()
            .map(|x| x.unwrap())
            .zip(learning_rates.iter())
            .max_by(|((_, acc_a), _), ((_, acc_b), _)| acc_a.partial_cmp(acc_b).unwrap())
            .map(|((nsnv, _), lr)| (nsnv, lr))
            .expect("Cannot determine the highest accuracy.");
        log::info!(
            "The best parameter is #SNVs={:?}, learning rate={:?}.",
            nsnv_acc_max,
            lr
        );

        let fpara_best = dout.get_fname_wgt_best_para();
        let para_best_string = wgt_io::create_best_para_buf(nsnv_acc_max, *lr);
        wgt_io::write_file_string(&fpara_best, para_best_string);
        log::info!("Write wgt para to {:?}", fpara_best);

        let wgt_best_string = wgt_io::create_best_wgt_buf(dout, nsnv_acc_max, *lr, None);
        let fwgt_best = dout.get_fname_wgt_best();
        wgt_io::write_file_string(&fwgt_best, wgt_best_string);
        log::info!("Write wgt to {:?}", fwgt_best);
    }
}

pub fn run_boosting_integrate_cv(
    dout: &DoutFile,
    dfile: &mut DatasetFile,
    //boost_method: BoostMethod,
    boost_params_types: &BoostParamsTypes,
    is_resume: bool,
    is_write_loss: bool,
    prune_snv: Option<f64>,
    make_major_a2_train: bool,
    mem: Option<usize>,
    cross_vali: Option<usize>,
    seed: Option<u64>,
) {
    match cross_vali {
        Some(cvn) => {
            // FIXME: write down sample id

            let n = match dfile.sample_buf() {
                Some(sample_buf) => {
                    let (_, n) = sample::make_use_samples_buf(Some(sample_buf), dfile.fins_genot());
                    n
                }
                None => genot_io::compute_num_sample(dfile.fins_genot()),
            };

            // create cv samples
            let prop_train_cvn1 = 0.8;
            let sample_idx_cvs: Vec<(Vec<usize>, Vec<usize>)> =
                genetics::sample::split_samples_random(n, cvn, Some(prop_train_cvn1), seed);

            let samples_use =
                sample::samples_id_use_samples(dfile.sample_buf(), dfile.fins_genot());

            log::debug!("created samples for cross validation.");

            // run each
            for cvi in 0..cvn {
                let (sample_idx_cvi, sample_idx_val_cvi) = &sample_idx_cvs[cvi];

                // create cv-specific dir
                let dout_cv = dout.dout_cv(cvi);
                let sample_string = sample::extract_sample_cvi(&sample_idx_cvi, &samples_use);
                let sample_val_string =
                    sample::extract_sample_cvi(&sample_idx_val_cvi, &samples_use);

                let sample_buf = sample::sample_string_to_buf(sample_string);
                let sample_val_buf = sample::sample_string_to_buf(sample_val_string);

                // renew sample buf for each cv
                dfile.update_sample_buf(sample_buf, sample_val_buf);

                log::debug!("update sample buf.");

                run_boosting_integrate(
                    &dout_cv,
                    &dfile,
                    //boost_method,
                    &boost_params_types,
                    is_resume,
                    is_write_loss,
                    prune_snv,
                    make_major_a2_train,
                    mem,
                );
            }
        }
        None => {
            // Genoboost
            run_boosting_integrate(
                &dout,
                &dfile,
                //boost_method,
                &boost_params_types,
                is_resume,
                is_write_loss,
                prune_snv,
                make_major_a2_train,
                mem,
            );
        }
    }
}

// TODO: concat with run_boosting()
/// Integrate Genoboost. Select best model from add and non-add
pub fn run_boosting_integrate(
    dout: &DoutFile,
    dfile: &DatasetFile,
    //boost_method: BoostMethod,
    boost_params_types: &BoostParamsTypes,
    is_resume: bool,
    is_write_loss: bool,
    prune_snv: Option<f64>,
    make_major_a2_train: bool,
    mem: Option<usize>,
) {
    // TODO: if all wgt exceeds #SNVs, then exit here.
    //io_genot::check_valid_fin(fin, gfmt);

    assert_eq!(
        dfile.sample_val_buf().is_some(),
        boost_params_types.is_monitor()
    );

    let mut nsnvss_acc_use: Vec<Vec<Option<(usize, f64)>>> = vec![];
    for boost_params in boost_params_types.clone().into_iter() {
        log::info!("boost type: {:?}", boost_params.boost_type());
        // check fwgt does not exist.
        if !is_resume {
            for boost_param_lr in boost_params.clone().into_iter() {
                let learning_rate = boost_param_lr.learning_rate();
                let dout_para =
                    dout.dout_para_integrate(learning_rate, boost_param_lr.boost_type());
                dout_para.check_file_wgt_not_exist();
            }
        }

        let start_time = Instant::now();
        // should be inside of boost type since prune_snv differs
        let (dataset, dataset_val) = load_dataset::load_dataset_boosting(
            dout,
            &dfile,
            &boost_params,
            prune_snv,
            make_major_a2_train,
            mem,
        );
        log::debug!("dataset_val exists: {:?}", dataset_val.is_some());
        log::info!("Created dataset: {} sec", start_time.elapsed().as_secs());

        dout.create_dout();

        let mut nsnvs_acc_use: Vec<Option<(usize, f64)>> = vec![];
        for boost_param_lr in boost_params.clone().into_iter() {
            let learning_rate = boost_param_lr.learning_rate();
            log::info!("learning rate: {:?}", learning_rate);
            let dout_para = dout.dout_para_integrate(learning_rate, boost_param_lr.boost_type());
            dout_para.create_dout_para();
            let mut writer = dout_para.bufwriter_fwgt_append();
            let nsnv_acc_use = if boost_param_lr.batch_way().is_none() {
                log::info!("Run boosting");
                boosting_train::boosting(
                    &mut writer,
                    &boost_param_lr,
                    &dataset,
                    is_resume,
                    is_write_loss,
                    Some(&dout_para),
                )
            } else {
                log::info!("Run boosting batch");
                boosting_train::boosting_batch(
                    &mut writer,
                    &boost_param_lr,
                    &dataset,
                    dataset_val.as_ref(),
                    is_resume,
                    is_write_loss,
                    Some(&dout_para),
                    None,
                    false,
                )
            };
            // TODO: should write all acc not only best nsnv
            //if let Some(nsnv_acc_use) = nsnv_acc_use {
            //    let facc = wgt_boost::io::get_file_acc(&dout_para);
            //    let nsnv_acc_string = create_acc_string(nsnv_acc_use);
            //    write_file_string(&facc, nsnv_acc_string);
            //}
            nsnvs_acc_use.push(nsnv_acc_use);
        }
        nsnvss_acc_use.push(nsnvs_acc_use);
    }

    if boost_params_types.is_monitor() {
        log::debug!("nsnvss_acc_use {:?}", nsnvss_acc_use);
        best_wgt_integrate(nsnvss_acc_use, boost_params_types, dout);
    }
}

fn best_wgt_integrate(
    nsnvss_acc_use: Vec<Vec<Option<(usize, f64)>>>,
    boost_params_types: &BoostParamsTypes,
    dout: &DoutFile,
) {
    let mut nsnv_max: Option<usize> = None;
    let mut lr_max: Option<f64> = None;
    let mut boost_type_max: Option<BoostType> = None;

    let mut acc_max: Option<f64> = None;
    for (boost_params, nsnvs_acc_use) in boost_params_types
        .clone()
        .into_iter()
        .zip(nsnvss_acc_use.iter())
    {
        for (boost_param_lr, nsnv_acc_use) in
            boost_params.clone().into_iter().zip(nsnvs_acc_use.iter())
        {
            //if nsnvs_acc_use.iter().any(|&x| x.is_none()) {
            // if any is None, do not create
            match nsnv_acc_use {
                None => {
                    log::info!("Do not create the most accurate weight file across learning rates since training was suspended in some learning rates.");
                    return;
                }
                Some(nsnv_acc_use) => {
                    let (nsnv, acc) = nsnv_acc_use;
                    if acc_max.is_none()
                        || (acc.partial_cmp(&acc_max.unwrap()).unwrap()
                            == std::cmp::Ordering::Greater)
                    {
                        acc_max = Some(*acc);
                        nsnv_max = Some(*nsnv);
                        lr_max = Some(boost_param_lr.learning_rate());
                        boost_type_max = Some(boost_param_lr.boost_type());
                    }
                }
            }
        }
    }

    if acc_max.is_none() {
        log::info!("Cannot determine the highest accuracy.");
        return;
    }
    //let acc_max = acc_max.unwrap();
    let nsnv_max = nsnv_max.unwrap();
    let lr_max = lr_max.unwrap();
    let boost_type_max = boost_type_max.unwrap();

    log::info!(
        "The best parameter is #SNVs={:?}, learning rate={:?}.",
        nsnv_max,
        lr_max
    );

    let fpara_best = dout.get_fname_wgt_best_para();
    let para_best_string = wgt_io::create_best_para_integrate_buf(nsnv_max, lr_max, boost_type_max);
    wgt_io::write_file_string(&fpara_best, para_best_string);

    let wgt_best_string = wgt_io::create_best_wgt_buf(dout, nsnv_max, lr_max, Some(boost_type_max));
    let fwgt_best = dout.get_fname_wgt_best();
    wgt_io::write_file_string(&fwgt_best, wgt_best_string);
}

pub fn run_boosting_score_cv(
    dout_score: &DoutScoreFile,
    dfile: &DatasetFile,
    is_every_para: bool,
    iterations_in: Option<&[usize]>,
    wgt_d_f: &WgtDoutOrFile,
    learning_rates: &[f64],
    use_iter: bool,
    cross_vali: Option<usize>,
    // fill_missing_in_dataset: bool,
    missing_score: Option<FillMissingScore>,
    fill_missing_group: Option<FillMissingGroup>,
    allow_nonexist_snv: bool,
    use_snv_pos: bool,
    mem: Option<usize>,
) {
    match cross_vali {
        Some(cvn) => {
            if wgt_d_f.is_file() {
                panic!("Use --dout_wgt for cross-validation");
            }

            for cvi in 0..cvn {
                let dout_wgt = wgt_d_f.dout_wgt();
                let dout_wgt_cv = dout_wgt.dout_cv(cvi);
                let wgt_d_f_para = WgtDoutOrFile::Dout(dout_wgt_cv);

                let dout_score_cv = dout_score.dir_score_cv(cvi);

                run_boosting_score(
                    &dout_score_cv,
                    dfile,
                    is_every_para,
                    iterations_in,
                    &wgt_d_f_para,
                    &learning_rates,
                    use_iter,
                    // fill_missing_in_dataset,
                    missing_score,
                    fill_missing_group,
                    allow_nonexist_snv,
                    use_snv_pos,
                    mem,
                );
            }
        }
        None => {
            run_boosting_score(
                &dout_score,
                dfile,
                is_every_para,
                iterations_in,
                wgt_d_f,
                &learning_rates,
                use_iter,
                // fill_missing_in_dataset,
                missing_score,
                fill_missing_group,
                allow_nonexist_snv,
                use_snv_pos,
                mem,
            );
        }
    }
}

/// if SNV in fin_wgt is not in fin, then
pub fn run_boosting_score(
    dout_score: &DoutScoreFile,
    dfile: &DatasetFile,
    is_every_para: bool,
    iterations_in: Option<&[usize]>,
    wgt_d_f: &WgtDoutOrFile,
    learning_rates: &[f64],
    use_iter: bool,
    // fill_missing_in_dataset: bool,
    missing_score: Option<FillMissingScore>,
    fill_missing_group: Option<FillMissingGroup>,
    allow_nonexist_snv: bool,
    use_snv_pos: bool,
    mem: Option<usize>,
) {
    if !is_every_para {
        // score of the best para

        let dout_wgt = wgt_d_f.dout_wgt();

        let file_wgt_para = dout_wgt.get_fname_wgt_best();
        let exist_fwgt = textfile::exist_file(&file_wgt_para);
        if !exist_fwgt {
            log::info!("fwgt does not exist: {:?}.", &file_wgt_para);
            return;
        }

        let dout_score_para = dout_score.dout_score_para();

        boosting_score::run_boosting_score_para_best(
            &dout_score_para,
            dfile,
            &file_wgt_para,
            // fill_missing_in_dataset,
            missing_score,
            fill_missing_group,
            allow_nonexist_snv,
            use_snv_pos,
            mem,
        );
    } else {
        match wgt_d_f {
            WgtDoutOrFile::Dout(dout_wgt) => {
                // TODO: load dataset before lrs like classic score
                for learning_rate in learning_rates.iter() {
                    log::debug!("learning rate: {:?}", learning_rate);

                    let dout_score_lr = dout_score.dout_score_para_lr(*learning_rate);
                    log::debug!("output to: {:?}", &dout_score_lr);

                    let dout_wgt_para = dout_wgt.dout_para_lr(*learning_rate);
                    let file_wgt_para = dout_wgt_para.get_file_wgt();
                    let exist_fwgt = textfile::exist_file(&file_wgt_para);
                    if !exist_fwgt {
                        log::info!("fwgt does not exist: {:?}.", &file_wgt_para);
                        continue;
                    }

                    boosting_score::run_boosting_score_para(
                        &dout_score_lr,
                        dfile,
                        iterations_in.unwrap(),
                        &file_wgt_para,
                        use_iter,
                        // fill_missing_in_dataset,
                        allow_nonexist_snv,
                        use_snv_pos,
                        missing_score,
                        // missing_to_mode,
                        // missing_to_mean,
                        fill_missing_group,
                        mem,
                    );
                }
            }
            WgtDoutOrFile::File(file_wgt) => {
                let dout_score_para = dout_score.dout_score_para();
                boosting_score::run_boosting_score_para(
                    &dout_score_para,
                    dfile,
                    iterations_in.unwrap(),
                    &file_wgt,
                    use_iter,
                    // fill_missing_in_dataset,
                    allow_nonexist_snv,
                    use_snv_pos,
                    missing_score,
                    // missing_to_mode,
                    // missing_to_mean,
                    fill_missing_group,
                    mem,
                );
            }
        };
    }
}

/*
// moved to integration test ./tests/boosting.rs
#[cfg(test)]
mod tests {
    use super::*;

    use std::env;
    //use std::path::Path;

    #[test]
    fn test_run_boosting() {
        let path = env::current_dir().unwrap();
        println!("current dir: {}", path.display());

        let fout = String::from("../../tests/boosting/result/toy1/unittest/toy1");
        //let fout = String::from("../../tests/boosting/result/toy1/tests_toy1");
        use io_rust::text;
        use std::fs;
        use std::path::Path;
        let fwgt = fout.to_owned() + ".wgt";
        if text::exist_file(&fwgt) {
            // delete file
            fs::remove_file(fwgt).unwrap_or_else(|why| {
                println!("! {:?}", why.kind());
            });
        }

        // tmp: better to put the file in this crate
        let fin = String::from("../../tests/boosting/data/toy1/toy1");
        //let fin = String::from("./data/boosting/1000g/1kg_phase1_all/1kg_phase1_all");

        //use std::path::Path;
        let fin_bim: String = path.to_str().unwrap().to_owned() + "/" + &fin + ".bim"; // + ".bim";
        println!("fin_bim {:?}", fin_bim);
        let path = Path::new(&fin_bim);
        println!("cano {:?}", path.canonicalize().unwrap());

        let t = 10;
        let boost_type = BoostType::Logit;

        run_boosting(
            &fout,
            &fin,
            t,
            &boost_type,
            None,
            None,
            None,
            None,
            false,
            None,
        );

        println!("Done!!");
    }

    #[test]
    fn test_run_boosting_withregcov() {
        let path = env::current_dir().unwrap();
        println!("current dir: {}", path.display());

        let fout = String::from("../../tests/boosting/result/toy1/unittest/toy1_regcov");
        //let fout = String::from("../../tests/boosting/result/toy1/tests_regcov_toy1");
        use io_rust::text;
        use std::fs;
        let fwgt = fout.to_owned() + ".wgt";
        if text::exist_file(&fwgt) {
            // delete file
            fs::remove_file(fwgt).unwrap_or_else(|why| {
                println!("! {:?}", why.kind());
            });
        }

        // tmp: better to put the file in this crate
        let fin = String::from("../../tests/boosting/data/toy1/toy1");
        //let fin = String::from("./data/boosting/1000g/1kg_phase1_all/1kg_phase1_all");

        let fin_wgt_cov = Some(String::from(
            "../../tests/boosting/data/toy1/toy1.regcov.wgtcov",
        ));

        let fin_cov = Some(String::from("../../tests/boosting/data/toy1/toy1.cov"));

        println!("fin_wgt_cov {:?}", fin_wgt_cov);

        let t = 10;
        let boost_type = BoostType::Logit;

        run_boosting(
            &fout,
            &fin,
            t,
            &boost_type,
            None,
            None,
            fin_cov.as_deref(),
            fin_wgt_cov.as_deref(),
            false,
            None,
        );

        println!("Done!!");
    }

    #[test]
    fn test_run_boosting_toy2_withregcov() {
        let path = env::current_dir().unwrap();
        println!("current dir: {}", path.display());

        let fout = String::from("../../tests/boosting/result/toy2/unittest/toy2_regcov");
        //let fout = String::from("../../tests/boosting/result/toy2/tests_regcov_toy2");
        use io_rust::text;
        use std::fs;
        let fwgt = fout.to_owned() + ".wgt";
        if text::exist_file(&fwgt) {
            // delete file
            fs::remove_file(fwgt).unwrap_or_else(|why| {
                println!("! {:?}", why.kind());
            });
        }

        // tmp: better to put the file in this crate
        let fin = String::from("../../tests/boosting/data/toy2/toy2.chrom%");
        //let fin = String::from("./data/boosting/1000g/1kg_phase1_all/1kg_phase1_all");

        /*
        //use std::path::Path;
        let fin_bim: String = path.to_str().unwrap().to_owned() + "/" + &fin + ".bim"; // + ".bim";
        println!("fin_bim {:?}", fin_bim);
        let path = Path::new(&fin_bim);
        println!("cano {:?}", path.canonicalize().unwrap());
        */

        let fin_wgt_cov = Some(String::from(
            "../../tests/boosting/data/toy2/toy2.regcov.wgtcov",
        ));

        let fin_cov = Some(String::from("../../tests/boosting/data/toy2/toy2.cov"));
        let fin_snv = Some(String::from("../../tests/boosting/data/toy2/toy2.snvs"));
        let fin_sample = Some(String::from(
            "../../tests/boosting/data/toy2/toy2.cv0.samples",
        ));

        println!("fin_wgt_cov {:?}", fin_wgt_cov);

        let t = 10;
        let boost_type = BoostType::Logit;

        run_boosting(
            &fout,
            &fin,
            t,
            &boost_type,
            fin_snv.as_deref(),
            fin_sample.as_deref(),
            fin_cov.as_deref(),
            fin_wgt_cov.as_deref(),
            false,
            None,
        );

        println!("Done!!");
    }

    #[test]
    fn test_run_boosting_score() {
        let path = env::current_dir().unwrap();
        println!("current dir: {}", path.display());

        let fout = String::from("../../tests/boosting/result/toy1/unittest/toy1_regcov");
        let fin_wgt = String::from("../../tests/boosting/data/toy1/toy1_for_score");
        //let fin_wgt = String::from("../../tests/boosting/result/toy1/unittest/toy1_for_score");

        // tmp: better to put the file in this crate
        let fin = String::from("../../tests/boosting/data/toy1/toy1");
        let fin_cov = Some(String::from("../../tests/boosting/data/toy1/toy1.cov"));

        let iters = [1, 2, 10];

        run_boosting_score(&fout, &fin, &iters, &fin_wgt, fin_cov.as_deref(), None);

        println!("Done!!");
    }

    #[test]
    fn test_run_boosting_score_withregcov() {
        let path = env::current_dir().unwrap();
        println!("current dir: {}", path.display());

        let fout = String::from("../../tests/boosting/result/toy1/unittest/toy1");
        let fin_wgt = String::from("../../tests/boosting/data/toy1/toy1_regcov_for_score");
        //let fin_wgt = String::from("../../tests/boosting/result/toy1/unittest/toy1_for_score");

        // tmp: better to put the file in this crate
        let fin = String::from("../../tests/boosting/data/toy1/toy1");
        //let fin = String::from("./data/boosting/1000g/1kg_phase1_all/1kg_phase1_all");

        let iters = [1, 2, 10];

        run_boosting_score(&fout, &fin, &iters, &fin_wgt, None, None);

        println!("Done!!");
    }

    #[test]
    fn test_run_boosting_score_toy2_withregcov() {
        let path = env::current_dir().unwrap();
        println!("current dir: {}", path.display());

        let fout = String::from("../../tests/boosting/result/toy2/unittest/toy2");
        //let fin_wgt = String::from("../../tests/boosting/result/toy2/unittest/toy2_regcov_for_score");
        let fin_wgt = String::from("../../tests/boosting/data/toy2/toy2_regcov_for_score");

        // tmp: better to put the file in this crate
        let fin = String::from("../../tests/boosting/data/toy2/toy2.chrom%");
        let fin_cov = Some(String::from("../../tests/boosting/data/toy2/toy2.cov"));
        let fin_sample = Some(String::from(
            "../../tests/boosting/data/toy2/toy2.cv0.samples",
        ));

        let iters = [1, 2, 10];

        run_boosting_score(
            &fout,
            &fin,
            &iters,
            &fin_wgt,
            fin_cov.as_deref(),
            fin_sample.as_deref(),
        );

        println!("Done!!");
    }
}
*/

pub fn run_aggregation(
    dout: &DoutFile,
    dfile: &DatasetFile,
    agg_params: &AggParam,
    agg_min_m: Option<usize>,
    is_write_agg_loss: bool,
    //make_major_a2_train: bool,
    mem: Option<usize>,
) {
    let start_time = Instant::now();

    // TODO: How to deal with fill_missing:?
    // TOCHECK: if not major_ref, raise error
    let (mut dataset, _) = load_dataset::load_dataset(
        dfile,
        Some(FillMissing::Ref),
        Some(FillMissingGroup::Ref),
        false,
        mem,
    );

    log::info!("Created dataset: {} sec", start_time.elapsed().as_secs());

    // modify agg_to_m in dataset.snvs().agg_to_m()
    // if skip-neg is Some()

    if agg_params.skip_neg() {
        log::info!("Skip negative SNVs in aggregation.");
        dataset.skip_neg_agg();
    } else {
        log::info!("Do not skip negative SNVs in aggregation.");
    };

    log::info!("Start aggregation: {} sec", start_time.elapsed().as_secs());

    aggregation_train::aggregation(
        &agg_params,
        &dataset,
        //dataset_val.as_ref(),
        //is_initial_only,
        &dout.dout_agg_para(), //Some(&dout_para),
        agg_min_m,
        //finitial_snv_interaction,
        is_write_agg_loss,
    );

    log::info!("Done aggregation: {} sec", start_time.elapsed().as_secs());
}
