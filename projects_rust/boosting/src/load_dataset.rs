pub use crate::boosting_param::BoostParamLrs;
pub use crate::dout_file::DoutFile;
use genetics::dataset_file::DatasetFile;
use genetics::{Dataset, FillMissing, FillMissingGroup};

pub fn load_dataset(
    dfile: &DatasetFile,
    fill_missing: Option<FillMissing>,
    fill_missing_group: Option<FillMissingGroup>,
    make_major_a2_train: bool,
    mem: Option<usize>,
) -> (Dataset, Option<Dataset>) {
    log::debug!("Load dataset from DatasetFile for training dataset.");
    let dataset = Dataset::new_datasetfile_training(
        dfile,
        false,
        None,
        fill_missing,
        fill_missing_group,
        make_major_a2_train,
        None,
        mem,
    );

    log::debug!("Load dataset from DatasetFile for validation dataset.");
    let dataset_val = if dfile.sample_val_buf().is_some() {
        let mem_val = mem.map(|x| x - dataset.genot().byte_self());
        let dataset_val = Dataset::new_datasetfile_training(
            dfile,
            true,
            None,
            fill_missing,
            fill_missing_group,
            make_major_a2_train,
            Some(dataset.snvs()),
            mem_val,
        );
        Some(dataset_val)
    } else {
        None
    };
    (dataset, dataset_val)
}

pub fn load_dataset_boosting(
    _dout: &DoutFile, // for prune
    dfile: &DatasetFile,
    boost_params: &BoostParamLrs,
    prune_snv: Option<f64>,
    make_major_a2_train: bool,
    mem: Option<usize>,
) -> (Dataset, Option<Dataset>) {
    let boost_param = boost_params.param_lr_none();
    let fill_missing_mode = boost_param.boost_type().fill_missing();
    let fill_missing = if fill_missing_mode {
        Some(FillMissing::Mode)
    } else {
        None
    };

    let fill_missing_group = Some(FillMissingGroup::Ref);

    // create dataset
    // extract snvs by loss function
    if let Some(_prop_prune_snv) = prune_snv {
        unimplemented!();
        // prune_snv_boosting()
    } else {
        load_dataset(
            dfile,
            fill_missing,
            fill_missing_group,
            make_major_a2_train,
            mem,
        )
        //let dataset = Dataset::new_datasetfile_training(
        //    dfile,
        //    false,
        //    None,
        //    fill_missing,
        //    make_major_a2_train,
        //    None,
        //    mem,
        //);

        //let dataset_val = if dfile.sample_val_buf().is_some() {
        //    // bug: sample_buf could exists even though sample_val does not exist.
        //    // -> solved
        //    // let dataset_val = if dfile.fin_sample_val().is_some() {
        //    // TODO: if fin for training and validation are different file,
        //    // need to align and sort snv

        //    let mem_val = mem.map(|x| x - dataset.genot().byte_self());
        //    let dataset_val = Dataset::new_datasetfile_training(
        //        dfile,
        //        true,
        //        None,
        //        fill_missing,
        //        make_major_a2_train,
        //        Some(dataset.snvs()),
        //        mem_val,
        //    );
        //    Some(dataset_val)
        //} else {
        //    None
        //};
        //(dataset, dataset_val)
    }
}

// To prune snv,
// 1. use --loss-for-prune to get loss and prune top n snvs on bash, or
// 2. [TODO] implement above in rust
//
// // old
// fn prune_snv_boosting() -> (Dataset, Option<Dataset>) {
//     log::info!("Prune SNVs by decreasing loss: {}", prop_prune_snv);
//     let start = Instant::now();

//     sample_val_buf.unwrap_or_else(|| panic!("Not Implemented"));
//     //.expect("Not Implemented");
//     //fin_sample_val.expect("Not Implemented");

//     // TODO: better
//     // to get m
//     let m;
//     {
//         //let m_in: usize = plink::compute_num_snv(fin, gfmt).unwrap();
//         let snvs_in = io_genot::load_snvs(fin, gfmt);
//         //(m, _) = snv::make_use_snvs(fin_snv, &snvs_in);
//         (m, _) = snv::make_use_snvs_buf(snv_buf, &snvs_in);
//     }

//     // TODO: depends on available memory
//     //let m_range = 200_000usize;
//     //let m_range = 400_000usize;
//     let m_range = 1_700_000usize;
//     //let m_range = 2_000_000usize;
//     //let m_range = 800_000usize;

//     let mut losss = Vec::new();
//     for snv_i_start in (0..m).step_by(m_range) {
//         let mut filt_snv = vec![false; m];
//         let snv_i_end = (snv_i_start + m_range).min(m);
//         let m_ = snv_i_end - snv_i_start;
//         filt_snv[snv_i_start..snv_i_end].fill(true);
//         assert_eq!(vec::count_true(&filt_snv), m_);

//         let dataset = Dataset::new(
//             fin,
//             gfmt,
//             phe_buf.as_deref(),
//             phe_name,
//             cov_name,
//             snv_buf.as_deref(),
//             sample_buf.as_deref(),
//             Some(&filt_snv),
//             fill_missing,
//             make_major_a2_train,
//             None,
//             mem,
//         );

//         let n = dataset.samples().samples_n();

//         let mut scores: Vec<f64> = vec![0.0; n];

//         let mut sample_weight = SampleWeight::new(
//             n,
//             boost_param.boost_type(),
//             boost_param.loss_func(),
//             boost_param.sample_weight_clip(),
//             boost_param.sample_weight_wls_clip(),
//         );
//         sample_weight.renew_sample_weight(&scores, dataset.samples().phe_unwrap());

//         let mut wgts = WgtBoosts::new(boost_param.boost_type());
//         let _ = boosting_train::boosting_covs(
//             &mut wgts,
//             &dataset,
//             &mut scores,
//             &mut sample_weight,
//             0,
//             None,
//             &mut [0.0f64; 0],
//         );

//         let mut loss = LossStruct::new(boost_param.boost_type(), m_);

//         boosting_train::loss::calculate_loss_gt(
//             &mut loss,
//             &dataset,
//             //dataset.genot(),
//             &sample_weight,
//             //dataset.samples().phe_unwrap(),
//             &boost_param,
//             &HashSet::<usize>::new(),
//             //use_adjloss,
//         );

//         losss.push(loss);
//     }

//     let mut v = Vec::new();
//     for loss_ in losss {
//         //let LossStruct::ModelFree(v_in, _) = loss_;
//         let v_in = loss_.inner();
//         v.extend(v_in);
//         //v.extend(loss_.inner());
//     }

//     let loss = LossStruct::new_vec(boost_param.boost_type(), v);
//     let mut writer_loss = wgt_boost::io::bufwriter_floss(dout, 0);
//     let snvs = Snvs::new_plink(fin, gfmt);
//     loss.write_writer(&mut writer_loss, &snvs);

//     let use_snvs_loss = loss.search_topprop(prop_prune_snv);

//     log::debug!("created use_snvs_loss {}", use_snvs_loss.len());

//     let dataset = Dataset::new_boost_training(
//         fin,
//         gfmt,
//         phe_buf.as_deref(),
//         phe_name,
//         cov_name,
//         snv_buf.as_deref(),
//         sample_buf.as_deref(),
//         fill_missing,
//         Some(&use_snvs_loss),
//         None,
//         make_major_a2_train,
//         mem,
//     );

//     log::debug!(
//         "It took {} seconds to prune Dataset.",
//         start.elapsed().as_secs()
//     );

//     //dataset_ext = dataset.extract_snvs(&use_snvs_loss);
//     (dataset, None)
// }
