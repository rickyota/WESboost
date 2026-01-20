//! Application of **Genoboost** for research.
//! Input plink file to run Genoboost.
//!
//!
//!
//! Input is one of following
//! 1. plink2 + fin_phe (phe and cov)
//! 2. plink2 (phe) + fin_phe (cov) : phe_name is in psam
//! 3. plink2 (cov and phe in .psam) : phe_buf is None
//! 4. plink1 + fin_phe (cov and phe)
//! 5. plink1 (phe) + fin_phe (cov) : phe_name is None
//!
//!
//!
//!
//!
// TODO: assume FID==IID for plink1 or plink2
// TODO: ensure the same para when resuming
// TODO: (optional) write down extract snvs from top_snvs
// samples indicate both fid and iid
// trimming sample weights: only use samples with large weights on choosing SNVs: Friedman, J., Hastie, T. and Tibshirani, R. (2000) ‘Additive logistic regression: a statistical view of boosting (With discussion and a rejoinder by the authors)’, Annals of statistics, 28(2), pp. 337–407. doi:10.1214/aos/1016218223.
//

use clap::{ArgGroup, Args, Parser, Subcommand};
use indoc::indoc;
//use crate::boosting::{BoostMethod, BoostParam, IterationNumber};
use boosting::{
    self, AggParam, BoostParamCommon, BoostParamLrs, DoutFile, DoutScoreFile, WgtDoutOrFile,
};
use genetics::{dataset_file::DatasetFile, GenotFormatArg};
use std::{path::PathBuf, time::Instant};

#[derive(Debug, Parser)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
    // globa=true makes you able to `-- train --verbose`
    #[arg(long, global = true, help = "Number of threads")]
    threads: Option<usize>,
    #[arg(long, global = true, help = "Verbose")]
    verbose: bool,
    //#[arg(long, global = true, help = "Memory [GB]")]
    #[arg(long, global = true, help = "Memory [MB]")]
    memory: Option<usize>,
}

#[derive(Debug, Subcommand)]
enum Commands {
    #[command(about = "train")]
    Train(TrainArgs),
}

//#[command(group(ArgGroup::new("iter_or_snv").required(true).args(["iter_snv", "iter"])))]

#[derive(Debug, Args)]
//#[command(group(ArgGroup::new("iter_or_snv").args(["iter_snv", "iter"])))]
//#[command(group(ArgGroup::new("ld").args(["ldr2", "ldradius"])))]
struct TrainArgs {
    #[arg(long)]
    dir: String,
    #[arg(long, help = "Should be major-ref.")]
    file_genot: Vec<String>,
    //file_genot: String,
    #[arg(long, value_enum)]
    genot_format: GenotFormatArg,
    //#[arg(long)]
    //boost_type: String,
    #[arg(long)]
    file_sample: Option<String>,
    #[arg(long)]
    file_sample_val: Option<String>,
    // option for covs and phes are in .psam
    #[arg(long)]
    file_phe: Option<String>,
    // option for phe in plink1
    #[arg(long)]
    phe: Option<String>,
    // parse later
    #[arg(long)]
    cov: Option<String>,
    // for rare var
    //#[arg(long)]
    //file_score_start: Option<String>,
    #[arg(long)]
    agg_top: Option<usize>,
    //#[arg(long)]
    //file_cov: Option<String>,
    #[arg(long)]
    covway: Option<String>,
    //#[arg(long)]
    //interactionway: Option<String>,
    #[arg(long, help = "1-column; 'id' or 'id:REF>ALT' ")]
    file_snv: Option<String>,
    #[arg(long)]
    file_snv_funct: Option<String>,
    #[arg(long)]
    file_snv_agg: Option<String>,
    // TODO: multiple loss_func
    #[arg(long, help = "logisitc, pval, error")]
    loss_func: String,
    #[arg(long, help = "write loss for all range in agg")]
    write_agg_loss: bool,
    #[arg(long)]
    group_freq_thre: Option<f64>,
    #[arg(long)]
    skip_neg: bool,
    //loss_func: Option<String>,
    //#[arg(long)]
    //iter_snv: Option<usize>,
    //#[arg(long)]
    //file_group_snv: Option<String>,
    //#[arg(long)]
    //iter: Option<usize>,
    //#[arg(long, value_parser, num_args = 1.., value_delimiter = ' ')]
    //learning_rates: Option<Vec<f64>>,
    //#[arg(long)]
    //batch: Option<String>,
    //#[arg(long)]
    //batchinteraction: Option<String>,
    //#[arg(long)]
    //clip_sample_weight: Option<String>,
    //#[arg(long)]
    //clip_sample_wls_weight: Option<String>,
    //#[arg(long)]
    //eps: Option<String>,
    //#[arg(long)]
    //effeps: Option<String>,
    //#[arg(long)]
    //prior: Option<String>,
    //#[arg(long)]
    //h2snv: Option<f64>,
    // TODO: allow hg19, hg38 in genoboost.rs
    // default value depnds on hg19 or hg38
    // format 6:28000000-34000000
    //#[arg(long)]
    //mhc_region: Option<String>,
    //#[arg(long)]
    //mafthrecommon: Option<f64>,
    //#[arg(long)]
    //acc_metric: Option<String>,
    //#[arg(long)]
    //use_adjloss: bool,
    //#[arg(long)]
    //use_const_for_loss: bool,
    //  not good for multi-allelic
    //  plink2 should be major-ref
    //#[arg(
    //    long,
    //    help = "Set major allele in training dataset as a2 allele. Otherwise, set ref allele as a2 allele."
    //)]
    //major_a2_train: bool,
    //#[arg(long, help = "LD radius for interactionLD")]
    //ldradius: Option<usize>,
    //#[arg(long, help = "LD r2 for interactionLD")]
    //ldr2: Option<f64>,
    //#[arg(long)]
    //resume: bool,
    //#[arg(long)]
    //is_initial_only: bool,
    //#[arg(long)]
    //write_loss: bool,
    //#[arg(long)]
    //file_initial_snvs: Option<String>,
}

fn main() {
    #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
    {
        // or use _enabled!
        if is_x86_feature_detected!("avx2") {
            log::info!("Able to use SIMD.")
        } else {
            log::info!("Not able to use SIMD since avx2 is not detected.")
        }
    }
    #[cfg(not(any(target_arch = "x86", target_arch = "x86_64")))]
    {
        log::info!("Not able to use SIMD since arch is not x86 or x86_64.")
    }

    let start = Instant::now();

    let cli = Cli::parse();
    log::info!("cli: {:?}", cli);

    if cli.verbose {
        std::env::set_var("RUST_LOG", "debug");
    } else {
        std::env::set_var("RUST_LOG", "info");
    }
    env_logger::init();

    if let Some(n_threads) = cli.threads {
        rayon::ThreadPoolBuilder::new()
            .num_threads(n_threads)
            .build_global()
            .unwrap();
    };
    // otherwise, use default thread number
    log::debug!("num_thread set: {}", rayon::current_num_threads());

    let mem = cli.memory.map(|x| x * 1024 * 1024);
    //let mem = cli.memory.map(|x| x * 1024 * 1024 * 1024);
    log::debug!("Memory : {:?} Byte", mem);

    match cli.command {
        Commands::Train(args) => {
            let dout = DoutFile::new(PathBuf::from(args.dir));

            let fins_genot = args.file_genot.iter().map(|x| PathBuf::from(x)).collect();
            let fins_genot = args.genot_format.to_genot_file(fins_genot);
            let fin_phe = args.file_phe.map(|x| PathBuf::from(x));
            let fin_sample = args.file_sample.map(|x| PathBuf::from(x));
            let fin_sample_val = args.file_sample_val.map(|x| PathBuf::from(x));
            let fin_snv = args.file_snv.map(|x| PathBuf::from(x));

            let mut dfile = DatasetFile::new(
                fins_genot,
                fin_phe,
                args.phe,
                args.cov,
                fin_snv,
                fin_sample,
                fin_sample_val,
            );
            //dfile.update_file_snv_funct(args.file_snv_funct.map(|x| PathBuf::from(x)));
            // TODO: move to boosting_res
            dfile.update_file_agg_snv(args.file_snv_agg.map(|x| PathBuf::from(x)));
            dfile.reads();
            let dfile = dfile;
            dfile.check_valid_fin();

            // TODO
            //let acc_metric = Some("cov-adjusted-pseudo-r2".to_string());

            let agg_param = AggParam::default()
                .set_loss_func(&args.loss_func)
                .set_group_freq_thre(args.group_freq_thre)
                .set_skip_neg(args.skip_neg);

            //let boost_param_common = BoostParamCommon::default()
            //    .set_loss_func("logistic")
            //    .set_is_monitor(is_monitor)
            //    .set_sample_weight_clip(args.clip_sample_weight.as_deref())
            //    .set_sample_weight_wls_clip(args.clip_sample_wls_weight.as_deref())
            //    .set_eps(args.eps.as_deref())
            //    .set_cov_way(args.covway.as_deref())
            //    .set_interaction_way(args.interactionway.as_deref())
            //    .set_batch_way(args.batch.as_deref())
            //    .set_batch_interaction_way(args.batchinteraction.as_deref())
            //    .set_eff_eps(args.effeps.as_deref())
            //    .set_prior(args.prior.as_deref(), args.h2snv)
            //    .set_mhc_region(args.mhc_region.as_deref())
            //    .set_maf_threshold_logit_common(args.mafthrecommon)
            //    .set_ld_criteria(args.ldr2, args.ldradius)
            //    .set_acc_metric(acc_metric.as_deref());

            //let boost_param_common = if args.iter.is_some() || args.iter_snv.is_some() {
            //    boost_param_common.set_iteration(args.iter, args.iter_snv)
            //} else {
            //    if !is_monitor {
            //        panic!("Either iter or iter_snv must be indicated when not monitoring.");
            //    }
            //    boost_param_common
            //};

            //let boost_params = BoostParamLrs::default()
            //    .set_boost_type(&*args.boost_type)
            //    .set_learning_rates(learning_rates)
            //    .set_boost_param_common(boost_param_common);

            //agg_param.check();
            log::debug!("boost_params {:?}", agg_param);

            log::info!("dout {:?}", dout);
            log::info!("agg_params {:?}", agg_param);

            // Genoboost
            //crate::boosting::run_boosting(
            crate::boosting::run_aggregation(
                &dout,
                &dfile,
                &agg_param,
                args.agg_top,
                //args.resume,
                args.write_agg_loss,
                //args.is_initial_only,
                //None, //prune_snv,
                //args.major_a2_train,
                mem,
            );
        }
    }

    let end = start.elapsed();
    log::info!("It took {} seconds.", end.as_secs());
    log::info!("Done!!");
}
