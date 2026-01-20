//! Dataset
//! It is important that all vector values can be extract without constructing. ex. mafs: Vec<f64> not Vec<Info>
//! Must use fn to access data so that it is easy to use Trait
//! For now, do not use generics, use in dataset_future.rs
//!
//!

//pub mod genot_io;
//pub mod samples;
//pub mod snvs;

use cmatrix::dense::BaseCMatrix;
use rayon::iter::{ParallelBridge, ParallelIterator};
use std::collections::HashSet;
use std::str::FromStr;
use std::time::Instant;

use crate::dataset_file::DatasetFile;
use crate::genot::prelude::*;
use crate::sample::{phe, BasePhe};
use crate::{assoc, genot_io, SnvId};
use crate::{sample, snv, vec, wgt::WgtTrait, GenotFiles, Samples, Snvs, Wgts};

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum FillMissing {
    Mode,
    Ref,
}

impl FillMissing {
    /// On generating genot, only REF should be used.
    /// For FillMissingScore::Mode and FillMissingScore::Mean, treat them in calculating score.
    pub fn from_fill_missing_score_generate_genot(
        fill_missing_score: Option<FillMissingScore>,
    ) -> Option<Self> {
        match fill_missing_score {
            Some(FillMissingScore::Ref) => Some(FillMissing::Ref),
            Some(FillMissingScore::Mode) | Some(FillMissingScore::Mean) => None,
            None => None,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum FillMissingScore {
    Mode,
    Mean,
    Ref,
}

impl FromStr for FillMissingScore {
    type Err = String;
    fn from_str(str: &str) -> Result<Self, Self::Err> {
        match str {
            "mode" => Ok(FillMissingScore::Mode),
            "mean" => Ok(FillMissingScore::Mean),
            "ref" => Ok(FillMissingScore::Ref),
            _ => Err(format!("Unknown FillMissingScore: {}", str)),
        }
    }
}

/// For loading group
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum FillMissingGroup {
    Ref,
}

impl FromStr for FillMissingGroup {
    type Err = String;
    fn from_str(str: &str) -> Result<Self, Self::Err> {
        match str {
            "ref" => Ok(FillMissingGroup::Ref),
            _ => Err(format!("Unknown FillMissingGroup: {}", str)),
        }
    }
}

// 8 x 1 bit
//type B8 = u8;
#[derive(Clone)]
pub struct Dataset {
    genot: Genot,
    snvs: Snvs,
    samples: Samples,
}

impl Dataset {
    // Do not use R since buf should be used twice and File::open() cannot be copied...
    //pub fn new<R: std::io::Read + Copy>(

    // add assertion for training
    pub fn new_datasetfile_training(
        dfile: &DatasetFile, // should be immutable for simplicity
        use_sample_val: bool,
        filt_snv: Option<&[bool]>,
        fill_missing: Option<FillMissing>,
        //fill_missing_mode: bool,
        fill_missing_group: Option<FillMissingGroup>,
        make_major_a2_train: bool,
        snvs_train: Option<&Snvs>, //for use_missing of vali
        //fill_missing: bool,
        //filt_snv: Option<&[bool]>,
        //snvs_train: Option<&Snvs>, //for use_missing of vali
        //make_major_a2_train: bool,
        mem: Option<usize>,
    ) -> Self {
        //if fill_missing_mode && use_sample_val && snvs_train.is_none() {
        if fill_missing.is_some_and(|x| x == FillMissing::Mode)
            && use_sample_val
            && snvs_train.is_none()
        {
            panic!("fill_missing and use_sample_val require snvs_train for maf");
        }

        let dataset = Self::new_datasetfile(
            dfile,
            use_sample_val,
            filt_snv,
            fill_missing,
            fill_missing_group,
            make_major_a2_train,
            snvs_train,
            mem,
        );

        // check if ys exist
        if dataset.samples().phe().is_none() {
            panic!("Could not load sample phenotype.");
        }

        dataset
    }

    pub fn new_datasetfile(
        dfile: &DatasetFile, // should be immutable for simplicity
        use_sample_val: bool,
        filt_snv: Option<&[bool]>,
        fill_missing: Option<FillMissing>,
        //fill_missing_mode: bool,
        fill_missing_group: Option<FillMissingGroup>,
        make_major_a2_train: bool,
        snvs_train: Option<&Snvs>,
        mem: Option<usize>,
    ) -> Self {
        let sample_buf = if use_sample_val {
            dfile.sample_val_buf()
        } else {
            dfile.sample_buf()
        };
        Self::new(
            dfile.fins_genot(),
            dfile.phe_buf(),
            dfile.phe_name(),
            dfile.cov_name(),
            dfile.snv_buf(),
            dfile.agg_snv_buf(),
            dfile.group_snv_buf(),
            sample_buf,
            filt_snv,
            fill_missing,
            fill_missing_group,
            make_major_a2_train,
            snvs_train,
            mem,
        )
    }

    /// Input is one of the following
    /// 1. plink2 + fin_phe
    /// 2. plink2 (phe) + fin_phe (cov) : phe_name is in psam
    /// 3. plink2 (cov and phe in .psam) : phe_buf is option
    /// 4. plink1 + fin_phe (cov and phe)
    /// 5. plink1 (phe) + fin_phe (cov) : phe_name is option
    fn new(
        fins_genot: &GenotFiles,
        //fin: &Path,
        //gfmt: GenotFormat,
        phe_buf: Option<&[u8]>,
        phe_name: Option<&str>,
        cov_name: Option<&str>,
        snv_buf: Option<&[u8]>,
        agg_snv_buf: Option<&[u8]>,
        group_snv_buf: Option<&[u8]>,
        sample_buf: Option<&[u8]>,
        // TODO: merge filt_snv and fin_snv into use_snvs,
        filt_snv: Option<&[bool]>,
        fill_missing: Option<FillMissing>,
        fill_missing_group: Option<FillMissingGroup>,
        //fill_missing_mode: bool,
        make_major_a2_train: bool,
        snvs_train: Option<&Snvs>, //for use_missing of vali
        mem: Option<usize>,
    ) -> Self {
        fins_genot.check_valid_open();
        //genot_io::check_valid_fin(fin, gfmt);

        // load snvs
        // len(file_snv_allele_idx_in) = m_in
        let (snvs_in, file_snv_allele_idx_in) = genot_io::load_snvs(fins_genot);
        //let snvs_in = genot_io::load_snvs(fins_genot);
        let (mut use_snvs, mut m_snv) = snv::make_use_snvs_buf(snv_buf, &snvs_in);
        log::debug!("m_snv: {}", m_snv);

        // filter snv only for use_snvs not for set_snvs
        if let Some(filt_snv_) = filt_snv {
            log::debug!("filt_snv before m: {}", m_snv);
            assert_eq!(filt_snv_.len(), m_snv);
            use_snvs = vec::and_in_bool_vec(&use_snvs, filt_snv_);
            m_snv = vec::count_true(&use_snvs);
            log::debug!("filt_snv m_snv: {}", m_snv);
        }

        // load agg snvs
        let (use_snvs_agg, m_agg) = snv::make_use_snvs_agg_buf(agg_snv_buf, &snvs_in);

        // load snv set
        //let (use_set_snvs, m_set) = snv::make_set_snvs_buf(set_snv_buf, &snvs_in);
        //assert_eq!(use_set_snvs.len(), use_snvs.len());
        //log::debug!("m_set: {}", m_set);

        // group is fixed, so create OR of snvs here
        let (group_to_m_in, m_group) = snv::make_group_to_m_in_buf(group_snv_buf, &snvs_in);
        //let set_snvs_use = snv::load_set_snvs_buf(set_snv_buf.unwrap());
        //let set_to_m_in = snv::load_set_to_m_in(set_snvs_use, &snvs_in);

        log::debug!("group_to_m_in: {:?}", group_to_m_in);

        // later
        //if m_set != 0 {
        //    //  set snvs exist
        //    use_snvs = vec::or_bool_vec(&use_snvs, &use_set_snvs);
        //    m_snv = vec::count_true(&use_snvs);
        //    log::debug!("m_snv after merge with snv set: {}", m_snv);
        //}

        if (m_snv == 0) && (m_agg == 0) && (m_group == 0) {
            panic!("Using SNVs are zero. Please check fin_snv.")
        }
        log::debug!("m_snv: {}", m_snv);
        log::debug!("m_group: {}", m_group);

        let (use_samples, n) = sample::make_use_samples_buf(sample_buf, fins_genot);
        if n == 0 {
            panic!("Using samples are zero. Please check fin_sample.")
        }
        log::debug!("n: {}", n);

        Self::new_use_vec(
            fins_genot,
            phe_buf,
            phe_name,
            cov_name,
            group_snv_buf,
            agg_snv_buf,
            Some(&use_snvs),
            Some(&use_snvs_agg),
            &file_snv_allele_idx_in,
            //Some(&use_set_snvs),
            Some(&use_samples),
            fill_missing,
            fill_missing_group,
            make_major_a2_train,
            snvs_train,
            mem,
            group_to_m_in,
        )
    }

    // TODO: clean: now loading Snvs before and here.
    // pass Option<Snvs>
    pub fn new_use_vec(
        fins_genot: &GenotFiles,
        phe_buf: Option<&[u8]>,
        phe_name: Option<&str>,
        cov_name: Option<&str>,
        group_snv_buf: Option<&[u8]>,
        agg_snv_buf: Option<&[u8]>,
        use_snvs: Option<&[bool]>,
        use_snvs_agg: Option<&[bool]>,
        file_snv_allele_idx_in: &[(usize, usize, usize)],
        //use_set_snvs: Option<&[bool]>,
        use_samples: Option<&[bool]>,
        fill_missing: Option<FillMissing>,
        fill_missing_group: Option<FillMissingGroup>,
        //fill_missing_mode: bool,
        make_major_a2_train: bool,
        snvs_train: Option<&Snvs>, //for fill_missing and make_major_a2_train of vali
        mem: Option<usize>,
        group_to_m_in: Option<Vec<Vec<usize>>>,
        // group_to_m_in: Option<HashMap<usize, Vec<usize>>>,
    ) -> Self {
        let start = Instant::now();

        //let use_snvs_or_set = match use_snvs {
        //    Some(use_snvs) => match use_set_snvs {
        //        // use both
        //        Some(use_set_snvs) => Some(vec::or_bool_vec(use_snvs, use_set_snvs)),
        //        // use_snvs only
        //        None => Some(use_snvs.to_vec()),
        //    },
        //    // use all snvs
        //    None => None,
        //};

        let use_snvs_or_agg: Option<Vec<bool>> = match (use_snvs, use_snvs_agg) {
            (Some(use_snvs), Some(use_snvs_agg)) => Some(vec::or_bool_vec(use_snvs, use_snvs_agg)),
            (Some(x), None) | (None, Some(x)) => Some(x.to_vec()),
            (None, None) => None,
        };

        let m = match use_snvs_or_agg {
            Some(ref x) => vec::count_true(&x),
            None => genot_io::compute_num_snv(fins_genot),
        };

        // let m = match use_snvs {
        //     Some(x) => vec::count_true(x),
        //     None => genot_io::compute_num_snv(fins_genot),
        // };

        let n = match use_samples {
            Some(x) => vec::count_true(x),
            None => genot_io::compute_num_sample(fins_genot),
        };

        // do not fill missing here since validation dataset requires maf to fill missing; below
        let genot = genot_io::load::generate_genot(
            fins_genot,
            m,
            n,
            // use_snvs,
            use_snvs_or_agg.clone().as_deref(),
            //use_snvs,
            file_snv_allele_idx_in,
            group_to_m_in,
            use_samples,
            None,
            // false,
            fill_missing_group,
            mem,
        );
        //genot_io::load::generate_genot(fin_genot, m, n, use_snvs, use_samples, false, mem);

        let samples = Samples::new_plink(
            fins_genot,
            phe_buf,
            phe_name,
            cov_name,
            use_samples,
            Some(n),
            true,
        );

        let snvs = Snvs::new_plink_use_snvs_and_group(
            fins_genot,
            use_snvs_or_agg.clone().as_deref(),
            // use_snvs,
            None,
            agg_snv_buf,
            group_snv_buf,
        );
        // Snvs::new_plink_use_snvs_and_set(fin_genot, use_snvs, None, set_snv_buf.unwrap());
        //let snvs = Snvs::new_plink_use_snvs(fin_genot, use_snvs, Some(m_snv));

        log::debug!(
            "It took {} seconds to create Dataset.",
            start.elapsed().as_secs()
        );

        let mut dataset = Dataset {
            genot,
            snvs,
            samples,
        };
        dataset.check();

        dataset.compute_maf();

        if make_major_a2_train {
            if let Some(snvs_train) = snvs_train {
                // for val
                dataset.set_major_a2(Some(snvs_train));
            } else {
                // for training
                dataset.set_major_a2(None);
            }
        }

        // fill missing to mode here
        // for training; use maf assuming hwe
        // for validation; use the training maf
        //if fill_missing {
        if fill_missing.is_some() {
            match fill_missing.unwrap() {
                FillMissing::Mode => {
                    if let Some(snvs_train) = snvs_train {
                        // for val
                        // unwrap to raise error when None
                        let mafs = snvs_train.freqs_a1().unwrap();
                        dataset.fill_missing_mode_maf(Some(mafs));
                    } else {
                        // for training
                        dataset.fill_missing_mode_maf(None);
                    }
                }
                FillMissing::Ref => {
                    // In most cases, assume REF is major
                    dataset.fill_missing_ref();
                }
            }
        }

        //if make_major_a2_train {
        //    if let Some(snvs_train) = snvs_train {
        //        // for val
        //        dataset.set_major_a2(Some(snvs_train));
        //    } else {
        //        // for training
        //        dataset.set_major_a2(None);
        //    }
        //}

        dataset
    }

    // for test
    pub fn new_field_phe(genot: Genot, samples: Samples, snvs: Snvs) -> Self {
        Self {
            genot,
            snvs,
            samples,
        }
    }

    /// for prune SNVs by loss
    pub fn extract_snvs(self, use_snvs: &[bool]) -> Self {
        let Dataset {
            genot,
            snvs,
            samples,
        } = self;

        let genot_use = genot.extract_snvs(use_snvs);
        let snvs_use = snvs.extract_snvs(use_snvs);

        let dataset = Dataset {
            genot: genot_use,
            snvs: snvs_use,
            samples,
        };
        dataset.check();
        dataset
    }

    pub fn genot(&self) -> &Genot {
        &self.genot
    }
    pub fn genot_mut(&mut self) -> &mut Genot {
        &mut self.genot
    }

    pub fn genot_split_group(&self) -> (GenotRef, GenotRef) {
        let m_snv = self.snvs().snvs_snv_n();
        let (g_snv, g_group) = self.genot.split_genot(m_snv);
        assert_eq!(g_group.m(), self.snvs().group_n());

        (g_snv, g_group)
    }

    pub fn snvs(&self) -> &Snvs {
        &self.snvs
    }
    pub fn snvs_mut(&mut self) -> &mut Snvs {
        &mut self.snvs
    }

    pub fn snv_ids_group(&self) -> Vec<SnvId> {
        let m_snv = self.snvs().snvs_snv_n();
        let snv_ids = self.snvs().snv_ids();
        let snv_ids_group = snv_ids[m_snv..].to_vec();
        // snv_ids[0..m_snv].to_vec()
        assert_eq!(snv_ids_group.len(), self.snvs().group_n());
        snv_ids_group
    }

    pub fn samples(&self) -> &Samples {
        &self.samples
    }

    // in SnvTrait
    //fn snvs(&self) -> &Snvs {
    //    &self.snvs
    //}
    //fn samples(&self) -> &Samples {
    //    &self.samples
    //}
    // move to impl DatasetBiNew{}
    fn check(&self) {
        // check if samples_n, snvs_n are the same
    }

    pub fn fill_missing_mode_maf(&mut self, mafs: Option<&Vec<f64>>) {
        let mafs_v: Vec<f64>;
        let mafs = if let Some(mafs) = mafs {
            // for validation
            mafs
        } else {
            // for training
            // if not using clone(), genot_mut() is error
            mafs_v = self.snvs().freqs_a1().unwrap().clone();
            &mafs_v
            // self.snvs().mafs().unwrap()
        };
        //let genot = self.genot_mut();
        self.genot_mut()
            .iter_snv_mut()
            .zip(mafs.iter())
            .par_bridge()
            .for_each(|(mut g_snv, maf)| g_snv.fill_missing_mode_maf(*maf));
        //.for_each(|(mut g_snv, maf)| io_genot::load::fill_missing_maf(&mut g_snv, *maf));
    }

    pub fn fill_missing_ref(&mut self) {
        self.genot_mut()
            .iter_snv_mut()
            .par_bridge()
            .for_each(|mut g_snv| g_snv.fill_missing_ref());
        //.for_each(|(mut g_snv, maf)| io_genot::load::fill_missing_maf(&mut g_snv, *maf));
    }

    pub fn set_major_a2(&mut self, snvs: Option<&Snvs>) {
        // 1. set major as a2 in snvs
        // 2. same in genot
        // 3. update maf

        let mafs_v: Vec<f64>;
        let mafs = if let Some(snvs) = snvs {
            // for validation
            snvs.freqs_a1().unwrap()
        } else {
            // for training
            mafs_v = self.snvs().freqs_a1().unwrap().clone();
            &mafs_v
        };

        // reverse snvs
        self.snvs_mut()
            .snv_ids_mut()
            .iter_mut()
            .zip(mafs.iter())
            .filter(|(_, maf)| **maf > 0.5)
            .par_bridge()
            .for_each(|(snv_id, _)| snv_id.reverse_alleles());

        // reverse genot
        self.genot_mut()
            .iter_snv_mut()
            .zip(mafs.iter())
            .filter(|(_, maf)| **maf > 0.5)
            .par_bridge()
            .for_each(|(mut g_snv, _)| g_snv.reverse_allele());

        // update maf
        self.compute_maf();

        // check all maf<0.5 for training?
    }

    pub fn compute_maf(&mut self) {
        let genot = self.genot();
        let m = genot.m();
        let mut mafs = vec![f64::NAN; m];
        mafs.iter_mut()
            .zip(genot.iter_snv())
            .par_bridge()
            .for_each(|(maf, g_snv)| *maf = g_snv.maf());

        self.snvs_mut().set_maf(mafs)
    }

    pub fn skip_neg_agg(&mut self) {
        // skip negative SNVs in aggregation
        let agg_to_m = self.snvs_mut().agg_to_m();

        if agg_to_m.is_none() {
            return;
        }

        let agg_to_m = agg_to_m.unwrap().clone();

        // vector of m in agg_to_m
        let ms_use: Vec<usize> = agg_to_m.iter().flat_map(|x| x.iter()).cloned().collect();

        let phe = self.samples().phe().unwrap();
        let case_n = phe.count();
        let genot_matrix = self.genot().genot_inner();

        // dominant
        let digit = 0usize;

        let sign_use: Vec<bool> = ms_use
            .iter()
            .map(|x| {
                assoc::judge_sign_odds_ratio(&genot_matrix.as_cbits_ref(*x, digit), phe, case_n)
            })
            .collect();

        // snv with negative sign
        let ms_exclude = sign_use
            .iter()
            .zip(ms_use.iter())
            .filter(|(sign, _)| !**sign)
            .map(|(_, m)| *m)
            .collect::<HashSet<usize>>();

        log::info!(
            "Excluding {} SNVs with negative sign in aggregation.",
            ms_exclude.len()
        );

        // TODO
        // unimplemented!("output agg_to_m with excluded SNVs");

        let agg_to_m = agg_to_m
            .into_iter()
            .map(|x| {
                x.into_iter()
                    .filter(|m| !ms_exclude.contains(m))
                    .collect::<Vec<usize>>()
            })
            .collect::<Vec<Vec<usize>>>();

        self.snvs_mut().set_agg_to_m(Some(agg_to_m))
    }

    pub fn new_datasetfile_score<W: WgtTrait>(
        dfile: &DatasetFile,
        wgts: &mut [W],
        //wgts: &[WgtBoost],
        // fill_missing_in_dataset: bool,
        //fill_missing: bool,
        fill_missing_score: Option<FillMissingScore>,
        fill_missing_group: Option<FillMissingGroup>,
        allow_nonexist_snv: bool,
        use_snv_pos: bool,
        mem: Option<usize>,
    ) -> Self {
        Self::new_score(
            dfile.fins_genot(),
            //dfile.fin(),
            //dfile.gfmt(),
            dfile.phe_buf(),
            dfile.cov_name(),
            // dfile.group_snv_buf(),
            dfile.sample_buf(),
            wgts,
            // fill_missing_in_dataset,
            //fill_missing,
            fill_missing_score,
            fill_missing_group,
            allow_nonexist_snv,
            use_snv_pos,
            mem,
        )
    }

    // for boosting
    // merge partly to new()
    pub fn new_score<W: WgtTrait>(
        fins_genot: &GenotFiles,
        //fin: &Path,
        //gfmt: GenotFormat,
        phe_buf: Option<&[u8]>,
        //fin_phe: Option<&Path>,
        //phe_name: Option<&str>,
        cov_name: Option<&str>,
        // group_snv_buf: Option<&[u8]>,
        sample_buf: Option<&[u8]>,
        //extract_sample_buf: Option<&[u8]>,
        //fin_sample: Option<&Path>,
        //fin_cov: Option<&Path>,
        wgts: &mut [W],
        //wgts: &[WgtBoost],
        // fill_missing_in_dataset: bool,
        //use_missing: bool, // use WgtBoosts and use wgts.use_missing()
        // fill_missing: bool,
        fill_missing_score: Option<FillMissingScore>,
        fill_missing_group: Option<FillMissingGroup>,
        allow_nonexist_snv: bool,
        use_snv_pos: bool,
        mem: Option<usize>,
        // group_to_m_in: Option<Vec<Vec<usize>>>,
    ) -> Self {
        //let m_in: usize = io_genot::compute_num_snv(fin, gfmt).unwrap();
        //log::debug!("m_in {}", m_in);
        //let n_in: usize = io_genot::compute_num_sample(fin, gfmt).unwrap();
        //log::debug!("n_in {}", n_in);

        let (use_samples, n) = sample::make_use_samples_buf(sample_buf, fins_genot);
        //let (n, use_samples) = sample::make_use_samples(fin_sample, fin, gfmt);
        if n == 0 {
            panic!("Using samples are zero. Please check fin_sample.")
        }

        let samples = Samples::new_plink(
            fins_genot,
            phe_buf,
            None,
            cov_name,
            Some(&use_samples),
            Some(n),
            false,
        );

        // set genotype index in wgt
        //let use_missing = true;
        let genot = genot_io::load_score::generate_genot_for_score_boosting(
            fins_genot,
            //gfmt,
            wgts,
            n,
            Some(&use_samples),
            // fill_missing,
            // fill_missing_in_dataset,
            allow_nonexist_snv,
            use_snv_pos,
            fill_missing_score,
            fill_missing_group,
            mem,
        );

        Dataset {
            genot,
            // unnecessary since index is in WgtKInd
            snvs: Snvs::new_empty(),
            samples,
        }
    }

    /// For multiple wgts
    pub fn new_datasetfile_score_genetics(
        dfile: &DatasetFile,
        wgts: &mut [Wgts],
        // fill_missing_in_dataset: bool,
        fill_missing_group: Option<FillMissingGroup>,
        allow_nonexist_snv: bool,
        use_snv_pos: bool,
        mem: Option<usize>,
    ) -> Self {
        Self::new_score_genetics(
            dfile.fins_genot(),
            dfile.phe_buf(),
            dfile.cov_name(),
            dfile.sample_buf(),
            dfile.freq_buf(),
            wgts,
            // fill_missing_in_dataset,
            fill_missing_group,
            allow_nonexist_snv,
            use_snv_pos,
            mem,
        )
    }

    // Usually not fill missing here. fill_missing_in_dataset is for backward compatibility.
    // for genetics::score()
    pub fn new_score_genetics(
        fins_genot: &GenotFiles,
        //fin: &Path,
        //gfmt: GenotFormat,
        phe_buf: Option<&[u8]>,
        cov_name: Option<&str>,
        sample_buf: Option<&[u8]>,
        freq_buf: Option<&[u8]>,
        wgts: &mut [Wgts],
        //wgts: &[WgtBoost],
        // fill_missing_in_dataset: bool,
        fill_missing_group: Option<FillMissingGroup>,
        allow_nonexist_snv: bool,
        use_snv_pos: bool,
        mem: Option<usize>,
    ) -> Self {
        //let m_in: usize = io_genot::compute_num_snv(fin, gfmt).unwrap();
        //log::debug!("m_in {}", m_in);
        //let n_in: usize = io_genot::compute_num_sample(fin, gfmt).unwrap();
        //log::debug!("n_in {}", n_in);

        //let (n, use_samples) = sample::make_use_samples_buf(sample_buf, fin, gfmt);
        let (use_samples, n) = sample::make_use_samples_buf(sample_buf, fins_genot);
        //let (n, use_samples) = sample::make_use_samples(fin_sample, fin, gfmt);
        if n == 0 {
            panic!("Using samples are zero. Please check fin_sample.")
        }

        // TODO: make (string, string) as key of hashmap
        //let samples_id = plink::load_samples_id(fin, &use_samples);

        let samples = Samples::new_plink(
            fins_genot,
            phe_buf,
            None,
            cov_name,
            Some(&use_samples),
            Some(n),
            false,
        );

        // set genotype index in wgt
        // TODO: argparse
        //let use_missing = false;
        // TOFIX: fill_missing using MAF (in wgt or) provided.
        //let fill_missing_in_test = true;
        let genot = genot_io::load_score::generate_genot_for_score_multiwgts(
            fins_genot,
            wgts,
            freq_buf,
            n,
            Some(&use_samples),
            // fill_missing_in_dataset,
            fill_missing_group,
            allow_nonexist_snv,
            use_snv_pos,
            mem,
        );

        // add freq to wgts
        //let freqs_in = genot_io::load_freq(freq_buf);
        // index
        // is_reversed

        // TMP
        //for snv in genot.iter_snv() {
        //    println!("snv {:?}", &snv.vals()[..10]);
        //}

        Dataset {
            genot,
            // unnecessary since index is in WgtKInd
            snvs: Snvs::new_empty(),
            samples,
        }
    }
}

#[cfg(test)]
mod tests {
    //use super::*;
}
