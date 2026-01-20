//! since for plink2, columns are not fixed and troublesome to detect them every time.
//! when creating fgenot class, detect columns etc.
//
// TODO: split loading from file and others for test

pub mod load;
pub mod load_score;

use crate::{sample, Snvs};
use crate::{textfile, vec, SnvId};

use serde::Deserialize;
use std::collections::{HashMap, HashSet};
use std::ffi::{OsStr, OsString};
use std::io::{self, BufRead};
use std::path::{Path, PathBuf};

// 8 x 1 bit
//type B8 = u8;
// 4 x 2 bit
//type B8_2 = u8;

//pub type GenotFiles = Vec<GenotFile>;

#[derive(Clone, PartialEq, Eq, Hash, Debug)]
pub struct GenotFiles {
    files: Vec<GenotFile>,
}

impl GenotFiles {
    pub fn new(files: Vec<GenotFile>) -> Self {
        let fins_genot = Self { files };
        fins_genot.check();
        fins_genot
    }

    pub fn new_plink1(files: Vec<PathBuf>) -> Self {
        let files = files.into_iter().map(|x| GenotFile::Plink1(x)).collect();
        let fins_genot = Self { files };
        fins_genot.check();
        fins_genot
    }

    pub fn new_plink2(files: Vec<PathBuf>) -> Self {
        let files = files.into_iter().map(|x| GenotFile::Plink2(x)).collect();
        let fins_genot = Self { files };
        fins_genot.check();
        fins_genot
    }

    pub fn new_plink2vzs(files: Vec<PathBuf>) -> Self {
        let files = files.into_iter().map(|x| GenotFile::Plink2Vzs(x)).collect();
        let fins_genot = Self { files };
        fins_genot.check();
        fins_genot
    }

    pub fn check(&self) {
        // Check if all files are the same format
        if self.files.len() == 0 {
            panic!("No genot files were given.");
        }
        let fin0 = &self.files[0];
        for fin_genot in self.files() {
            // compare enum type only
            if std::mem::discriminant(fin0) != std::mem::discriminant(fin_genot) {
                panic!(
                    "All files should be the same format: {:?} vs {:?}",
                    fin0, fin_genot
                );
            }
        }
    }

    pub fn files(&self) -> &Vec<GenotFile> {
        &self.files
    }

    /// Check genot files can be open
    /// If chrom is split,
    ///     fin_chrom.bim, .fam, .bed exists
    /// allow some of chorom not exist
    pub fn check_valid_open(&self) {
        for fin_genot in self.files().iter() {
            fin_genot.check_valid_open();
        }
    }
}

// multi-allelic
// computationally heavy
pub fn compute_num_snv(fins_genot: &GenotFiles) -> usize {
    let (snvs, _) = load_snvs(fins_genot);
    snvs.len()
}

// multi-allelic
// computationally heavy
pub fn compute_num_snv_file(fin_genot: &GenotFile, file_i: Option<usize>) -> usize {
    let (snvs, _) = load_snvs_file(fin_genot, file_i);
    snvs.len()
}

/// assume only part of chrom files could exist
// pub fn compute_num_snv(fins_genot: &GenotFiles) -> usize {
//     fins_genot
//         .files()
//         .iter()
//         .map(|fin_genot| compute_num_snv_file(fin_genot))
//         .sum()
// }

// expand multi-allelics snvs into multiple bi-allelics
//pub fn load_snvs(fins_genot: &GenotFiles) -> Vec<SnvId> {
pub fn load_snvs(fins_genot: &GenotFiles) -> (Vec<SnvId>, Vec<(usize, usize, usize)>) {
    let mut snvs: Vec<SnvId> = vec![];
    let mut files_snv_allele_idx: Vec<(usize, usize, usize)> = vec![];
    for (file_i, fin_genot) in fins_genot.files().iter().enumerate() {
        let (mut snv_file, mut file_snv_allele_idx) = load_snvs_file(fin_genot, Some(file_i));
        snvs.append(&mut snv_file);
        files_snv_allele_idx.append(&mut file_snv_allele_idx);
    }
    if snvs.len() == 0 {
        panic!("No snvs were loaded from fins_genot");
    }

    let snvs_rs: HashSet<String> = HashSet::from_iter(
        snvs.iter()
            .map(|x| x.ida().to_owned())
            .collect::<Vec<String>>(),
    );

    if snvs.len() != snvs_rs.len() {
        let mut map: HashSet<String> = HashSet::new();

        let mut dups = vec![];
        for snv_id in snvs.iter() {
            if !map.insert(snv_id.ida().to_string()) {
                dups.push(snv_id.ida().to_string());
            }
        }

        panic!("Duplicate snv id in plink2: {:?}", dups);
    }

    // files_snv_allele_idx could be empty for plink1
    (snvs, files_snv_allele_idx)
}

pub fn compute_num_sample(fins_genot: &GenotFiles) -> usize {
    compute_num_sample_file(&fins_genot.files()[0])
}

pub fn load_samples_id(fins_genot: &GenotFiles, use_samples: Option<&[bool]>) -> Vec<String> {
    load_samples_id_file(&fins_genot.files()[0], use_samples)
}

pub fn load_ys_buf(
    fins_genot: &GenotFiles,
    phe_buf: Option<&[u8]>,
    phe_name: Option<&str>,
    sample_id_to_n: &HashMap<String, usize>,
) -> Option<Vec<bool>> {
    load_ys_buf_file(&fins_genot.files()[0], phe_buf, phe_name, sample_id_to_n)
}

#[derive(Clone, PartialEq, Eq, Hash, Debug)]
pub enum GenotFile {
    Plink1(PathBuf),
    Plink2(PathBuf),
    Plink2Vzs(PathBuf),
}

impl GenotFile {
    pub fn file(&self) -> &PathBuf {
        match self {
            Self::Plink1(fin) | Self::Plink2(fin) | Self::Plink2Vzs(fin) => fin,
        }
    }

    pub fn sample_file(&self) -> PathBuf {
        match self {
            Self::Plink1(x) => add_ext(x, "fam"),
            Self::Plink2(x) | Self::Plink2Vzs(x) => add_ext(x, "psam"),
        }
    }

    // TODO: move to GenotFiles?
    pub fn sample_file_exist_chrom(&self) -> Option<PathBuf> {
        Some(self.sample_file())
        //if !self.judge_split_chrom() {
        //    Some(self.sample_file(None))
        //} else {
        //    let chrom = self.choose_a_chrom_exist();
        //    chrom.map(|x| self.sample_file(Some(&x)))
        //}
    }

    pub fn snv_file(&self) -> PathBuf {
        match self {
            Self::Plink1(x) => add_ext(x, "bim"),
            Self::Plink2(x) => add_ext(x, "pvar"),
            Self::Plink2Vzs(x) => add_ext(x, "pvar.zst"),
        }
    }

    pub fn genotype_file(&self) -> PathBuf {
        match self {
            Self::Plink1(x) => add_ext(x, "bed"),
            Self::Plink2(x) | Self::Plink2Vzs(x) => add_ext(x, "pgen"),
        }
    }

    // pub fn replace_file_chrom(&self, chrom: Option<&Chrom>) -> PathBuf {
    pub fn replace_file_chrom(&self, chrom: Option<&str>) -> PathBuf {
        let fin = self.file();
        match chrom {
            Some(chrom) => {
                let mut dname = fin.parent().unwrap().to_path_buf();
                let fname = fin.file_name().unwrap().to_owned().into_string().unwrap();
                let fname = fname.replace("%", &chrom.to_string());
                dname.push(fname);
                dname
            }
            None => fin.to_owned(),
        }
    }

    pub fn check_valid_open(&self) {
        textfile::check_open_file(&self.sample_file());
        textfile::check_open_file(&self.snv_file());
        textfile::check_open_file(&self.genotype_file());
    }

    // /// Check genot files can be open
    // /// If chrom is split,
    // ///     fin_chrom.bim, .fam, .bed exists
    // /// allow some of chorom not exist
    //pub fn check_valid_open(&self) {
    //    //genot_io::check_valid_fin(self.fin(), self.gfmt());
    //    if !self.judge_split_chrom() {
    //        textfile::check_open_file(&self.sample_file(None));
    //        textfile::check_open_file(&self.snv_file(None));
    //        textfile::check_open_file(&self.genotype_file(None));
    //    } else {
    //        for chrom_i in Chrom::variants().iter() {
    //            let mut is_exists = vec![false, false, false];
    //            is_exists[0] = textfile::able_open_file(&self.sample_file(Some(chrom_i)));
    //            is_exists[1] = textfile::able_open_file(&self.snv_file(Some(chrom_i)));
    //            is_exists[2] = textfile::able_open_file(&self.genotype_file(Some(chrom_i)));

    //            // check if all elements of is_exists are the same
    //            // allow all true or all false
    //            if !is_exists.iter().all(|v| *v == is_exists[0]) {
    //                if !is_exists[0] {
    //                    panic!(
    //                        "Sample file does not exist: {:?}",
    //                        self.sample_file(Some(chrom_i))
    //                    );
    //                } else if !is_exists[1] {
    //                    panic!(
    //                        "Snv file does not exist: {:?}",
    //                        self.snv_file(Some(chrom_i))
    //                    );
    //                } else if !is_exists[2] {
    //                    panic!(
    //                        "Genotype file does not exist: {:?}",
    //                        self.genotype_file(Some(chrom_i))
    //                    );
    //                }
    //            }
    //        }
    //    }
    //}

    // pub fn judge_split_chrom(&self) -> bool {
    //     let fin = match self {
    //         Self::Plink1(fin) | Self::Plink2(fin) | Self::Plink2Vzs(fin) => fin,
    //     };

    //     let fname = fin.file_name().unwrap().to_owned().into_string().unwrap();
    //     fname.contains("%")
    // }

    // /// get smallest number of chrom whose plink exists.
    // /// mainly used for get .fam
    // /// return if fin is not split chrom or all chrom does not exist
    //pub fn choose_a_chrom_exist(&self) -> Option<Chrom> {
    //    if !self.judge_split_chrom() {
    //        None
    //    } else {
    //        for chrom_i in Chrom::variants().iter() {
    //            let fin_fam = self.sample_file(Some(chrom_i));
    //            if textfile::able_open_file(&fin_fam) {
    //                return Some(chrom_i.clone());
    //            }
    //        }
    //        None
    //    }
    //}
}

// TO DELETE
//#[derive(Copy, Clone, PartialEq, Eq, Hash, Debug)]
//pub enum GenotFormat {
//    Plink1,
//    Plink2,
//    Plink2Vzs,
//}
//
//impl GenotFormat {
//    // use match; try not to use this
//    //fn is_plink2(self) -> bool {
//    //    match self {
//    //        Self::Plink => false,
//    //        Self::Plink2 | Self::Plink2Vzs => true,
//    //    }
//    //}
//    //fn is_plink(self) -> bool {
//    //    match self {
//    //        Self::Plink => true,
//    //        Self::Plink2 | Self::Plink2Vzs => false,
//    //    }
//    //}
//}

// use GenotFile::judge_split_chrom()
// error:
////pub fn judge_split_chrom<P: AsRef<Path>>(fin: P) -> bool {
//pub fn judge_split_chrom(fin: &Path) -> bool {
//    //fin.contains("%")
//    let fname = fin.file_name().unwrap().to_owned().into_string().unwrap();
//    fname.contains("%")
//}

// use GenotFile::replace_file_chrom()
//pub fn replace_fname(fin: &Path, chrom: Option<&Chrom>) -> String {
//pub fn replace_fname(fin: &Path, chrom: Option<&Chrom>) -> PathBuf {
//    match chrom {
//        Some(chrom) => {
//            let mut dname = fin.parent().unwrap().to_path_buf();
//            let fname = fin.file_name().unwrap().to_owned().into_string().unwrap();
//            let fname = fname.replace("%", &chrom.to_string());
//            dname.push(fname);
//            dname
//        }
//        None => fin.to_owned(),
//    }
//    //fin.replace("%", &chrom.to_string())
//}

//pub fn fname_chrom(fin: &Path, chrom: Option<&Chrom>) -> PathBuf {
//    replace_fname(fin, chrom)
//}

// move to mod path
// pathbuf ver. of fout+"_abc"
pub fn add_name(f: &Path, ext: impl AsRef<OsStr>) -> PathBuf {
    let mut fnew: OsString = f.into();
    fnew.push(ext.as_ref());
    fnew.into()
}

// move to mod path
// https://internals.rust-lang.org/t/pathbuf-has-set-extension-but-no-add-extension-cannot-cleanly-turn-tar-to-tar-gz/14187/11
pub fn add_ext(f: &Path, ext: impl AsRef<OsStr>) -> PathBuf {
    let mut fnew: OsString = f.into();
    fnew.push(".");
    fnew.push(ext.as_ref());
    fnew.into()
}

// use GenotFile::sample_file()
//pub fn fname_plinks_sample(fin: &Path, gfmt: GenotFormat, chrom: Option<&Chrom>) -> PathBuf {
//    let f = replace_fname(fin, chrom);
//    let ext = match gfmt {
//        GenotFormat::Plink1 => "fam",
//        GenotFormat::Plink2 | GenotFormat::Plink2Vzs => "psam",
//    };
//    add_ext(&f, ext)
//}

//pub fn fname_sample(fin: &Path, chrom: Option<&Chrom>) -> PathBuf {
//    let f = replace_fname(fin, chrom);
//    add_ext(&f, "fam")
//    //replace_fname(fin, chrom) + ".fam"
//}

// use GenotFile::snv_file()
//pub fn fname_plinks_snv(fin: &Path, gfmt: GenotFormat, chrom: Option<&Chrom>) -> PathBuf {
//    let f = replace_fname(fin, chrom);
//    let ext = match gfmt {
//        GenotFormat::Plink1 => "bim",
//        GenotFormat::Plink2 => "pvar",
//        GenotFormat::Plink2Vzs => "pvar.zst",
//    };
//    add_ext(&f, ext)
//}

// use GenotFile::genotype_file()
//pub fn fname_plinks_genot(fin: &Path, gfmt: GenotFormat, chrom: Option<&Chrom>) -> PathBuf {
//    let f = replace_fname(fin, chrom);
//    let ext = match gfmt {
//        GenotFormat::Plink1 => "bed",
//        GenotFormat::Plink2 | GenotFormat::Plink2Vzs => "pgen",
//    };
//    add_ext(&f, ext)
//}

//// TODO: move to plink1
//pub fn bed_per_snv_size(n: usize) -> usize {
//    (n + 3) / 4
//}
//
//// TODO: move to plink1
//pub fn calculate_bed_size_genotype(m: usize, n: usize) -> usize {
//    m * bed_per_snv_size(n)
//    //m * ((n + 3) / 4)
//}
//
//// TODO: move to plink1
//pub fn calculate_bed_size(m: usize, n: usize) -> usize {
//    3 + calculate_bed_size_genotype(m, n)
//    //3 + m * ((n + 3) / 4)
//}

///// get smallest number of chrom whose plink exists.
///// mainly used for get .fam
///// return if fin is not split chrom or all chrom does not exist
//pub fn choose_a_chrom_exist(fin: &Path, gfmt: GenotFormat) -> Option<Chrom> {
//    if !judge_split_chrom(fin) {
//        None
//    } else {
//        for chrom_i in Chrom::variants().iter() {
//            let fin_fam = fname_plinks_sample(fin, gfmt, Some(chrom_i));
//            if textfile::able_open_file(&fin_fam) {
//                return Some(chrom_i.clone());
//            }
//        }
//        None
//    }
//}

// use GenotFile::sample_file_exist_chrom()
//pub fn fname_fam_exist_chrom(fin: &Path, gfmt: GenotFormat) -> Option<PathBuf> {
//    if !judge_split_chrom(fin) {
//        Some(fname_plinks_sample(fin, gfmt, None))
//    } else {
//        let chrom = choose_a_chrom_exist(fin, gfmt);
//        chrom.map(|chrom_| fname_plinks_sample(fin, gfmt, Some(&chrom_)))
//    }
//}

// return None for missing chrom file
//pub fn compute_num_sample_chrom(fin_genot: &GenotFile) -> Option<usize> {
pub fn compute_num_sample_file(fin_genot: &GenotFile) -> usize {
    let fin_sample = fin_genot.sample_file();

    let linen = textfile::compute_num_line_text(&fin_sample, None).unwrap();

    let has_header = match fin_genot {
        GenotFile::Plink1(_) => false,
        GenotFile::Plink2(_) | GenotFile::Plink2Vzs(_) => has_fam_header_plink2(&fin_sample),
    };

    if has_header {
        linen - 1
    } else {
        linen
    }
}

//// return None for missing chrom file
//pub fn compute_num_sample_chrom(
//    fin: &Path,
//    gfmt: GenotFormat,
//    chrom: Option<&Chrom>,
//) -> Option<usize> {
//    let fin_sample = fname_plinks_sample(&fin, gfmt, chrom);
//    // replace chrom if any
//    //let fin = replace_fname(fin, chrom);
//    //let fin_sample = fname_plinks_sample(&fin, gfmt, None);
//    let linen = match gfmt {
//        GenotFormat::Plink1 => textfile::compute_num_line_text(&fin_sample, None),
//        GenotFormat::Plink2 | GenotFormat::Plink2Vzs => {
//            textfile::compute_num_line_text(&fin_sample, None)
//        }
//    };
//    let has_header = match gfmt {
//        GenotFormat::Plink1 => false,
//        GenotFormat::Plink2 | GenotFormat::Plink2Vzs => has_bim_header_plink2(&fin_sample, None),
//    };
//
//    let linen = linen.ok();
//
//    if has_header {
//        linen.map(|x| x - 1)
//    } else {
//        linen
//    }
//}

//pub fn compute_num_sample(fin_genot: &GenotFile) -> Option<usize> {
//pub fn compute_num_sample(fin_genot: &GenotFile) -> usize {
//compute_num_sample_file(fin_genot)
//if !fin_genot.judge_split_chrom() {
//    compute_num_sample_chrom(fin_genot, None)
//} else {
//    let chrom = fin_genot.choose_a_chrom_exist();
//    match chrom {
//        Some(chrom_) => compute_num_sample_chrom(fin_genot, Some(&chrom_)),
//        None => None,
//    }
//}
//}

///// assume only part of chrom files could exist
//pub fn compute_num_sample(fin: &Path, gfmt: GenotFormat) -> Option<usize> {
//    // use just one of valid chrom
//    if !judge_split_chrom(fin) {
//        compute_num_sample_chrom(fin, gfmt, None)
//        //let fin_fam = fname_plinks_sample(fin, gfmt, None);
//        //textfile::compute_num_line_text(&fin_fam, None)
//    } else {
//        let chrom = choose_a_chrom_exist(fin, gfmt);
//        match chrom {
//            Some(chrom_) => compute_num_sample_chrom(fin, gfmt, Some(&chrom_)),
//            None => None,
//        }
//    }
//}

pub fn compute_num_snv_allele(file_snv_allele_idx_in: &[(usize, usize, usize)]) -> usize {
    file_snv_allele_idx_in.len()
}

pub fn compute_num_snv_allele_file(
    file_snv_allele_idx_in: &[(usize, usize, usize)],
    file_i: usize,
) -> usize {
    file_snv_allele_idx_in
        .iter()
        .filter(|x| x.0 == file_i)
        .count()
}

// not for multi-allelic snvs
// return None for missing chrom file
//pub fn compute_num_snv_file_chrom(fin_genot: &GenotFile) -> Option<usize> {
pub fn compute_num_snv_file_bi_allelic(fin_genot: &GenotFile) -> usize {
    match fin_genot {
        GenotFile::Plink2(_) | GenotFile::Plink2Vzs(_) => {
            panic!("use compute_num_snv_allele_file()");
        }
        _ => (),
    }

    let fin_snv = fin_genot.snv_file();
    let linen = match fin_genot {
        GenotFile::Plink1(_) | GenotFile::Plink2(_) => {
            textfile::compute_num_line_text(&fin_snv, None)
        }
        GenotFile::Plink2Vzs(_) => textfile::compute_num_line_text(&fin_snv, Some("zst")),
    };
    //let linen = linen.ok();
    let linen = linen.unwrap();

    let header_n = match fin_genot {
        GenotFile::Plink1(_) => 0,
        GenotFile::Plink2(_) => header_line_n_bim_plink2(&fin_snv, None),
        GenotFile::Plink2Vzs(_) => header_line_n_bim_plink2(&fin_snv, Some("zst")),
    };
    linen - header_n

    //if has_header {
    //    linen - 1
    //    //linen.map(|x| x - 1)
    //    //Some(linen - 1)
    //} else {
    //    linen
    //    //Some(linen)
    //}
}

/*
pub fn compute_num_snv_file_chrom(
    fin: &Path,
    gfmt: GenotFormat,
    chrom: Option<&Chrom>,
) -> Option<usize> {
    let fin_snv = fname_plinks_snv(&fin, gfmt, chrom);
    // replace chrom if any
    //let fin = replace_fname(fin, chrom);
    //let fin_snv = fname_plinks_snv(&fin, gfmt, None);
    let linen = match gfmt {
        GenotFormat::Plink1 => textfile::compute_num_line_text(&fin_snv, None),
        GenotFormat::Plink2 => textfile::compute_num_line_text(&fin_snv, None),
        GenotFormat::Plink2Vzs => textfile::compute_num_line_text(&fin_snv, Some("zst")),
    };
    let has_header = match gfmt {
        GenotFormat::Plink1 => false,
        GenotFormat::Plink2 => has_bim_header_plink2(&fin_snv, None),
        GenotFormat::Plink2Vzs => has_bim_header_plink2(&fin_snv, Some("zst")),
    };

    let linen = linen.ok();

    if has_header {
        linen.map(|x| x - 1)
        //Some(linen - 1)
    } else {
        linen
        //Some(linen)
    }
}
*/

///// assume only part of chrom files could exist
//pub fn compute_num_snv(fin_genot: &GenotFile) -> Option<usize> {
//    if !fin_genot.judge_split_chrom() {
//        compute_num_snv_file_chrom(fin_genot, None)
//    } else {
//        let mut num_snv = 0;
//        for chrom_i in Chrom::variants().iter() {
//            let num_snv_chrom = compute_num_snv_file_chrom(fin_genot, Some(chrom_i));
//            let num_snv_chrom = match num_snv_chrom {
//                Some(v) => v,
//                None => 0,
//            };
//            num_snv += num_snv_chrom;
//        }
//
//        if num_snv == 0 {
//            return None;
//        } else {
//            return Some(num_snv);
//        }
//    }
//}

/// assume only part of chrom files could exist
//pub fn compute_num_snv(fin: &Path, gfmt: GenotFormat) -> Option<usize> {
//    if !judge_split_chrom(fin) {
//        compute_num_snv_file_chrom(&fin, gfmt, None)
//        //compute_num_snv_file(&fin, gfmt)
//    } else {
//        let mut num_snv = 0;
//        //for chrom_i in 1..=22 {
//        for chrom_i in Chrom::variants().iter() {
//            let num_snv_chrom = compute_num_snv_file_chrom(fin, gfmt, Some(chrom_i));
//            //let fin_bim = get_fname_bim(fin, Some(chrom_i));
//            //let num_snv_chrom = match text::compute_num_line(&fin_bim) {
//            let num_snv_chrom = match num_snv_chrom {
//                Some(v) => v,
//                None => 0,
//            };
//            num_snv += num_snv_chrom;
//        }
//
//        if num_snv == 0 {
//            return None;
//        } else {
//            return Some(num_snv);
//        }
//    }
//}

// use GenotFile::check_valid_fin()
//// TODO: shoudld use `check_open_file()`
//// here, should be check_exist
///// Check ...
///// If chrom is split,
/////     fin_chrom.bim, .fam, .bed exists
///// allow some of chorom not exist
//pub fn check_valid_fin(fin: &Path, gfmt: GenotFormat) {
//    if !judge_split_chrom(fin) {
//        // check three files exist
//        textfile::check_exist_file(&fname_plinks_sample(fin, gfmt, None));
//        textfile::check_exist_file(&fname_plinks_snv(fin, gfmt, None));
//        textfile::check_exist_file(&fname_plinks_genot(fin, gfmt, None));
//    } else {
//        for chrom_i in Chrom::variants().iter() {
//            let mut is_exists = vec![false, false, false];
//            is_exists[0] = textfile::exist_file(&fname_plinks_sample(fin, gfmt, Some(chrom_i)));
//            is_exists[1] = textfile::exist_file(&fname_plinks_snv(fin, gfmt, Some(chrom_i)));
//            is_exists[2] = textfile::exist_file(&fname_plinks_genot(fin, gfmt, Some(chrom_i)));
//
//            // check if all elements are the same
//            if !is_exists.iter().all(|v| *v == is_exists[0]) {
//                if !is_exists[0] {
//                    panic!(
//                        "Fam file does not exist: {:?}",
//                        fname_plinks_sample(fin, gfmt, Some(chrom_i))
//                    );
//                    //panic!("Fam file does not exist even though ~ file exists.");
//                } else if !is_exists[1] {
//                    panic!(
//                        "Bim file does not exist: {:?}",
//                        fname_plinks_snv(fin, gfmt, Some(chrom_i))
//                    );
//                } else if !is_exists[2] {
//                    panic!(
//                        "Bed file does not exist: {:?}",
//                        fname_plinks_genot(fin, gfmt, Some(chrom_i))
//                    );
//                }
//            }
//        }
//    }
//}

//// move to mod plink
////  TODO: should unwrap here?
///// return bed_size if valid, error otherwise
////pub fn check_valid_bed(fin: &str, n: usize, m: usize) -> Result<usize, String> {
//pub fn check_valid_bed(
//    //fin: &Path,
//    //gfmt: GenotFormat,
//    fin_genot: &GenotFile,
//    chrom: Option<&Chrom>,
//    m: usize,
//    n: usize,
//) -> Result<usize, Box<dyn Error>> {
//    let fin_bed = fin_genot.genotype_file(chrom);
//    //let fin_bed = fname_plinks_genot(fin, gfmt, chrom);
//    // check if open
//    let mut reader = File::open(fin_bed)?;
//
//    // check if size is correct
//    let f_end: usize = reader.seek(SeekFrom::End(0)).unwrap() as usize;
//    log::debug!("file end {}", f_end);
//    let bed_size = plink::calculate_bed_size(m, n);
//    if f_end != bed_size {
//        return Err(format!(
//            "File size of .bed is wrong: {} vs correct {}.",
//            f_end, bed_size
//        )
//        .into());
//    }
//
//    // check if the first 3 bytes are correct.
//    reader.seek(SeekFrom::Start(0)).unwrap();
//    let mut buf: Vec<u8> = Vec::with_capacity(n);
//    unsafe {
//        buf.set_len(3);
//    }
//    reader.read_exact(&mut buf).unwrap();
//    //log::debug!("{:?}", buf);
//    if buf != vec![0x6cu8, 0x1b, 0x01] {
//        return Err("Magic number of .bed file is wrong.".into());
//    }
//    Ok(bed_size)
//}

fn load_fam(fin_genot: &GenotFile) -> Vec<u8> {
    match fin_genot {
        GenotFile::Plink1(_) => {
            // TODO: all right? -> panic should be in load_ys_buf_option()
            panic!("Do not use --phe-name for plink1.");
        }
        GenotFile::Plink2(_) | GenotFile::Plink2Vzs(_) => {
            let fin_fam = fin_genot.sample_file_exist_chrom().unwrap();
            textfile::read_file_to_end(&fin_fam, None).unwrap()
        }
    }
}

//fn load_fam(fin: &Path, gfmt: GenotFormat) -> Vec<u8> {
//    match gfmt {
//        GenotFormat::Plink1 => {
//            panic!("Do not use --phe-name for plink1.");
//        }
//        GenotFormat::Plink2 | GenotFormat::Plink2Vzs => {
//            let fin_fam = fname_fam_exist_chrom(fin, gfmt).unwrap();
//            textfile::read_file_to_end(&fin_fam, None).unwrap()
//            //&phe_buf_v[..]
//            //phe_buf_v
//        }
//    }
//}

/// If phe_buf = None, phe_name = None, gfmt != plink1, then return None
// TODO: what if fin.fam includes phe=9?
// -> let user exclude samples with --fin-sample
pub fn load_ys_buf_file(
    fin_genot: &GenotFile,
    phe_buf: Option<&[u8]>,
    phe_name: Option<&str>,
    sample_id_to_n: &HashMap<String, usize>,
) -> Option<Vec<bool>> {
    let phe_buf_v: Vec<u8>;
    // TODO: or loading phe and samples_id separately
    let valss = if let Some(phe_name) = phe_name {
        let phe_buf = match phe_buf {
            None => {
                log::debug!("load_ys_buf(); phe_name exists and fphe does not exist.");
                // phe is in .psam
                phe_buf_v = load_fam(fin_genot);
                &phe_buf_v[..]
            }
            Some(x) => {
                // if phe_name is in phe_buf, use the col
                // elif phe in .psam, use the col
                if textfile::coli_of_header_buf(x, phe_name).is_some() {
                    log::debug!(
                        "load_ys_buf(); phe_name exists and fphe exists and phe_name is in column."
                    );
                    x
                } else {
                    log::debug!("load_ys_buf(); phe_name exists and fphe exists but phe_name is not in column.");
                    phe_buf_v = load_fam(fin_genot);
                    //phe_buf_v = load_fam(fin, gfmt);
                    if textfile::coli_of_header_buf(&phe_buf_v[..], phe_name).is_some() {
                        &phe_buf_v[..]
                    } else {
                        panic!("phe_name is not in fin-phe or psam file.");
                    }
                }
            }
        };
        //TODO
        let col_iid = 0;
        let col_phe = textfile::coli_of_header_buf(phe_buf, phe_name)
            .expect("phe_name is not in fin-phe or psam.");
        textfile::load_table_cols_buf(&phe_buf[..], &[col_iid, col_phe], true)
    } else {
        match fin_genot {
            GenotFile::Plink1(_) => {
                //TODO
                let col_iid = 1;
                let col_y = 5;
                let fin_fam = fin_genot.sample_file_exist_chrom().unwrap();
                //let fin_fam = fname_fam_exist_chrom(fin, gfmt).unwrap();
                textfile::load_table_cols(&fin_fam, &[col_iid, col_y], false)
            }
            GenotFile::Plink2(_) | GenotFile::Plink2Vzs(_) => {
                // TODO: Why?? -> for score??
                return None;
                //panic!("Use fin_phe for plink2.")
            }
        }
    };

    Some(load_ys_buf_vals(valss, sample_id_to_n))

    //let vals = sample::vals_align_id(&valss[1], &valss[0], sample_id_to_n);

    //log::debug!("vals[0]: {}", vals[0]);

    //let uniq = vec::uniq_string(&vals);
    //let code_type = if uniq == HashSet::from_iter([String::from("0"), String::from("1")]) {
    //    "01"
    //} else if uniq == HashSet::from_iter([String::from("1"), String::from("2")]) {
    //    "12"
    //} else {
    //    panic!("Unknown coding {:?}", uniq);
    //};

    //fn decode_phe(val: &str, code_type: &str) -> bool {
    //    let code: u8 = (*val).parse::<u8>().unwrap();
    //    if code_type == "01" {
    //        match code {
    //            0 | 1 => code != 0,
    //            z => panic!("Unknown phenotype included: {}.", z),
    //        }
    //    } else if code_type == "12" {
    //        match code {
    //            1 | 2 => (code - 1) != 0,
    //            z => panic!("Unknown phenotype included: {}.", z),
    //        }
    //    } else {
    //        panic!("Unknown code_type {:?}", code_type);
    //    }
    //}

    //let ys = vals.iter().map(|x| decode_phe(x, code_type)).collect();

    //Some(ys)
}

fn load_ys_buf_vals(valss: Vec<Vec<String>>, sample_id_to_n: &HashMap<String, usize>) -> Vec<bool> {
    let vals = sample::vals_align_id(&valss[1], &valss[0], sample_id_to_n);

    log::debug!("vals[0]: {}", vals[0]);

    let uniq = vec::uniq_clone(&vals);
    // let code_type = if uniq == HashSet::from_iter([String::from("0"), String::from("1"), String::from("-9")]) {
    let code_type = if uniq == HashSet::from_iter([String::from("0"), String::from("1")]) {
        "01"
    } else if uniq == HashSet::from_iter([String::from("1"), String::from("2")]) {
        "12"
    } else {
        panic!("Unknown coding {:?}", uniq);
    };

    fn decode_phe(val: &str, code_type: &str) -> bool {
        let code: u8 = (*val).parse::<u8>().unwrap();
        if code_type == "01" {
            match code {
                0 | 1 => code != 0,
                z => panic!("Unknown phenotype included: {}.", z),
            }
        } else if code_type == "12" {
            match code {
                1 | 2 => (code - 1) != 0,
                z => panic!("Unknown phenotype included: {}.", z),
            }
        } else {
            panic!("Unknown code_type {:?}", code_type);
        }
    }

    let ys = vals.iter().map(|x| decode_phe(x, code_type)).collect();

    ys
    //Some(ys)
}

//pub fn load_ys_buf_unwrap(
//    fin_genot: &GenotFile,
//    //fin: &Path,
//    //gfmt: GenotFormat,
//    phe_buf: Option<&[u8]>,
//    phe_name: Option<&str>,
//    sample_id_to_n: &HashMap<String, usize>,
//    //use_samples: &[bool],
//) -> Vec<bool> {
//    load_ys_buf(fin_genot, phe_buf, phe_name, sample_id_to_n).unwrap()
//}

// option for missing chrom
// No multi-allelic in plink1
fn load_snvs_text_plink1(fin_genot: &GenotFile) -> Vec<SnvId> {
    let fin_bim = fin_genot.snv_file();
    let buf = textfile::read_file_to_end(&fin_bim, None);
    let buf = buf.unwrap();

    // original order is chrom, rs, None, pos, A1, A2
    // rs, chrom, pos, A1, A2
    let cols = [1usize, 0, 3, 4, 5];
    let vss = textfile::load_table_cols_buf(&buf[..], &cols, false);

    let mut snvs: Vec<SnvId> = vec![];
    //let vss = textfile::load_table_cols(&fin_bim, &cols, false);
    for vi in 0..vss[0].len() {
        snvs.push(SnvId::new(
            vss[0][vi].clone(),
            vss[1][vi].clone(),
            &vss[2][vi],
            vss[3][vi].clone(),
            vss[4][vi].clone(),
        ));
    }
    snvs
    //Some(snvs)
}

#[derive(Clone, Debug, Deserialize)]
struct SnvPlink2In {
    #[serde(alias = "ID", alias = "#ID")]
    id: String,
    #[serde(alias = "CHROM", alias = "#CHROM")]
    chrom: String,
    #[serde(alias = "POS", alias = "#POS")]
    pos: String,
    #[serde(alias = "REF")]
    refa: String,
    #[serde(alias = "ALT")]
    alta: String,
}

impl SnvPlink2In {
    pub fn new(id: String, chrom: String, pos: String, refa: String, alta: String) -> Self {
        if textfile::isin_string(",", &alta) {
            panic!("Cannot have multiple alt alleles.");
        }
        SnvPlink2In {
            id,
            chrom,
            pos,
            refa,
            alta,
        }
    }
}

#[derive(Clone, Debug, Deserialize)]
struct FreqPlinkIn {
    #[serde(alias = "ID", alias = "#ID")]
    id: String,
    #[serde(alias = "CHROM", alias = "#CHROM")]
    chrom: String,
    #[serde(alias = "POS", alias = "#POS")]
    pos: String,
    #[serde(alias = "REF")]
    refa: String,
    #[serde(alias = "ALT")]
    alta: String,
    #[serde(alias = "ALT_FREQS")]
    alt_freq: f64,
    #[serde(alias = "REF_FREQ")]
    ref_freq: Option<f64>,
}

// this will fail for a file without header but contain string "ID" as name, but should be fine
// reading to end for comment lines at the top of file
// -> in plink2, any header should start with '#'
pub fn has_bim_header_plink2(fin_bim: &Path, compress: Option<&str>) -> bool {
    let buf = textfile::read_file_to_end(fin_bim, compress).unwrap();
    buf[0] == b'#'
    //if buf[0] == b'#' {
    //    return true;
    //}
    //textfile::isin_header("ID", &buf[..])
}

pub fn header_line_n_bim_plink2(fin_bim: &Path, compress: Option<&str>) -> usize {
    let buf = textfile::read_file_to_end(fin_bim, compress).unwrap();

    let mut n_line = 0;
    for line in buf.lines() {
        if !line.unwrap().starts_with('#') {
            break;
        }
        n_line += 1;
    }
    n_line

    //let mut rdr = csv::ReaderBuilder::new()
    //    .delimiter(b'\t')
    //    .from_reader(&buf[..]);
    //let mut n = 0;
    //for result in rdr.records() {
    //    let record = result.unwrap();
    //    if record[0] == "ID" {
    //        return n;
    //    }
    //    n += 1;
    //}
    //n
}

fn load_snvs_tsv(
    fin_genot: &GenotFile,
    compress: Option<&str>,
    file_i: Option<usize>,
) -> (Vec<SnvId>, Vec<(usize, usize, usize)>) {
    //) -> Vec<SnvId> {
    let fin_bim = fin_genot.snv_file();

    let mut snvs_in: Vec<SnvPlink2In> = vec![];
    let buf = textfile::read_file_to_end(&fin_bim, compress);
    let buf = buf.unwrap();

    if has_bim_header_plink2(&fin_bim, compress) {
        // Exclude lines starting with '##'
        let cursor = io::Cursor::new(&buf[..]);
        let mut start_index = 0;

        for line in cursor.lines().filter_map(|l| l.ok().map(|l| l)) {
            if !line.starts_with("##") {
                break;
            }
            // Update start_index to track the next line's starting byte position
            start_index += line.len() + 1; // +1 for '\n'
        }

        let mut rdr = csv::ReaderBuilder::new()
            .delimiter(b'\t')
            .from_reader(&buf[start_index..]);
        //.from_reader(&buf[..]);

        for result in rdr.deserialize() {
            let record: SnvPlink2In =
                result.unwrap_or_else(|_| panic!("Error while reading: {:?}", fin_bim));
            snvs_in.push(record)
        }
    } else {
        //.has_headers(false)
        // then should have exactly same cols as SnvIn, so avoid using Struct
        let mut rdr = csv::ReaderBuilder::new()
            .has_headers(false)
            .delimiter(b'\t')
            .from_reader(&buf[..]);

        // original order is chrom, rs, None, pos, A1, A2
        // rs, chrom, pos, A1, A2
        let cols = [1usize, 0, 3, 4, 5];

        for result in rdr.records() {
            let record = result.unwrap_or_else(|_| panic!("Error while reading: {:?}", fin_bim));
            //println!("{:?}", record);
            snvs_in.push(SnvPlink2In::new(
                record[cols[0]].to_string(),
                record[cols[1]].to_string(),
                record[cols[2]].to_string(),
                record[cols[3]].to_string(),
                record[cols[4]].to_string(),
            ))
        }
    }

    let mut snvs: Vec<SnvId> = vec![];

    let mut file_snv_allele_idx: Vec<(usize, usize, usize)> = vec![];
    for (snv_i, record) in snvs_in.into_iter().enumerate() {
        let alt_alleles = record.alta.split(",").collect::<Vec<&str>>();

        for (alt_i, alt) in alt_alleles.iter().enumerate() {
            let snv = SnvId::new(
                record.id.clone(),
                record.chrom.clone(),
                &record.pos,
                alt.to_string(),
                record.refa.clone(),
            );
            snvs.push(snv);
            match file_i {
                Some(file_i) => file_snv_allele_idx.push((file_i, snv_i, alt_i)),
                None => file_snv_allele_idx.push((0, snv_i, alt_i)),
                //None => file_snv_allele_idx.push((usize::MAX, snv_i, alt_i)),
            }
            //file_snv_allele_idx.push((usize::MAX, snv_i, alt_i));
        }
    }

    return (snvs, file_snv_allele_idx);

    //for record in snvs_in {
    //    let snv = SnvId::new(
    //        record.id.clone(),
    //        &record.chrom,
    //        &record.pos,
    //        record.alta.clone(),
    //        record.refa.clone(),
    //    );
    //    snvs.push(snv);
    //}
    //snvs
    //Some(snvs)
}

// option for missing chrom
fn load_snvs_file(
    fin_genot: &GenotFile,
    file_i: Option<usize>,
) -> (Vec<SnvId>, Vec<(usize, usize, usize)>) {
    match fin_genot {
        GenotFile::Plink1(_) => (load_snvs_text_plink1(fin_genot), vec![]),
        GenotFile::Plink2(_) => load_snvs_tsv(fin_genot, None, file_i),
        GenotFile::Plink2Vzs(_) => load_snvs_tsv(fin_genot, Some("zst"), file_i),
    }
}

// option for missing chrom
//fn load_snvs_chrom(fin: &Path, gfmt: GenotFormat, chrom: Option<&Chrom>) -> Option<Vec<SnvId>> {
//    match gfmt {
//        GenotFormat::Plink1 => load_snvs_text_plink(fin, gfmt, chrom),
//        GenotFormat::Plink2 => load_snvs_tsv(fin, gfmt, chrom, None),
//        GenotFormat::Plink2Vzs => load_snvs_tsv(fin, gfmt, chrom, Some("zst")),
//    }
//}

//// TODO: why not use_snvs here?
//// since match snv_index is not simple?
//pub fn load_snvs(fins_genot: &GenotFiles) -> Vec<SnvId> {
//    let mut snvs: Vec<SnvId> = vec![];
//    //for chrom_i in Chrom::variants().iter() {
//    for fin_genot in fins_genot.0.iter() {
//        let v = load_snvs_file(fin_genot);
//        //let v = load_snvs_chrom(fin_genot, Some(chrom_i));
//        if let Some(mut snvs_i) = v {
//            snvs.append(&mut snvs_i);
//        }
//        // else continue
//    }
//    // all files returned None
//    if snvs.len() == 0 {
//        panic!("No snvs were loaded from fins_genot");
//    }
//    snvs
//
//    //if !fin_genot.judge_split_chrom() {
//    //    load_snvs_chrom(fin_genot, None)
//    //        .unwrap_or_else(|| panic!("Could not load snvs from fin: {:?}", fin_genot.file()))
//    //} else {
//    //    let mut snvs: Vec<SnvId> = vec![];
//    //    for chrom_i in Chrom::variants().iter() {
//    //        let v = load_snvs_chrom(fin_genot, Some(chrom_i));
//    //        if let Some(mut snvs_i) = v {
//    //            snvs.append(&mut snvs_i);
//    //        }
//    //        // else continue
//    //    }
//    //    // all chrom returned None
//    //    if snvs.len() == 0 {
//    //        panic!("No snvs were loaded from {:?}", fin_genot.file());
//    //    }
//    //    snvs
//    //}
//}

fn load_freq_tsv(
    freq_buf: &[u8],
    //freq_buf: Option<&[u8]>,
    //compress: Option<&str>,
) -> Snvs {
    let buf = freq_buf;
    //let fin_bim = fin_genot.snv_file(chrom);
    ////let fin_bim = fname_plinks_snv(fin, gfmt, chrom);
    ////println!("fin_bim {:?}", fin_bim);

    let mut snvs_in: Vec<FreqPlinkIn> = vec![];
    //let buf = textfile::read_file_to_end(&fin_bim, compress);
    //if buf.is_err() {
    //    return None;
    //}
    //let buf = buf.unwrap();

    //if has_bim_header_plink2(&fin_bim, compress) {
    let mut rdr = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .from_reader(&buf[..]);

    for result in rdr.deserialize() {
        let record: FreqPlinkIn = result.unwrap_or_else(|_| panic!("Error while reading freq"));
        snvs_in.push(record)
    }
    //} else {
    //    //.has_headers(false)
    //    // then should have exactly same cols as SnvIn, so avoid using Struct
    //    let mut rdr = csv::ReaderBuilder::new()
    //        .has_headers(false)
    //        .delimiter(b'\t')
    //        .from_reader(&buf[..]);

    //    // original order is chrom, rs, None, pos, A1, A2
    //    // rs, chrom, pos, A1, A2
    //    let cols = [1usize, 0, 3, 4, 5];

    //    for result in rdr.records() {
    //        let record = result.unwrap_or_else(|_| panic!("Error while reading: {:?}", fin_bim));
    //        println!("{:?}", record);
    //        snvs_in.push(SnvPlink2In::new(
    //            record[cols[0]].to_string(),
    //            record[cols[1]].to_string(),
    //            record[cols[2]].to_string(),
    //            record[cols[3]].to_string(),
    //            record[cols[4]].to_string(),
    //        ))
    //    }
    //}

    let mut snvs: Vec<SnvId> = vec![];
    let mut freqs_alt: Vec<f64> = vec![];
    let mut freqs_ref: Vec<Option<f64>> = vec![];
    // let mut freqs_ref: Vec<f64> = vec![];

    // TOFIX: multiple alt_freqs
    for record in snvs_in {
        let snv = SnvId::new(
            record.id.clone(),
            record.chrom.clone(),
            &record.pos,
            record.alta.clone(),
            record.refa.clone(),
        );
        snvs.push(snv);
        freqs_alt.push(record.alt_freq);
        freqs_ref.push(record.ref_freq);
    }
    // TOFIX: correct?
    let freqs_ref = if freqs_ref.iter().any(|x| x.is_none() || x.unwrap().is_nan()) {
        None
    } else {
        Some(freqs_ref.into_iter().map(|x| x.unwrap()).collect())
        // Some(freqs_ref)
    };

    Snvs::new_snv(snvs, Some(freqs_alt), freqs_ref)
    //Some(snvs)
}

//pub fn load_freq(freq_buf: Option<&[u8]>) -> Option<Snvs> {
pub fn load_freq(freq_buf: &[u8]) -> Snvs {
    // TODO: compressed ver.

    //match freq_buf {
    //    Some(x)=>load_freq(x),
    //    None=>None,
    //}
    load_freq_tsv(freq_buf)

    //let freq: Option<Snvs> = match freq_buf {
    //    Some(x) => x,
    //    None => None,
    //};

    // cov_name.map(|x| Covs::new_buf(phe_buf, fin_genot, x, &sample_id_to_n));

    //if !fin_genot.judge_split_chrom() {
    //    load_snvs_chrom(fin_genot, None)
    //        .unwrap_or_else(|| panic!("Could not load snvs from fin: {:?}", fin_genot.file()))
    //} else {
    //    let mut snvs: Vec<SnvId> = vec![];
    //    for chrom_i in Chrom::variants().iter() {
    //        let v = load_snvs_chrom(fin_genot, Some(chrom_i));
    //        if let Some(mut snvs_i) = v {
    //            snvs.append(&mut snvs_i);
    //        }
    //        // else continue
    //    }
    //    // all chrom returned None
    //    if snvs.len() == 0 {
    //        panic!("No snvs were loaded from {:?}", fin_genot.file());
    //    }
    //    snvs
    //}
}

//// TODO: why not use_snvs here?
//// since match snv_index is not simple?
//pub fn load_snvs(fin: &Path, gfmt: GenotFormat) -> Vec<SnvId> {
//    if !judge_split_chrom(fin) {
//        //load_snvs_chrom(fin, gfmt, None).unwrap()
//
//        load_snvs_chrom(fin, gfmt, None)
//            .unwrap_or_else(|| panic!("Could not load snvs from fin: {:?}", fin))
//    } else {
//        let mut snvs: Vec<SnvId> = vec![];
//        for chrom_i in Chrom::variants().iter() {
//            let v = load_snvs_chrom(fin, gfmt, Some(chrom_i));
//            if let Some(mut snvs_i) = v {
//                snvs.append(&mut snvs_i);
//            }
//            // else continue
//        }
//        // all chrom returned None
//        if snvs.len() == 0 {
//            panic!("No snvs were loaded from {:?}", fin);
//        }
//        snvs
//    }
//}

#[derive(Clone, Debug, Deserialize)]
struct SamplePlink2In {
    // ng
    //#[serde(rename = ("FID","FID"))]
    //#[serde(alias = "FID", alias = "#FID")]
    //fid: Option<String>,
    #[serde(alias = "ID", alias = "IID", alias = "#ID", alias = "#IID")]
    iid: String,
}

impl SamplePlink2In {
    pub fn new(iid: String) -> Self {
        SamplePlink2In {
            //fid: Some(fid),
            iid,
        }
    }
}

// this will fail for a file without header but contain string "ID" as name, but should be fine
pub fn has_fam_header_plink2(fin_fam: &Path) -> bool {
    let buf = textfile::read_file_to_end(fin_fam, None).unwrap();
    textfile::isin_header("ID", &buf[..])
}

// TODO: iid col number
// ex. fid, iid cols = [0usize, 1];
fn load_samples_id_tsv(fin_genot: &GenotFile, use_samples: Option<&[bool]>) -> Vec<String> {
    let fin_fam = fin_genot.sample_file_exist_chrom().unwrap();
    let mut samples_in: Vec<SamplePlink2In> = vec![];
    if has_fam_header_plink2(&fin_fam) {
        let mut rdr = csv::ReaderBuilder::new()
            .delimiter(b'\t')
            .from_path(&fin_fam)
            .unwrap();
        // cannot load properly
        //.quote(b'#')

        for result in rdr.deserialize() {
            let record: SamplePlink2In =
                result.unwrap_or_else(|_| panic!("Error while reading: {:?}", fin_fam));
            samples_in.push(record)
        }
    } else {
        //.has_headers(false)
        // then should have exactly same cols as SnvIn, so avoid using Struct
        let mut rdr = csv::ReaderBuilder::new()
            .has_headers(false)
            .delimiter(b'\t')
            .from_path(&fin_fam)
            .unwrap();

        for result in rdr.records() {
            let record = result.unwrap_or_else(|_| panic!("Error while reading: {:?}", fin_fam));
            //println!("{:?}", record);
            // the first line sould be FID or IID and both are fine
            samples_in.push(SamplePlink2In::new(record[0].to_string()))
        }
    }

    //let n = use_samples.iter().filter(|&v| *v).count();
    //let mut samples_id: Vec<(String, String)> = Vec::with_capacity(n);

    let mut samples_id: Vec<String> = vec![];

    if use_samples.is_some() {
        for (record, b) in samples_in.into_iter().zip(use_samples.unwrap().iter()) {
            if *b {
                samples_id.push(record.iid);
            }
        }
    } else {
        for record in samples_in.into_iter() {
            samples_id.push(record.iid);
        }
    }

    samples_id
}

fn load_samples_id_text_plink1(fin_genot: &GenotFile, use_samples: Option<&[bool]>) -> Vec<String> {
    let fin_fam = fin_genot.sample_file_exist_chrom().unwrap();

    // fid, iid
    //let cols = [0usize, 1];
    // assume fid=iid
    let col_iid = 1;
    let cols = [col_iid];
    let samples_in: Vec<Vec<String>> = textfile::load_table_cols(&fin_fam, &cols, false);
    if use_samples.is_some() {
        assert_eq!(samples_in[0].len(), use_samples.unwrap().len());
    }

    //let n = use_samples.iter().filter(|&v| *v).count();
    //let mut sample_id_to_n: HashMap<String, usize> = HashMap::with_capacity(n);
    //let mut samples_id: Vec<(String, String)> = Vec::with_capacity(n);

    let mut samples_id: Vec<String> = vec![];

    if use_samples.is_some() {
        for (n_in_i, v) in use_samples.unwrap().iter().enumerate() {
            if *v {
                samples_id.push(samples_in[0][n_in_i].clone());
            }
        }
    } else {
        for n_in_i in 0..samples_in[0].len() {
            samples_id.push(samples_in[0][n_in_i].clone());
        }
    }

    samples_id
}

pub fn load_samples_id_file(fin_genot: &GenotFile, use_samples: Option<&[bool]>) -> Vec<String> {
    match fin_genot {
        GenotFile::Plink1(_) => load_samples_id_text_plink1(fin_genot, use_samples),
        GenotFile::Plink2(_) | GenotFile::Plink2Vzs(_) => {
            load_samples_id_tsv(fin_genot, use_samples)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_load_samples_id_plink1() {
        let fin = PathBuf::from("../../test/data/toy1/genot");
        let fin_genot = GenotFile::Plink1(fin);

        let v = load_samples_id_file(&fin_genot, None);
        assert_eq!(v[0], "1");
    }

    #[test]
    fn test_load_samples_id_plink2() {
        let fin = PathBuf::from("../../test/data/toy3/genot");
        let fin_genot = GenotFile::Plink2(fin);

        let v = load_samples_id_file(&fin_genot, None);
        assert_eq!(v[0], "1");
    }

    /*     #[test]
    fn test_load_ys() {
        let fin = String::from("../../test/data/toy1/genot");

        let n = compute_num_sample(&fin).unwrap();
        let (_, use_samples) = make_use_samples(None, &fin, n);

        let ys = load_ys_b8(&fin, n, &use_samples);
        let ys_expect = vec![true, true, false, false, false];
        //let ys_expect = vec![false, true, false, false, true]; -> this is sex

        for (i, exp) in ys_expect.iter().enumerate() {
            assert_eq!(bit::bget(&ys, i), *exp);
        }
    } */

    #[test]
    fn test_compute_num_sample_file() {
        let fin = PathBuf::from("../../test/data/toy3/genot");

        let fin_genot = GenotFile::Plink1(fin.clone());
        let m = compute_num_sample_file(&fin_genot);
        assert_eq!(m, 10);

        let fin_genot = GenotFile::Plink2(fin.clone());
        let m = compute_num_sample_file(&fin_genot);
        assert_eq!(m, 10);

        let fin_genot = GenotFile::Plink2Vzs(fin.clone());
        let m = compute_num_sample_file(&fin_genot);
        assert_eq!(m, 10);
    }

    #[test]
    fn test_compute_num_snv_file_toy3() {
        let fin = PathBuf::from("../../test/data/toy3/genot");

        let fin_genot = GenotFile::Plink1(fin.clone());
        let m = compute_num_snv_file(&fin_genot, None);
        assert_eq!(m, 3);

        let fin_genot = GenotFile::Plink2(fin.clone());
        let m = compute_num_snv_file(&fin_genot, None);
        assert_eq!(m, 3);

        let fin_genot = GenotFile::Plink2Vzs(fin.clone());
        let m = compute_num_snv_file(&fin_genot, None);
        assert_eq!(m, 3);
    }

    #[test]
    fn test_compute_num_snv_file_toy8() {
        let fin = PathBuf::from("../../test/data/toy8/genot.ref");

        let fin_genot = GenotFile::Plink2(fin.clone());
        let m = compute_num_snv_file(&fin_genot, None);
        assert_eq!(m, 15);

        let fin_genot = GenotFile::Plink2Vzs(fin.clone());
        let m = compute_num_snv_file(&fin_genot, None);
        assert_eq!(m, 15);
    }

    #[test]
    fn test_compute_num_snv_toy7() {
        let fins = vec![
            PathBuf::from("../../test/data/toy7/genot_1.ref"),
            PathBuf::from("../../test/data/toy7/genot_2.ref"),
        ];

        let fin_genot = GenotFiles::new_plink2(fins.clone());
        let m = compute_num_snv(&fin_genot);
        assert_eq!(m, 10);

        let fin_genot = GenotFiles::new_plink2vzs(fins.clone());
        let m = compute_num_snv(&fin_genot);
        assert_eq!(m, 10);
    }

    #[test]
    fn test_load_ys_buf_vals_01() {
        let valss = vec![
            "C,A,B,D"
                .split(",")
                .map(|x| x.to_string())
                .collect::<Vec<String>>(),
            "0,1,1,0"
                .split(",")
                .map(|x| x.to_string())
                .collect::<Vec<String>>(),
        ];

        let sample_id_to_n = HashMap::from([("B".to_string(), 0), ("C".to_string(), 1)]);
        //let mut sample_id_to_n = HashMap::new();
        //sample_id_to_n.insert("B".to_string(), 0);
        //sample_id_to_n.insert("C".to_string(), 1);

        let ys = load_ys_buf_vals(valss, &sample_id_to_n);

        // phe of [B, C]
        let ys_ans = vec![true, false];

        assert_eq!(ys, ys_ans);
    }

    #[test]
    fn test_load_ys_buf_vals_12() {
        let valss = vec![
            "C,A,B,D"
                .split(",")
                .map(|x| x.to_string())
                .collect::<Vec<String>>(),
            "1,2,2,1"
                .split(",")
                .map(|x| x.to_string())
                .collect::<Vec<String>>(),
        ];

        let sample_id_to_n = HashMap::from([("B".to_string(), 0), ("C".to_string(), 1)]);
        //let mut sample_id_to_n = HashMap::new();
        //sample_id_to_n.insert("B".to_string(), 0);
        //sample_id_to_n.insert("C".to_string(), 1);

        let ys = load_ys_buf_vals(valss, &sample_id_to_n);

        // phe of [B, C]
        let ys_ans = vec![true, false];

        assert_eq!(ys, ys_ans);
    }
}
