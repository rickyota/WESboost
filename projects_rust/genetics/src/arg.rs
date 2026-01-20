use clap::ValueEnum;
use std::path::PathBuf;

use crate::GenotFiles;

#[derive(Copy, Clone, PartialEq, Eq, Debug, ValueEnum)]
pub enum GenotFormatArg {
    Plink1,
    Plink2,
    Plink2Vzs,
}

impl GenotFormatArg {
    pub fn to_genot_file(self, fin: Vec<PathBuf>) -> GenotFiles {
        match self {
            GenotFormatArg::Plink1 => GenotFiles::new_plink1(fin),
            GenotFormatArg::Plink2 => GenotFiles::new_plink2(fin),
            GenotFormatArg::Plink2Vzs => GenotFiles::new_plink2vzs(fin),
        }
    }
}
