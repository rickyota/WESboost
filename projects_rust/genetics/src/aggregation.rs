use crate::{genot::prelude::*, textfile};
// use genetics::genot::BaseGenot;
use std::{io::Write, path::Path};

use crate::Dataset;

pub fn group_freq(fgroup_snv_freq: &Path, dataset: &Dataset) {
    log::debug!("genot: {:?}", dataset.genot().m());
    log::debug!("snvs: {:?}", dataset.snvs().snv_ids());
    log::debug!("snv_ids group: {:?}", dataset.snv_ids_group());

    let (_, g_group) = dataset.genot_split_group();

    let freqs = g_group
        .iter_snv()
        .map(|g_snv| g_snv.maf_group())
        .collect::<Vec<f64>>();
    log::debug!("freqs: {:?}", freqs);

    let str_header = "groupid\tfreq".to_owned();
    let strings = dataset
        .snv_ids_group()
        .iter()
        .zip(freqs.iter())
        .map(|(snv_id, freq)| format!("{}\t{}", snv_id.id(), freq))
        .collect::<Vec<String>>();
    let str_out = format!("{}\n{}", str_header, strings.join("\n"));

    let mut writer = textfile::bufwriter(&fgroup_snv_freq);
    writer.write_all(&str_out.as_bytes()).unwrap();
}
