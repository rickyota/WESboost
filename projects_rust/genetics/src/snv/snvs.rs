//! Snvs
//! It is important that all vector values can be extract without constructing. ex. mafs: Vec<f64> not Vec<Info>
//! Must use fn to access data so that it is easy to use Trait

//use crate::{genot_io, snv, vec};
//use crate::{GenotFormat, SnvId};
//use std::path::Path;
use crate::{genot_io, snv, vec, GenotFiles, SnvId};

use super::group_index::AggId;

#[derive(Clone, Default)]
pub struct Snvs {
    // Regard both group and snv as SnvId
    snv_ids: Vec<SnvId>,
    // freq of a1
    // For freqs_snvs input, freq of alt
    freqs_a1: Option<Vec<f64>>,
    // maf: Option<Vec<f64>>,
    // freq of a2
    // For freqs_snvs input, freq of ref
    freqs_a2: Option<Vec<f64>>,
    // mafs_outside: Option<> // for missing value and pass maf file as argument
    // snvs_n+group_n = len(snv_ids)
    // not always len(snv_ids) since sometimes snv_ids are empty // <- when!?!?
    // -> old?; now, it seems that snvs_n is the number of snvs in snv_ids
    // -> say, using freqs_a1 but not using snv_ids??
    snvs_n: usize,
    group_n: usize,
    //qc: Option<SnvsQc>, //later
    agg_ids: Option<Vec<AggId>>,
    agg_to_m: Option<Vec<Vec<usize>>>,
}

impl Snvs {
    pub fn snvs_snv_n(&self) -> usize {
        self.snvs_n
    }

    // TODO: rename to snvs_n()
    pub fn snvs_both_n(&self) -> usize {
        self.snvs_n + self.group_n
    }

    // pub fn snvs_n(&self) -> usize {
    //     self.snvs_n
    // }

    pub fn group_n(&self) -> usize {
        self.group_n
    }

    pub fn snvs(&self) -> &Self {
        &self
    }
    pub fn freqs_a1(&self) -> Option<&Vec<f64>> {
        self.freqs_a1.as_ref()
    }
    pub fn freqs_a2(&self) -> Option<&Vec<f64>> {
        self.freqs_a2.as_ref()
    }
    pub fn snv_ids(&self) -> &[SnvId] {
        &self.snv_ids
    }

    pub fn snv_ids_mut(&mut self) -> &mut [SnvId] {
        &mut self.snv_ids
    }

    pub fn agg_ids(&self) -> Option<&Vec<AggId>> {
        self.agg_ids.as_ref()
    }

    pub fn agg_to_m(&self) -> Option<&Vec<Vec<usize>>> {
        self.agg_to_m.as_ref()
    }

    pub fn new_plink_use_snvs(
        fins_genot: &GenotFiles,
        use_snvs: Option<&[bool]>,
        m: Option<usize>, // avoid loading again
    ) -> Self {
        let m = match m {
            Some(x) => x,
            None => match use_snvs {
                Some(x) => vec::count_true(x),
                None => genot_io::compute_num_snv(&fins_genot),
            },
        };

        //let m_in: usize = genot_io::compute_num_snv(fin, gfmt).unwrap();

        let (snvs_in, _) = genot_io::load_snvs(fins_genot);
        let snv_ids = match use_snvs {
            Some(x) => snv::extract_snvs_consume(snvs_in, x, m),
            None => snvs_in,
        };
        //let use_snvs = vec![true; m_in];
        //let snv_indexs = snv::extract_snvs_consume(snvs_in, &use_snvs, m_in);
        let snvs = Self::new_from_snv_ids(snv_ids);
        snvs
    }

    pub fn new_plink_use_snvs_and_group(
        fins_genot: &GenotFiles,
        use_snvs_or_agg: Option<&[bool]>,
        m: Option<usize>, // avoid loading again
        agg_snv_buf: Option<&[u8]>,
        group_snv_buf: Option<&[u8]>,
    ) -> Self {
        // log::debug!("use_snvs_or_agg: {:?}", use_snvs_or_agg);
        log::debug!("m: {:?}", m);

        // TOFIX: add 'is_snv' and 'is_agg' to SnvId
        let snvs = Self::new_plink_use_snvs(fins_genot, use_snvs_or_agg, m);
        log::debug!("snvs m: {:?}", snvs.snvs_snv_n());

        // TOFIX: exlcude unloaded snvs in group
        // TOFIX: update group_ids to SnvId in Snvs
        // TOFIX: return group_to_m
        let mut group_snv_ids = snv::load_group_snvs_buf(group_snv_buf);
        let m_group = group_snv_ids.len();
        log::debug!("snvs m_group: {:?}", m_group);

        let mut snv_ids = snvs.snv_ids().to_vec();
        snv_ids.append(&mut group_snv_ids);
        // TODO: clean
        let mut snvs = Self::new_from_snv_ids_group(snv_ids, m_group);
        // let mut snvs = Self::new_from_snv_ids(snv_ids);
        // snvs.set_group_n(m_group);

        assert_eq!(snvs.snvs_snv_n() + snvs.group_n(), snvs.snv_ids().len());

        let (agg_snv, agg_to_m) = snv::make_agg_to_m_buf(agg_snv_buf, snvs.snv_ids());
        snvs.set_agg_ids(Some(agg_snv));
        snvs.set_agg_to_m(Some(agg_to_m));

        snvs
    }

    pub fn set_group_n(&mut self, group_n: usize) {
        self.group_n = group_n;
    }

    pub fn set_agg_ids(&mut self, agg_ids: Option<Vec<AggId>>) {
        self.agg_ids = agg_ids;
    }

    pub fn set_agg_to_m(&mut self, agg_to_m: Option<Vec<Vec<usize>>>) {
        self.agg_to_m = agg_to_m;
    }

    pub fn new_from_snv_ids(snv_ids: Vec<SnvId>) -> Self {
        Self::new_snv(snv_ids, None, None)
    }

    pub fn new_from_snv_ids_group(snv_ids: Vec<SnvId>, m_group: usize) -> Self {
        // Self::new_snv(snv_ids, None, None)

        let snvs_n = snv_ids.len() - m_group;
        Self {
            snv_ids,
            freqs_a1: None,
            freqs_a2: None,
            snvs_n,
            group_n: m_group,
            agg_ids: None,
            agg_to_m: None,
        }
    }

    // tmp
    pub fn new_empty() -> Self {
        Self::new_snv(vec![], None, None)
        //Self {
        //    snv_indexs: vec![],
        //    mafs: None,
        //    snvs_n: 0,
        //}
    }

    // TODO: new_check() check len
    pub fn new_snv(
        snv_ids: Vec<SnvId>,
        freqs_a1: Option<Vec<f64>>,
        freqs_a2: Option<Vec<f64>>,
    ) -> Self {
        let snvs_n = snv_ids.len();
        Self {
            snv_ids,
            freqs_a1,
            freqs_a2,
            snvs_n,
            group_n: 0,
            agg_ids: None,
            agg_to_m: None,
        }
    }

    // or set_maf(&self, mafs)->Self{}
    // -> this way makes difficult to set_maf if in Dataset{Snvs}
    pub fn set_maf(&mut self, mafs: Vec<f64>) {
        self.freqs_a1 = Some(mafs);
    }

    // for snv only, not group
    pub fn extract_snvs(self, use_snvs: &[bool]) -> Self {
        let Self {
            snv_ids,
            freqs_a1,
            freqs_a2,
            snvs_n,
            group_n,
            agg_ids,
            agg_to_m,
        } = self;

        assert_eq!(snvs_n, use_snvs.len());

        let snv_ids_use: Vec<SnvId> = snv_ids
            .iter()
            .zip(use_snvs.iter())
            .filter(|(_, b)| **b)
            .map(|(x, _)| x.clone())
            .collect();

        let freqs_a1_use: Option<Vec<f64>> = freqs_a1.map(|x| {
            x.iter()
                .zip(use_snvs.iter())
                .filter(|(_, b)| **b)
                .map(|(x, _)| x.clone())
                .collect()
        });

        let freqs_a2_use: Option<Vec<f64>> = freqs_a2.map(|x| {
            x.iter()
                .zip(use_snvs.iter())
                .filter(|(_, b)| **b)
                .map(|(x, _)| x.clone())
                .collect()
        });

        let n = snv_ids_use.len();

        Self {
            snv_ids: snv_ids_use,
            freqs_a1: freqs_a1_use,
            freqs_a2: freqs_a2_use,
            snvs_n: n,
            group_n,
            agg_ids,
            agg_to_m,
        }
    }

    // pub fn use_snvs_chrom(&self, chrom: &Chrom) -> Vec<bool> {
    pub fn use_snvs_chrom(&self, chrom: &str) -> Vec<bool> {
        self.snv_ids().iter().map(|x| x.chrom() == chrom).collect()
    }

    //pub fn extract_chrom_indexs(&self, chrom: &Chrom) -> Vec<usize> {
    //    self.snv_ids()
    //        .iter()
    //        .enumerate()
    //        .filter(|(_, x)| x.chrom() == chrom)
    //        .map(|(i, _)| i)
    //        .collect()
    //}

    pub fn positions(&self) -> Vec<usize> {
        self.snv_ids().iter().map(|x| x.pos()).collect()
    }

    // use extrac_chrom_indexs() and positions()
    //pub fn positions(&self, chrom: &Chrom) -> Vec<usize> {
    //    self.snv_ids()
    //        .iter()
    //        .filter(|x| x.chrom() == chrom)
    //        .map(|x| x.pos())
    //        .collect()
    //}
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_extract_snvs() {
        let snv_ids = vec![SnvId::default(); 3];
        let snvs = Snvs::new_from_snv_ids(snv_ids);
        let use_snvs = [true, false, true];
        let snvs_use = snvs.extract_snvs(&use_snvs);

        //let snv_use_exp = Snvs::new_from_snv_index(vec![SnvId::default(); 2]);

        assert_eq!(snvs_use.snv_ids().len(), 2);
    }
}
