use super::load;
use crate::dataset::FillMissingScore;
use crate::genot::prelude::*;
use crate::snv;
use crate::FillMissing;
use crate::FillMissingGroup;
use crate::Genot;
use crate::GenotFiles;
use crate::Snvs;
use crate::WgtKind;
//use crate::GenotFormat;
use crate::SnvId;
use crate::WgtTrait;
use crate::Wgts;
use crate::{genot_io, vec};
use rayon::prelude::*;
use std::collections::HashMap;
use std::hash::Hash;
//use std::path::Path;

// TODO: merge to load_genotypes_for_score_multiwgs
// mainly for boosting?
// create allow_nonexist_snv=false ver. in boosting
// refer to generate_genot_for_score_multiwgts()
// pub fn generate_genot_for_score<W: WgtTrait>(
pub fn generate_genot_for_score_boosting<W: WgtTrait>(
    fins_genot: &GenotFiles,
    wgts: &mut [W],
    //wgts: &mut [Wgt],
    n: usize,
    use_samples: Option<&[bool]>,
    //fill_missing: bool,
    // deprecated: for backward compatibility
    // fill_missing_mode_in_dataset: bool,
    allow_nonexist_snv: bool,
    use_snv_pos: bool,
    fill_missing_score: Option<FillMissingScore>,
    fill_missing_group: Option<FillMissingGroup>,
    mem: Option<usize>,
) -> Genot {
    let (snvs_in, file_snv_allele_idx_in) = genot_io::load_snvs(fins_genot);

    // For rough filtering. Snvs with flipped alleles are contained. They will not be used since wgts knows it.
    let use_snvs = create_wgt_to_genotype_index(&snvs_in, wgts, use_snv_pos);

    let m = vec::count_true(&use_snvs);
    // if m == 0 {
    //     panic!("No variants to be loaded from genotype file.")
    // }

    // length is m_in
    // ?? -> m
    let snvs_use: Vec<SnvId> = vec::extract_if_into_iter(snvs_in.clone(), &use_snvs);
    let snvs_use = Snvs::new_from_snv_ids(snvs_use);

    //let allow_nonexist_snv = false;
    // To check unconsistent alleles
    let is_reversed_hashmap =
        set_wgts_index_reversed_hashmap(wgts, &snvs_use, allow_nonexist_snv, use_snv_pos);

    //     // this happens for group only
    // if is_reversed_hashmap.len() == 0 {
    //     panic!("No variants to be loaded after matching alelles.")
    // }

    // length is m
    let is_reversed =
        is_reversed_to_vec_bool(&is_reversed_hashmap, snvs_use.snv_ids(), use_snv_pos);
    log::debug!("Reversed snvs: {}", vec::count_true(&is_reversed));

    // extract group from wgt
    let snvs_group: Vec<SnvId> = wgts
        .iter()
        .filter(|wgt| wgt.kind().is_snv_single())
        .filter(|wgt| wgt.kind().snv_id().is_group())
        .map(|wgt| wgt.kind().snv_id().clone())
        .collect();
    let snvs_group = dedup_snvs_group(&snvs_group);

    let group_to_m_in = if snvs_group.len() == 0 {
        None
    } else {
        Some(snv::make_group_to_m_in_vec(&snvs_group, &snvs_in))
    };

    // let mut group_mi = m;
    let group_id_to_mi: HashMap<&str, usize> = snvs_group
        .iter()
        .enumerate()
        .map(|(group_id_i, group_id)| (group_id.id(), m + group_id_i))
        .collect();

    for wgt in wgts
        .iter_mut()
        .filter(|wgt| wgt.kind().is_snv_single())
        .filter(|wgt| wgt.kind().snv_id().is_group())
    {
        let mi = group_id_to_mi[wgt.kind().snv_id().id()];
        wgt.set_snv_index_check(Some(mi));
    }

    log::info!("Temporary: Skip index reversed hashmap for group");

    if m == 0 && snvs_group.len() == 0 {
        panic!("No variants to be loaded from genotype file.")
    }

    // TOFIX: should not fill_missing_mode in score
    let mut g = load::generate_genot(
        fins_genot,
        m,
        n,
        Some(&use_snvs),
        &file_snv_allele_idx_in,
        group_to_m_in,
        // Some(group_to_m_in),
        // None,
        use_samples,
        //false,
        //fill_missing_mode_in_dataset,
        // FillMissing:Ref if FillMissingScore::Ref
        // None otherwise
        FillMissing::from_fill_missing_score_generate_genot(fill_missing_score),
        fill_missing_group,
        mem,
    );

    reverse_genotypes_rayon(&mut g.as_genot_mut(), &is_reversed);

    g
}

fn integrate_use_snvs(uses_snvs: &[Vec<bool>]) -> Vec<bool> {
    let mut use_snvs = uses_snvs[0].clone();
    for use_snv_ in uses_snvs.iter() {
        vec::or_bool_vec_mut(&mut use_snvs, use_snv_);
    }
    use_snvs
}

fn dedup_snvs_group(snvs_group: &[SnvId]) -> Vec<SnvId> {
    let mut hash_id = HashMap::new();

    for snv in snvs_group.iter() {
        hash_id.insert(snv.id(), snv.clone());
    }

    hash_id.into_values().collect()
}

//fn integrate_is_reversed(iss_reversed: &[Vec<bool>]) -> Vec<bool> {
//    let mut is_reversed = iss_reversed[0].clone();
//    for is_reversed_ in iss_reversed.iter() {
//        vec::or_bool_vec_mut(&mut is_reversed, is_reversed_);
//    }
//    is_reversed
//}

//fn check_consistent_is_reversed(
//    is_reversed: &[bool],
//    iss_reversed: &[Vec<bool>],
//    useds_snvs: &[Vec<bool>],
//) {
//    for (is_reversed_, used_snv_) in iss_reversed.iter().zip(useds_snvs.iter()) {
//        for (snv_i, &used_snv_val) in used_snv_.iter().enumerate() {
//            if used_snv_val {
//                if is_reversed_[snv_i] != is_reversed[snv_i] {
//                    panic!("Alleles are reversed within fwgt in dwgt.");
//                }
//            }
//        }
//    }
//}

fn integrate_is_reversed_hashmap_check(
    iss_reversed_hashmap: Vec<HashMap<String, bool>>,
    //iss_reversed_hashmap: &[HashMap<String, bool>],
) -> HashMap<String, bool> {
    let mut is_reversed: HashMap<String, bool> = HashMap::new();
    //let mut is_reversed = iss_reversed[0].clone();
    for is_reversed_wgt in iss_reversed_hashmap.into_iter() {
        // check no conflict
        for (vid, &is_rev) in is_reversed_wgt.iter() {
            if let Some(&is_rev_) = is_reversed.get(vid) {
                if is_rev != is_rev_ {
                    panic!("Alleles are reversed within fwgt in dwgt.");
                }
            }
        }
        // merge
        is_reversed.extend(is_reversed_wgt.into_iter());
        //vec::or_bool_vec_mut(&mut is_reversed, is_reversed_);
    }
    is_reversed
}

fn is_reversed_to_vec_bool(
    is_reversed_hashmap: &HashMap<String, bool>,
    snvs: &[SnvId],
    use_snv_pos: bool,
) -> Vec<bool> {
    let is_reversed = snvs
        .iter()
        .map(|snv| {
            let vida = snv.vida(use_snv_pos);
            if let Some(&is_rev) = is_reversed_hashmap.get(vida) {
                is_rev
            } else {
                false
            }
        })
        .collect::<Vec<bool>>();
    is_reversed
}

pub fn merge_freq(snvs_in: &mut Snvs, freq: Snvs, use_snv_pos: bool) {
    // map:
    // freq_in vida -> index in snvs_in
    // freq_in vidar -> index in snvs_in
    let mut vid_to_m = HashMap::with_capacity(snvs_in.snvs_both_n());

    for (si, s) in snvs_in.snv_ids().iter().enumerate() {
        vid_to_m.insert(s.vida(use_snv_pos), si);
        vid_to_m.insert(s.vidar(use_snv_pos), si);
    }

    // println!("{:?}", vid_to_m);

    // flip
    let mut freq_in_snvs: Vec<f64> = vec![f64::NAN; snvs_in.snvs_both_n()];

    // assume freq.freq_alt() is Some
    let freqs_ref = freq.freqs_a2();

    for (si, (s, f)) in freq
        .snv_ids()
        .iter()
        .zip(freq.freqs_a1().unwrap().iter())
        .enumerate()
    {
        // println!("{} {}", s, s.vida(use_snv_pos));

        let m = match vid_to_m.get(s.vida(use_snv_pos)) {
            Some(&m) => m,
            // vid in freq not in Snvs
            None => continue,
        };

        // if !vid_to_m.contains_key(s.vida(use_snv_pos)) {
        //     // vid in freq not in Snvs
        //     continue;
        // }
        // let m = vid_to_m[s.vida(use_snv_pos)];

        let s_in = snvs_in.snv_ids()[m].clone();
        // println!("{} {}", s, &s_in);

        if s.is_rev(&s_in, use_snv_pos).unwrap() {
            // if rev
            // ref: wgts_index_reversed_hashmap_snv()
            match freqs_ref {
                Some(freq_ref) => freq_in_snvs[m] = freq_ref[si],
                // assume bi-allelic
                None => freq_in_snvs[m] = 1.0 - *f,
            }
        } else {
            // if match
            // This considers use_snv_pos since use_snv_pos is in vid_to_m
            freq_in_snvs[m] = *f;
        }

        // should never violate this ==
        // if s == &s_in {
        //     if s.is_rev(&s_in, use_snv_pos) {
        //         // if rev
        //         // ref: wgts_index_reversed_hashmap_snv()
        //         freq_in_snvs[m] = 1.0 - *f;
        //     } else {
        //         // if match
        //         // This considers use_snv_pos since use_snv_pos is in vid_to_m
        //         freq_in_snvs[m] = *f;
        //     }
        // }
        // else ignore
    }

    if freq_in_snvs.iter().any(|x| x.is_nan()) {
        panic!("Some variants in fin_genot has not registered in fin_freq")
    }

    snvs_in.set_maf(freq_in_snvs)
}

pub fn generate_genot_for_score_multiwgts(
    fins_genot: &GenotFiles,
    wgts_multi: &mut [Wgts],
    //wgts: &mut [Wgt],
    freq_buf: Option<&[u8]>,
    n: usize,
    use_samples: Option<&[bool]>,
    // fill_missing_in_test: bool, // for backward compatibility
    fill_missing_group: Option<FillMissingGroup>,
    allow_nonexist_snv: bool,
    use_snv_pos: bool,
    mem: Option<usize>,
) -> Genot {
    let (snvs_in, file_snv_allele_idx_in) = genot_io::load_snvs(fins_genot);
    let mut snvs_in = Snvs::new_from_snv_ids(snvs_in);

    match freq_buf {
        Some(x) => {
            let freqs = genot_io::load_freq(x);
            merge_freq(&mut snvs_in, freqs, use_snv_pos)
        }
        None => {}
    };

    let mut uses_snvs = Vec::new();
    for wgts in wgts_multi.iter_mut() {
        // length is m_in
        let use_snvs = create_wgt_to_genotype_index(snvs_in.snv_ids(), wgts.wgts(), use_snv_pos);
        uses_snvs.push(use_snvs);
    }
    let use_snvs = integrate_use_snvs(&uses_snvs);

    let m = vec::count_true(&use_snvs);
    // if m == 0 {
    //     panic!("Using snvs are zero. Please check weight and genotype file.")
    // }

    let snvs_use = snvs_in.clone().extract_snvs(&use_snvs);
    // let snvs_use: Vec<SnvId> = vec::extract_if_into_iter(snvs_in, &use_snvs);

    // log::info!("Temporary: Skip index reversed hashmap for group");
    let mut iss_reversed_hashmap = Vec::new();
    for wgts in wgts_multi.iter_mut() {
        // length is m
        // wgt vida -> is_reversed
        let is_reversed_hashmap = set_wgts_index_reversed_hashmap(
            wgts.wgts_mut(),
            &snvs_use,
            allow_nonexist_snv,
            use_snv_pos,
        );
        iss_reversed_hashmap.push(is_reversed_hashmap);
    }

    let is_reversed_hashmap_integ = integrate_is_reversed_hashmap_check(iss_reversed_hashmap);
    if is_reversed_hashmap_integ.len() == 0 {
        // this happens for group only
        // panic!("No variants to be loaded after matching alelles.")
    }

    let is_reversed =
        is_reversed_to_vec_bool(&is_reversed_hashmap_integ, &snvs_use.snv_ids(), use_snv_pos);
    log::debug!("Reversed snvs: {}", vec::count_true(&is_reversed));

    let (snvs_group, group_to_m_in) = set_group_index(wgts_multi, snvs_in, m);

    log::info!("Temporary: Skip index reversed hashmap for group");

    if m == 0 && snvs_group.len() == 0 {
        panic!("No variants to be loaded from genotype file.")
    }

    // TOFIX: should not fill_missing_mode in score
    let mut g = load::generate_genot(
        fins_genot,
        m,
        n,
        Some(&use_snvs),
        &file_snv_allele_idx_in,
        group_to_m_in,
        // None,
        use_samples,
        None,
        // fill_missing_in_test,
        fill_missing_group,
        mem,
    );

    // TMP
    //for snv in g.iter_snv() {
    //    println!("snv bfr reverse {:?}", &snv.vals()[..10]);
    //}

    reverse_genotypes_rayon(&mut g.as_genot_mut(), &is_reversed);

    g
}

fn set_group_index(
    wgts_multi: &mut [Wgts],
    snvs_in: Snvs,
    m: usize,
) -> (Vec<SnvId>, Option<Vec<Vec<usize>>>) {
    // extract group from wgt
    let snvs_group: Vec<SnvId> = wgts_multi
        .iter()
        .map(|wgts| {
            wgts.wgts()
                .iter()
                .filter(|wgt| wgt.kind().is_snv_single())
                .filter(|wgt| wgt.kind().snv_id().is_group())
                .map(|wgt| wgt.kind().snv_id().clone())
        })
        .flatten()
        .collect();
    let snvs_group = dedup_snvs_group(&snvs_group);

    let group_to_m_in = if snvs_group.len() == 0 {
        None
    } else {
        Some(snv::make_group_to_m_in_vec(&snvs_group, snvs_in.snv_ids()))
    };

    // let mut group_mi = m;
    let group_id_to_mi: HashMap<&str, usize> = snvs_group
        .iter()
        .enumerate()
        .map(|(group_id_i, group_id)| (group_id.id(), m + group_id_i))
        .collect();

    for wgts in wgts_multi.iter_mut() {
        for wgt in wgts
            .wgts_mut()
            .iter_mut()
            .filter(|wgt| wgt.kind().is_snv_single())
            .filter(|wgt| wgt.kind().snv_id().is_group())
        {
            // TOFIX: fail when wgt is not single snv? -> see load_score::94
            // -> fixed
            let mi = group_id_to_mi[wgt.kind().snv_id().id()];
            wgt.set_snv_index_check(Some(mi));
        }
    }
    (snvs_group, group_to_m_in)
}

/// use in score calculation
/// Judge the corresponding SNVs by chrom, pos, A1, A2.
/// Reversed wgts should have reversed genotypes
/// Flipped genotype do nth.
///
/// genotype_index should be in corresponding order to fin
///     (wgt_to_m, use_snvs, is_flipped): len: .. , m_in, m
///
/// The reason why is_flipped and set_index is here is they require snvs
///
fn create_wgt_to_genotype_index<W: WgtTrait>(
    snvs_in: &[SnvId],
    wgts: &[W],
    use_snv_pos: bool,
) -> Vec<bool> {
    let m_in = snvs_in.len();

    let mut use_snvs = vec![false; m_in];

    // map
    // snv_in vida -> m_in_i
    // snv_in vidar -> m_in_i
    let mut sid_in_to_m_in = HashMap::with_capacity(m_in);

    for (si, s) in snvs_in.iter().enumerate() {
        sid_in_to_m_in.insert(s.vida(use_snv_pos), si);
        sid_in_to_m_in.insert(s.vidar(use_snv_pos), si);
        //sid_in_to_m_in.insert(s.sid(), si);
        //sid_in_to_m_in.insert(s.sida(), si);
    }

    // create use_snvs
    // Even though snv's alelles do not match, use_snvs is true, but no problem.
    for wgt in wgts.iter() {
        match wgt.kind() {
            WgtKind::Snv(snv_wgt) => {
                let snv_id = snv_wgt.snv_id();
                if snv_id.is_group() {
                    // snv only in group is not loaded as snv
                    continue;
                    // snv_id.group_ids().unwrap().iter().for_each(|s| {
                    //     snv_set_use_snvs(s, use_snv_pos, &sid_in_to_m_in, &mut use_snvs);
                    // });
                } else {
                    let vida = snv_id.vida(use_snv_pos);
                    if let Some(v) = sid_in_to_m_in.get(&vida) {
                        use_snvs[*v] = true;
                    }
                }
                // skip if snv in wgt not in snvs_in
            }
            WgtKind::SnvInteraction(snv_inter_wgt) => {
                let (snv_id_1, snv_id_2) = snv_inter_wgt.snv_ids();
                //let vid = wgt.kind().snv_index().vid(use_snv_pos);
                let vida_1 = snv_id_1.vida(use_snv_pos);
                let vida_2 = snv_id_2.vida(use_snv_pos);
                if let Some(v) = sid_in_to_m_in.get(&vida_1) {
                    use_snvs[*v] = true;
                }
                if let Some(v) = sid_in_to_m_in.get(&vida_2) {
                    use_snvs[*v] = true;
                }
                // skip if snv in wgt not in snvs_in
            }
            // if not snv, skip
            WgtKind::Cov(..) => {}
        }

        //if wgt.kind().is_snv() {
        //    let vid = wgt.kind().snv_index().vid(use_snv_pos);
        //    if let Some(v) = sid_in_to_m_in.get(&vid) {
        //        use_snvs[*v] = true;
        //    }
        //}
    }

    use_snvs
}

// TO DELETE
/// use in score calculation
/// Judge the corresponding SNVs by chrom, pos, A1, A2.
/// Reversed wgts should have reversed genotypes
/// Flipped genotype do nth.
///
/// genotype_index should be in corresponding order to fin
///     (wgt_to_m, use_snvs, is_flipped): len: .. , m_in, m
///
/// The reason why is_flipped and set_index is here is they require snvs
///
//fn create_wgt_to_genotype_index<W: WgtTrait>(
//    fin: &Path,
//    gfmt: GenotFormat,
//    wgts: &[W],
//) -> Vec<bool> {
//    let m_in: usize = io_genot::compute_num_snv(fin, gfmt).unwrap();
//    let snvs_in = io_genot::load_snvs(fin, gfmt);
//
//    let mut use_snvs = vec![false; m_in];
//
//    // map: sid_in -> genotype index (m)
//    let mut sid_in_to_m_in = HashMap::with_capacity(m_in);
//
//    for (si, s) in snvs_in.iter().enumerate() {
//        sid_in_to_m_in.insert(s.to_sid(), si);
//        //sid_in_to_m_in.insert(s.sida(), si);
//    }
//
//    // create use_snvs
//    // Even though snv's alelles do not match, use_snvs is true, but no problem.
//    for wgt in wgts.iter() {
//        // if not snv, skip
//        if wgt.kind().is_snv() {
//            let sid = wgt.kind().snv_index().to_sid();
//            if let Some(v) = sid_in_to_m_in.get(&sid) {
//                use_snvs[*v] = true;
//            }
//        }
//    }
//
//    /*
//    let (m_in_to_m, m) = vec::create_hashmap_index_from_use(&use_snvs);
//
//    // needs to reverse genotype
//    let mut is_reversed = vec![false; m];
//
//    //  Not here, want to make load_dataset(&wgts) not &mut
//    // boosting.rs:wgt_sida_to_dataset_m()
//    // create wgt_to_m
//    // is_reversed will be overwritten many times but should be all right if all of the same in wgt
//    for wgt in wgts.iter_mut() {
//        if wgt.kind().is_snv() {
//            let sid = wgt.kind().snv_index().to_sid();
//            //log::debug!("sida {}", sida);
//            //log::debug!("m_in {:?}", sida_in_to_m_in.get(sida));
//            if let Some(v) = sid_in_to_m_in.get(&sid) {
//                let mi = m_in_to_m[v];
//
//                // clone is necessary
//                let snv_in = snvs_in[*v].clone();
//                // clone is necessary
//                let snv_wgt = wgt.kind().snv_index().clone();
//
//                // check alleles match
//                // otherwise alleles do not match
//                if &snv_wgt == &snv_in {
//                    wgt.set_snv_index(Some(mi));
//
//                    // rev or not
//                    // here since putting inside of if raises error
//                    let is_rev = snv_wgt.is_rev(&snv_in);
//
//                    if is_rev {
//                        log::debug!(
//                            "Alleles are reversed in wgt and fplink {:?}, {:?}",
//                            &snv_wgt,
//                            &snv_in
//                        );
//                    }
//
//                    is_reversed[mi] = is_rev;
//                    //is_reversed[mi] = false;
//                } else {
//                    log::info!(
//                        "Ignore SNV: alleles do not match in wgt and fplink {:?}, {:?}",
//                        snv_wgt,
//                        &snv_in
//                    );
//                }
//            } else {
//                wgt.set_snv_index(None);
//            }
//        }
//
//        // do nth for cov
//        // not found
//        //wgt_to_m.insert(wgti, None);
//    }
//    */
//
//    use_snvs
//}

// use set_wgts_index_reversed_hashmap()
//fn set_wgts_index_reversed<W: WgtTrait>(
//    wgts: &mut [W],
//    snvs: &[SnvId],
//    allow_nonexist_snv: bool,
//    use_snv_pos: bool,
//) -> (Vec<bool>, Vec<bool>) {
//    let m = snvs.len();
//
//    // map: sid_in -> genotype index (m)
//    let mut vid_to_m = HashMap::with_capacity(m);
//
//    for (si, s) in snvs.iter().enumerate() {
//        if vid_to_m.contains_key(s.vid(use_snv_pos)) {
//            panic!(
//                "Snv in genot is duplicated. If this is intentional, supress --use-snv_pos: {}",
//                s.vid(use_snv_pos)
//            )
//        }
//        vid_to_m.insert(s.vid(use_snv_pos), si);
//    }
//
//    // needs to reverse genotype
//    let mut is_reversed = vec![false; m];
//
//    // for check_reversed
//    let mut used_snvs = vec![false; m];
//
//    // TODO
//    //  Not here, want to make load_dataset(&wgts) not &mut
//    // boosting.rs:wgt_sida_to_dataset_m()
//    // create wgt_to_m
//    // is_reversed will be overwritten many times but should be all right if all of the same in wgt
//    for wgt in wgts.iter_mut() {
//        if wgt.kind().is_snv() {
//            let vid = wgt.kind().snv_index().vid(use_snv_pos);
//
//            //log::debug!("sida {}", sida);
//            //log::debug!("m_in {:?}", sida_in_to_m_in.get(sida));
//            if let Some(&mi) = vid_to_m.get(&vid) {
//                //let mi = m_in_to_m[v];
//
//                // TOFIX: bug; when alleles do not match, used_snvs should be false
//                // fixed in HashMap ver.
//                used_snvs[mi] = true;
//
//                // clone is necessary
//                let snv_in = snvs[mi].clone();
//                // clone is necessary
//                let snv_wgt = wgt.kind().snv_index().clone();
//
//                // check alleles match
//                if &snv_wgt == &snv_in {
//                    wgt.set_snv_index(Some(mi));
//
//                    // rev or not
//                    // here since putting inside of if raises error
//                    let is_rev = snv_wgt.is_rev(&snv_in, use_snv_pos);
//
//                    if is_rev {
//                        log::debug!(
//                            "Alleles are reversed in wgt and fplink: {:?}, {:?}",
//                            &snv_wgt,
//                            &snv_in
//                        );
//                    }
//
//                    is_reversed[mi] = is_rev;
//                    //is_reversed[mi] = false;
//                } else {
//                    //  alleles do not match
//                    log::info!(
//                        "Alleles do not match in wgt and fplink: {:?}, {:?}",
//                        snv_wgt,
//                        &snv_in
//                    );
//                    if allow_nonexist_snv {
//                        wgt.set_snv_index(None);
//                    } else {
//                        panic!("Alleles do not match in wgt and fplink. Use --allow-nonexist-snv.");
//                    }
//                }
//            } else {
//                log::info!(
//                    "SNV in wgt is not in fplink: {:?}",
//                    wgt.kind().snv_index().clone()
//                );
//                if allow_nonexist_snv {
//                    wgt.set_snv_index(None);
//                } else {
//                    panic!("SNV in wgt is not in fplink. Use --allow-nonexist-snv. --use-snv-pos might help as well.");
//                }
//            }
//        }
//
//        // do nth for cov
//        // not found
//        //wgt_to_m.insert(wgti, None);
//    }
//
//    (is_reversed, used_snvs)
//}

//fn set_wgts_index<W: WgtTrait>(
//    wgts: &mut [W],
//    snvs: &[SnvId],
//    allow_nonexist_snv: bool,
//    use_snv_pos: bool,
//    vid_to_m: &HashMap<String, usize>,
//) {
//    // Set index of snv in wgt
//    for wgt in wgts.iter_mut() {
//        if wgt.kind().is_snv() {
//            let vid = wgt.kind().snv_index().vid(use_snv_pos);
//
//            //log::debug!("sida {}", sida);
//            //log::debug!("m_in {:?}", sida_in_to_m_in.get(sida));
//            if let Some(&mi) = vid_to_m.get(vid) {
//                let snv_in = &snvs[mi];
//                // clone is necessary
//                //let snv_in = snvs[mi].clone();
//                let snv_wgt = wgt.kind().snv_index();
//                // clone is necessary
//                //let snv_wgt = wgt.kind().snv_index().clone();
//
//                // check alleles match
//                if &snv_wgt == &snv_in {
//                    wgt.set_snv_index(Some(mi));
//                } else {
//                    //  alleles do not match
//                    log::info!(
//                        "Alleles do not match in wgt and fplink: {:?}, {:?}",
//                        &snv_wgt,
//                        &snv_in
//                    );
//                    if allow_nonexist_snv {
//                        wgt.set_snv_index(None);
//                    } else {
//                        panic!("Alleles do not match in wgt and fplink. Use --allow-nonexist-snv.");
//                    }
//                }
//            } else {
//                log::info!(
//                    "SNV in wgt is not in fplink: {:?}",
//                    wgt.kind().snv_index().clone()
//                );
//                if allow_nonexist_snv {
//                    wgt.set_snv_index(None);
//                } else {
//                    panic!("SNV in wgt is not in fplink. Use --allow-nonexist-snv. --use-snv-pos might help as well.");
//                }
//            }
//        }
//    }
//}
//
//fn add_is_reversed(
//    is_reversed: &mut HashMap<String, bool>,
//    wgt_snv_index: &SnvId,
//    snvs: &[SnvId],
//    vid_to_m: &HashMap<String, usize>,
//    use_snv_pos: bool,
//) {
//    let vid = wgt_snv_index.vid(use_snv_pos);
//    //log::debug!("sida {}", sida);
//    //log::debug!("m_in {:?}", sida_in_to_m_in.get(sida));
//
//    if let Some(&mi) = vid_to_m.get(vid) {
//        // bug; when alleles do not match, used_snvs should be false
//        // fixed in HashMap ver.
//        //used_snvs[mi] = true;
//
//        let snv_in = &snvs[mi];
//        //let wgt_snv_index = wgt.kind().snv_index();
//
//        // check alleles match
//        if &wgt_snv_index == &snv_in {
//            //wgt.set_snv_index(Some(mi));
//
//            // rev or not
//            // here since putting inside of if raises error
//            let is_rev = wgt_snv_index.is_rev(&snv_in, use_snv_pos);
//
//            if is_rev {
//                log::debug!(
//                    "Alleles are reversed in wgt and fplink: {:?}, {:?}",
//                    &wgt_snv_index,
//                    &snv_in
//                );
//            }
//
//            is_reversed.insert(vid.to_string(), is_rev);
//            //is_reversed[mi] = is_rev;
//        }
//        // else, already done in set_wgts_index()
//    }
//}
//
//fn create_is_reversed_hashmap<W: WgtTrait>(
//    wgts: &[W],
//    snvs: &[SnvId],
//    use_snv_pos: bool,
//    vid_to_m: &HashMap<String, usize>,
//) -> HashMap<String, bool> {
//    // Create is_reversed
//    let mut is_reversed = HashMap::new();
//    for wgt in wgts.iter() {
//        match wgt.kind() {
//            WgtKind::Snv(wgt_snv_index, _, _) => {
//                //let wgt_snv_index = wgt.kind().snv_index();
//                add_is_reversed(&mut is_reversed, wgt_snv_index, snvs, vid_to_m, use_snv_pos)
//
//                //let vid = wgt.kind().snv_index().vid(use_snv_pos);
//
//                ////log::debug!("sida {}", sida);
//                ////log::debug!("m_in {:?}", sida_in_to_m_in.get(sida));
//                //if let Some(&mi) = vid_to_m.get(vid) {
//                //    // TOFIX: bug; when alleles do not match, used_snvs should be false
//                //    // fixed in HashMap ver.
//                //    //used_snvs[mi] = true;
//
//                //    let snv_in = &snvs[mi];
//                //    let snv_wgt = wgt.kind().snv_index();
//
//                //    // check alleles match
//                //    if &snv_wgt == &snv_in {
//                //        //wgt.set_snv_index(Some(mi));
//
//                //        // rev or not
//                //        // here since putting inside of if raises error
//                //        let is_rev = snv_wgt.is_rev(&snv_in, use_snv_pos);
//
//                //        if is_rev {
//                //            log::debug!(
//                //                "Alleles are reversed in wgt and fplink: {:?}, {:?}",
//                //                &snv_wgt,
//                //                &snv_in
//                //            );
//                //        }
//
//                //        is_reversed.insert(vid.to_string(), is_rev);
//                //        //is_reversed[mi] = is_rev;
//                //    }
//                //    // else, already done in set_wgts_index()
//                //}
//                // else, already done in set_wgts_index()
//            }
//            WgtKind::SnvInteraction(wgt_snv_index_1, _, wgt_snv_index_2, _) => {
//                add_is_reversed(
//                    &mut is_reversed,
//                    wgt_snv_index_1,
//                    snvs,
//                    vid_to_m,
//                    use_snv_pos,
//                );
//                add_is_reversed(
//                    &mut is_reversed,
//                    wgt_snv_index_2,
//                    snvs,
//                    vid_to_m,
//                    use_snv_pos,
//                );
//            }
//            _ => {}
//        }
//    }
//    is_reversed
//}

fn wgts_index_reversed_hashmap_snv(
    wgt_snv_id: &SnvId,
    snvs: &[SnvId],
    vida_to_m: &HashMap<String, usize>,
    allow_nonexist_snv: bool,
    use_snv_pos: bool,
) -> Option<(usize, bool)> {
    let vida = wgt_snv_id.vida(use_snv_pos);

    if let Some(&mi) = vida_to_m.get(vida) {
        let snv_id_in = &snvs[mi];

        // rev or not
        let is_rev = wgt_snv_id.is_rev(&snv_id_in, use_snv_pos).unwrap();

        if is_rev {
            log::debug!(
                "Alleles are reversed in wgt and fplink: {:?}, {:?}",
                &wgt_snv_id,
                &snv_id_in
            );
        }

        Some((mi, is_rev))
    } else {
        log::info!("SNV in wgt is not in fplink: {:?}", wgt_snv_id);
        if allow_nonexist_snv {
            None
        } else {
            panic!("SNV in wgt is not in fplink. Use --allow-nonexist-snv. --use-snv-pos might help as well.");
        }
    }
}

fn set_wgts_index_reversed_hashmap_snv<W: WgtTrait>(
    wgt: &mut W,
    is_reversed: &mut HashMap<String, bool>,
    wgt_snv_id: &SnvId,
    snvs: &Snvs,
    //snvs: &[SnvId],
    vida_to_m: &HashMap<String, usize>,
    allow_nonexist_snv: bool,
    use_snv_pos: bool,
) {
    let index_reverse = wgts_index_reversed_hashmap_snv(
        wgt_snv_id,
        snvs.snv_ids(),
        vida_to_m,
        allow_nonexist_snv,
        use_snv_pos,
    );

    // TODO: check consistency before inserting is_rev
    if let Some((mi, is_rev)) = index_reverse {
        wgt.set_snv_index_check(Some(mi));

        match snvs.freqs_a1() {
            Some(x) => {
                let f = x[mi];
                //let f = snvs.mafs().unwrap()[mi];
                let f = if is_rev {
                    // TOFIX: add --assume-bi-allelic to skip this unimplemented! for common variant pgs
                    // For now, comment below for common variant pgs
                    unimplemented!("multi-allelic; load REF_FREQ from freq");
                    1.0 - f
                } else {
                    f
                };
                wgt.set_freq(Some(f));
            }
            None => {}
        }

        update_is_reversed_hashmap(is_reversed, wgt_snv_id, is_rev, use_snv_pos);
    }
    // TODO: is this necessary?; yes if overwrite
    // else{
    //wgt.set_snv_index(None);
    //}
}

fn update_is_reversed_hashmap(
    is_reversed: &mut HashMap<String, bool>,
    wgt_snv_id: &SnvId,
    is_rev: bool,
    use_snv_pos: bool,
) {
    let vida = wgt_snv_id.vida(use_snv_pos).to_string();
    match is_reversed.get(&vida) {
        None => {
            is_reversed.insert(vida, is_rev);
        }
        Some(v) => {
            if *v != is_rev {
                panic!("Alleles are reversed within fwgt.");
            }
        }
    }
}

fn set_wgts_index_reversed_hashmap_snv_interaction<W: WgtTrait>(
    wgt: &mut W,
    is_reversed: &mut HashMap<String, bool>,
    wgt_snv_id_1: &SnvId,
    wgt_snv_id_2: &SnvId,
    snvs: &Snvs,
    //snvs: &[SnvId],
    vida_to_m: &HashMap<String, usize>,
    allow_nonexist_snv: bool,
    use_snv_pos: bool,
) {
    let index_reverse_1 = wgts_index_reversed_hashmap_snv(
        wgt_snv_id_1,
        snvs.snv_ids(),
        vida_to_m,
        allow_nonexist_snv,
        use_snv_pos,
    );

    let index_reverse_2 = wgts_index_reversed_hashmap_snv(
        wgt_snv_id_2,
        snvs.snv_ids(),
        vida_to_m,
        allow_nonexist_snv,
        use_snv_pos,
    );

    match (index_reverse_1, index_reverse_2) {
        (Some((mi_1, is_rev_1)), Some((mi_2, is_rev_2))) => {
            wgt.set_snv_index_interaction_check(Some(mi_1), Some(mi_2));
            update_is_reversed_hashmap(is_reversed, wgt_snv_id_1, is_rev_1, use_snv_pos);
            update_is_reversed_hashmap(is_reversed, wgt_snv_id_2, is_rev_2, use_snv_pos);
            //is_reversed.insert(wgt_snv_id_1.vid(use_snv_pos).to_string(), is_rev_1);
            //is_reversed.insert(wgt_snv_id_2.vid(use_snv_pos).to_string(), is_rev_2);
        }
        // do nth
        _ => {}
    }

    //if let Some((mi, is_rev)) = index_reverse_1 {
    //    wgt.set_snv_index(Some(mi));
    //    is_reversed.insert(wgt_snv_id_1.vid(use_snv_pos).to_string(), is_rev);
    //}
    //if let Some((mi, is_rev)) = index_reverse_2 {
    //    wgt.set_snv_index(Some(mi));
    //    is_reversed.insert(wgt_snv_id_2.vid(use_snv_pos).to_string(), is_rev);
    //}
}

//fn set_wgts_index_reversed_hashmap_snv_old<W: WgtTrait>(
//    wgt: &mut W,
//    is_reversed: &mut HashMap<String, bool>,
//    wgt_snv_id: &SnvId,
//    snvs: &[SnvId],
//    vid_to_m: &HashMap<String, usize>,
//    allow_nonexist_snv: bool,
//    use_snv_pos: bool,
//) {
//    let vid = wgt_snv_id.vid(use_snv_pos);
//
//    if let Some(&mi) = vid_to_m.get(vid) {
//        let snv_id_in = &snvs[mi];
//
//        // check alleles match
//        if &wgt_snv_id == &snv_id_in {
//            wgt.set_snv_index(Some(mi));
//
//            // rev or not
//            let is_rev = wgt_snv_id.is_rev(&snv_id_in, use_snv_pos);
//
//            if is_rev {
//                log::debug!(
//                    "Alleles are reversed in wgt and fplink: {:?}, {:?}",
//                    &wgt_snv_id,
//                    &snv_id_in
//                );
//            }
//
//            is_reversed.insert(vid.to_string(), is_rev);
//            //is_reversed[mi] = is_rev;
//        } else {
//            //  alleles do not match
//            log::info!(
//                "Alleles do not match in wgt and fplink: {:?}, {:?}",
//                &wgt_snv_id,
//                &snv_id_in
//            );
//            if allow_nonexist_snv {
//                wgt.set_snv_index(None);
//            } else {
//                panic!("Alleles do not match in wgt and fplink. Use --allow-nonexist-snv.");
//            }
//        }
//    } else {
//        log::info!(
//            "SNV in wgt is not in fplink: {:?}",
//            wgt.kind().snv_id().clone()
//        );
//        if allow_nonexist_snv {
//            wgt.set_snv_index(None);
//        } else {
//            panic!("SNV in wgt is not in fplink. Use --allow-nonexist-snv. --use-snv-pos might help as well.");
//        }
//    }
//}

fn set_wgts_index_reversed_hashmap_both<W: WgtTrait>(
    wgts: &mut [W],
    snvs: &Snvs,
    //snvs: &[SnvId],
    allow_nonexist_snv: bool,
    use_snv_pos: bool,
    vida_to_m: &HashMap<String, usize>,
) -> HashMap<String, bool> {
    // wgt vida -> is_reversed
    let mut is_reversed = HashMap::new();
    for wgt in wgts.iter_mut() {
        let wgt_kind = wgt.kind().clone();
        //match wgt.kind() {
        match wgt_kind {
            //WgtKind::Snv(wgt_snv_id, _, _) => {
            WgtKind::Snv(snv_wgt) => {
                let snv_id = snv_wgt.snv_id();
                //let wgt_snv_index = wgt.kind().snv_index();
                if snv_id.is_group() {
                    continue;
                    // snv_id.group_ids().unwrap().iter().for_each(|s| {
                    // log::info!("Skip index reversed hashmap for group");
                    // TODO: implement after enum SnvId
                    // set_wgts_index_reversed_hashmap_snv(
                    //     wgt,
                    //     &mut is_reversed,
                    //     s,
                    //     snvs,
                    //     vida_to_m,
                    //     allow_nonexist_snv,
                    //     use_snv_pos,
                    // );
                    // });
                } else {
                    set_wgts_index_reversed_hashmap_snv(
                        wgt,
                        &mut is_reversed,
                        &snv_id,
                        snvs,
                        vida_to_m,
                        allow_nonexist_snv,
                        use_snv_pos,
                    );
                }
            }
            //WgtKind::SnvInteraction(snv_id_1, _, snv_id_2, _) => {
            WgtKind::SnvInteraction(snv_inter_wgt) => {
                let (snv_id_1, snv_id_2) = snv_inter_wgt.snv_ids();
                set_wgts_index_reversed_hashmap_snv_interaction(
                    wgt,
                    &mut is_reversed,
                    &snv_id_1,
                    &snv_id_2,
                    snvs,
                    vida_to_m,
                    allow_nonexist_snv,
                    use_snv_pos,
                );
            }
            // cov
            _ => {}
        }
    }
    is_reversed
}

// TODO: split into two parts
fn set_wgts_index_reversed_hashmap<W: WgtTrait>(
    wgts: &mut [W],
    snvs: &Snvs,
    //snvs: &[SnvId],
    allow_nonexist_snv: bool,
    use_snv_pos: bool,
) -> HashMap<String, bool> {
    let m = snvs.snvs_both_n();
    //let m = snvs.len();

    // map:
    // snv_in vida -> genotype index (m)
    // snv_in vidar -> genotype index (m)
    let mut vid_to_m = HashMap::with_capacity(m);

    for (si, s) in snvs.snv_ids().iter().enumerate() {
        //if vid_to_m.contains_key(s.vida(use_snv_pos)) {
        //    panic!(
        //        "Variant ID in genot is duplicated. If this is intentional, supress --use-snv-pos: {}",
        //        s.vid(use_snv_pos)
        //    )
        //}

        if s.is_group() {
            continue;
        } else {
            if vid_to_m
                .insert(s.vida(use_snv_pos).to_string(), si)
                .is_some()
            {
                panic!(
                "Variant ID in genot is duplicated. If this is intentional, supress --use-snv-pos: {}",
                s.vida(use_snv_pos)
            )
            }
            if vid_to_m
                .insert(s.vidar(use_snv_pos).to_string(), si)
                .is_some()
            {
                panic!(
                "Variant ID in genot is duplicated. If this is intentional, supress --use-snv-pos: {}",
                s.vidar(use_snv_pos)
            )
            }
        }
    }

    // Simultaniously set wgt index and create is_reversed
    // since conditional branch is complicated
    let is_reversed = set_wgts_index_reversed_hashmap_both(
        wgts,
        snvs,
        allow_nonexist_snv,
        use_snv_pos,
        &vid_to_m,
    );

    // split ver.
    //// Set index of snv in wgt
    //set_wgts_index(wgts, snvs, allow_nonexist_snv, use_snv_pos, &vid_to_m);
    //// Set index of snv in wgt
    //let is_reversed = create_is_reversed_hashmap(wgts, snvs, use_snv_pos, &vid_to_m);

    is_reversed
}

// move to somewhere else (in GenotMut)
/// rayon ver.
fn reverse_genotypes_rayon(g: &mut GenotMut, is_reversed: &[bool]) {
    g.iter_snv_mut()
        .zip(is_reversed.iter())
        .par_bridge()
        .filter(|(_, b)| **b)
        .for_each(|(mut gsnv, _)| reverse_genotypes_snv(&mut gsnv));
}

// {0:2, 1:1, 2:0, 3:3}
const REVERSE_GENOTYPE_AR: [u8; 4] = [2, 1, 0, 3];

#[inline]
fn reverse_genotypes_snv(gsnv: &mut GenotSnvMut) {
    for ni in 0..gsnv.n() {
        let v = gsnv.get_val_unchecked(ni) as usize;
        gsnv.set_unchecked(REVERSE_GENOTYPE_AR[v], ni);
    }
}

#[cfg(test)]
mod tests {
    //use crate::snv::snv_index;

    use super::*;

    #[test]
    fn test_integrate_is_reversed_hashmap_check() {
        // 0,2 are flipped
        let tuples_1 = [("0".to_string(), true), ("1".to_string(), false)];
        let tuples_2 = [("1".to_string(), false), ("2".to_string(), true)];
        let tuples_3 = [("0".to_string(), true), ("3".to_string(), false)];

        let tuples = vec![tuples_1, tuples_2, tuples_3];
        let is_reversed = tuples
            .into_iter()
            .map(|x| x.into_iter().collect::<HashMap<String, bool>>())
            .collect::<Vec<_>>();

        let ans = integrate_is_reversed_hashmap_check(is_reversed);

        let exp = [
            ("0".to_string(), true),
            ("1".to_string(), false),
            ("2".to_string(), true),
            ("3".to_string(), false),
        ]
        .into_iter()
        .collect::<HashMap<String, bool>>();

        assert_eq!(ans, exp);
    }

    #[test]
    #[should_panic]
    fn test_integrate_is_reversed_hashmap_check_panic() {
        // 0,2 are flipped
        // but tuples_1 says snv0 is not flipped
        let tuples_1 = [("0".to_string(), false), ("1".to_string(), false)];
        let tuples_2 = [("1".to_string(), false), ("2".to_string(), true)];
        let tuples_3 = [("0".to_string(), true), ("3".to_string(), false)];

        let tuples = vec![tuples_1, tuples_2, tuples_3];
        let is_reversed = tuples
            .into_iter()
            .map(|x| x.into_iter().collect::<HashMap<String, bool>>())
            .collect::<Vec<_>>();

        let _ = integrate_is_reversed_hashmap_check(is_reversed);
    }

    #[test]
    fn test_merge_freq() {
        let snv_id_1 = SnvId::new(
            "rs1".to_owned(),
            "1".to_string(),
            "1",
            "A".to_owned(),
            "G".to_owned(),
        );
        let snv_id_2 = SnvId::new(
            "rs2".to_owned(),
            "1".to_string(),
            "2",
            "A".to_owned(),
            "G".to_owned(),
        );
        let snv_id_3 = SnvId::new(
            "rs3".to_owned(),
            "1".to_string(),
            "3",
            "A".to_owned(),
            "G".to_owned(),
        );

        let mut snvs_in = Snvs::new_from_snv_ids(vec![snv_id_1, snv_id_2, snv_id_3]);

        // rev
        let snv_id_1 = SnvId::new(
            "rs2".to_owned(),
            "1".to_string(),
            "2",
            "G".to_owned(),
            "A".to_owned(),
        );
        // old: rev & flipped
        let snv_id_2 = SnvId::new(
            "rs3".to_owned(),
            "1".to_string(),
            "3",
            "G".to_owned(),
            "A".to_owned(),
        );
        //let snv_id_2 = SnvId::new("rs3".to_owned(), "1".to_string(), "3", "C".to_owned(), "T".to_owned());
        // same
        let snv_id_3 = SnvId::new(
            "rs1".to_owned(),
            "1".to_string(),
            "1",
            "A".to_owned(),
            "G".to_owned(),
        );
        // not match
        let snv_id_4 = SnvId::new(
            "rs1".to_owned(),
            "1".to_string(),
            "1",
            "A".to_owned(),
            "C".to_owned(),
        );

        let freqs = vec![0.1f64, 0.2, 0.3, 0.4];
        let freqs_in = Snvs::new_snv(
            vec![snv_id_1, snv_id_2, snv_id_3, snv_id_4],
            Some(freqs),
            None,
        );

        merge_freq(&mut snvs_in, freqs_in, true);

        let freqs = snvs_in.freqs_a1().unwrap();
        let freqs_ans = vec![0.3, 0.9, 0.8];

        assert_eq!(freqs, &freqs_ans);
    }

    #[test]
    fn test_merge_freq_multi_allele() {
        let snv_id_1 = SnvId::new(
            "rs1".to_owned(),
            "1".to_string(),
            "1",
            "A".to_owned(),
            "G".to_owned(),
        );
        let snv_id_2 = SnvId::new(
            "rs2".to_owned(),
            "1".to_string(),
            "2",
            "A".to_owned(),
            "G".to_owned(),
        );
        let snv_id_2_2 = SnvId::new(
            "rs2".to_owned(),
            "1".to_string(),
            "2",
            "A".to_owned(),
            "T".to_owned(),
        );
        let snv_id_3 = SnvId::new(
            "rs3".to_owned(),
            "1".to_string(),
            "3",
            "A".to_owned(),
            "G".to_owned(),
        );

        let mut snvs_in = Snvs::new_from_snv_ids(vec![snv_id_1, snv_id_2, snv_id_2_2, snv_id_3]);

        // rev
        let snv_id_1 = SnvId::new(
            "rs2".to_owned(),
            "1".to_string(),
            "2",
            "G".to_owned(),
            "A".to_owned(),
        );
        let snv_id_1_2 = SnvId::new(
            "rs2".to_owned(),
            "1".to_string(),
            "2",
            "T".to_owned(),
            "A".to_owned(),
        );
        // old: rev & flipped
        let snv_id_2 = SnvId::new(
            "rs3".to_owned(),
            "1".to_string(),
            "3",
            "G".to_owned(),
            "A".to_owned(),
        );
        //let snv_id_2 = SnvId::new("rs3".to_owned(), "1".to_string(), "3", "C".to_owned(), "T".to_owned());
        // same
        let snv_id_3 = SnvId::new(
            "rs1".to_owned(),
            "1".to_string(),
            "1",
            "A".to_owned(),
            "G".to_owned(),
        );
        // not match
        let snv_id_4 = SnvId::new(
            "rs1".to_owned(),
            "1".to_string(),
            "1",
            "A".to_owned(),
            "C".to_owned(),
        );

        let freqs = vec![0.1f64, 0.25, 0.2, 0.3, 0.4];
        let freqs_ref = vec![0.8f64, 0.6, 0.7, 0.6, 0.5];
        let freqs_in = Snvs::new_snv(
            vec![snv_id_1, snv_id_1_2, snv_id_2, snv_id_3, snv_id_4],
            Some(freqs),
            Some(freqs_ref),
        );

        merge_freq(&mut snvs_in, freqs_in, true);

        let freqs = snvs_in.freqs_a1().unwrap();
        let freqs_ans = vec![0.3, 0.8, 0.6, 0.7];
        // let freqs_ans = vec![0.3, 0.9, 0.75, 0.8];

        assert_eq!(freqs, &freqs_ans);
    }

    #[test]
    fn test_merge_freq_multi_allele_use_snv_pos() {
        let snv_id_1 = SnvId::new(
            "rs1".to_owned(),
            "1".to_string(),
            "1",
            "A".to_owned(),
            "G".to_owned(),
        );
        let snv_id_2 = SnvId::new(
            "rs2".to_owned(),
            "1".to_string(),
            "2",
            "A".to_owned(),
            "G".to_owned(),
        );
        let snv_id_2_2 = SnvId::new(
            "rs2".to_owned(),
            "1".to_string(),
            "2",
            "A".to_owned(),
            "T".to_owned(),
        );
        let snv_id_3 = SnvId::new(
            "rs3".to_owned(),
            "1".to_string(),
            "3",
            "A".to_owned(),
            "G".to_owned(),
        );

        let mut snvs_in = Snvs::new_from_snv_ids(vec![snv_id_1, snv_id_2, snv_id_2_2, snv_id_3]);

        // rev
        let snv_id_1 = SnvId::new(
            "chr1:2".to_owned(),
            "1".to_string(),
            "2",
            "G".to_owned(),
            "A".to_owned(),
        );
        let snv_id_1_2 = SnvId::new(
            "chr1:2".to_owned(),
            "1".to_string(),
            "2",
            "T".to_owned(),
            "A".to_owned(),
        );
        // old: rev & flipped
        let snv_id_2 = SnvId::new(
            "chr1:3".to_owned(),
            "1".to_string(),
            "3",
            "G".to_owned(),
            "A".to_owned(),
        );
        //let snv_id_2 = SnvId::new("rs3".to_owned(), "1".to_string(), "3", "C".to_owned(), "T".to_owned());
        // same
        let snv_id_3 = SnvId::new(
            "chr1:1".to_owned(),
            "1".to_string(),
            "1",
            "A".to_owned(),
            "G".to_owned(),
        );
        // not match
        let snv_id_4 = SnvId::new(
            "chr1:1".to_owned(),
            "1".to_string(),
            "1",
            "A".to_owned(),
            "C".to_owned(),
        );

        let freqs = vec![0.1f64, 0.25, 0.2, 0.3, 0.4];
        let freqs_ref = vec![0.8f64, 0.6, 0.7, 0.6, 0.5];
        let freqs_in = Snvs::new_snv(
            vec![snv_id_1, snv_id_1_2, snv_id_2, snv_id_3, snv_id_4],
            Some(freqs),
            Some(freqs_ref),
        );

        merge_freq(&mut snvs_in, freqs_in, true);

        let freqs = snvs_in.freqs_a1().unwrap();
        let freqs_ans = vec![0.3, 0.8, 0.6, 0.7];
        // let freqs_ans = vec![0.3, 0.9, 0.75, 0.8];

        assert_eq!(freqs, &freqs_ans);
    }

    //#[test]
    //fn test_new() {
    //    let vec = vec![1, 2, 3, 0];
    //    let mut g = GenotSnv::new(&vec);
    //    reverse_genotypes_snv(&mut g.as_genot_snv_mut_snv());
    //    assert_eq!(g.vals(), vec![1, 0, 3, 2]);
    //}

    //#[test]
    //fn test_check_consistent_is_reversed() {
    //    let useds_snvs = vec![
    //        vec![true, true, false, false, false],
    //        vec![false, false, true, true, true],
    //        vec![false, true, true, false, false],
    //    ];

    //    // 0,2 are flipped
    //    let iss_reversed = vec![
    //        vec![true, false, false, false, false],
    //        vec![false, false, true, false, false],
    //        vec![false, false, true, false, false],
    //    ];

    //    let is_reversed = vec![true, false, true, false, false];

    //    check_consistent_is_reversed(&is_reversed, &iss_reversed, &useds_snvs);
    //}

    //#[should_panic]
    //#[test]
    //fn test_check_consistent_is_reversed_panic() {
    //    let useds_snvs = vec![
    //        vec![true, true, false, false, false],
    //        vec![false, false, true, true, true],
    //        vec![false, true, true, false, false],
    //    ];

    //    // 0,2th snv are flipped
    //    // but 1st wgt says 2nd snv is not flipped
    //    let iss_reversed = vec![
    //        vec![true, false, false, false, false],
    //        vec![false, false, false, false, false],
    //        vec![false, false, true, false, false],
    //    ];

    //    let is_reversed = vec![true, false, true, false, false];

    //    check_consistent_is_reversed(&is_reversed, &iss_reversed, &useds_snvs);
    //}

    /*
    #[test]
    fn test_integrate_use_snvs_reversed() {
        let uses_snvs = vec![
            vec![false, true, false, false, false],
            vec![false, false, true, false, false],
            vec![false, true, true, false, false],
        ];

        let iss_reversed = vec![vec![false], vec![true], vec![false, true]];

        let use_snvs_ans = vec![false, true, true, false, false];
        let is_reversed_ans = vec![false, true];

        let (use_snvs, is_reversed) = integrate_use_snvs_reversed(uses_snvs, iss_reversed);

        assert_eq!(use_snvs, use_snvs_ans);
        assert_eq!(is_reversed, is_reversed_ans);
    }

    #[should_panic]
    #[test]
    fn test_integrate_use_snvs_reversed_panic() {
        let uses_snvs = vec![
            vec![false, true, false, false, false],
            vec![false, false, true, false, false],
            vec![false, true, true, false, false],
        ];

        let iss_reversed = vec![
            vec![false],
            vec![true],
            // this is not consistent so should panic
            vec![false, false],
        ];

        let use_snvs_ans = vec![false, true, true, false, false];
        let is_reversed_ans = vec![false, true];

        let (use_snvs, is_reversed) = integrate_use_snvs_reversed(uses_snvs, iss_reversed);

        assert_eq!(use_snvs, use_snvs_ans);
        assert_eq!(is_reversed, is_reversed_ans);
    }
    */
}
