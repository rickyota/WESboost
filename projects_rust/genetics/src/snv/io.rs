// use rayon::iter::FromParallelIterator;

use super::group_index::AggId;
use super::SnvId;
use crate::{textfile, vec};

use std::collections::{HashMap, HashSet};
use std::path::Path;

/// Use rs to extract snvs.
/// TODO: Format of fin_snv is the same as plink `--extract`.
pub fn make_use_snvs_buf(
    extract_snv_buf: Option<&[u8]>,
    snvs_in: &Vec<SnvId>,
) -> (Vec<bool>, usize) {
    let m_in: usize = snvs_in.len();

    if extract_snv_buf.is_none() {
        let use_snvs = vec![true; m_in];
        return (use_snvs, m_in);
    }

    let snvs_use = load_snvs_use_buf(extract_snv_buf.unwrap());

    make_use_snvs_buf_vec(&snvs_use, snvs_in)
}

pub fn make_use_snvs_agg_buf(
    agg_snv_buf: Option<&[u8]>,
    snvs_in: &Vec<SnvId>,
) -> (Vec<bool>, usize) {
    let m_in: usize = snvs_in.len();

    if agg_snv_buf.is_none() {
        let use_snvs = vec![false; m_in];
        return (use_snvs, m_in);
    }

    let agg_snvs_use = load_agg_snvs_buf(agg_snv_buf);
    let snvs_in_agg = list_snvs_in_agg(&agg_snvs_use);

    make_use_snvs_buf_vec(&snvs_in_agg, snvs_in)
}

pub fn make_use_snvs_group_buf(
    group_snv_buf: Option<&[u8]>,
    snvs_in: &Vec<SnvId>,
) -> (Vec<bool>, usize) {
    let m_in: usize = snvs_in.len();

    if group_snv_buf.is_none() {
        let use_snvs = vec![false; m_in];
        return (use_snvs, 0);
        //let use_snvs = vec![true; m_in];
        //return (use_snvs, m_in);
    }

    //let group_snvs_use = load_group_snvs_buf(group_snv_buf.unwrap());
    let group_snvs_use = load_group_snvs_buf(group_snv_buf);
    let snvs_in_group = list_snvs_in_group(&group_snvs_use);

    make_use_snvs_buf_vec(&snvs_in_group, snvs_in)
    //make_use_snvs_buf_vec(&snvs_in_group, snvs_in)
}

// multi-allelic
// judge by sida if allele is registered,
// judge by rs otherwise
// not allow reversed allele
fn make_use_snvs_buf_vec(snvs_ma_use: &Vec<SnvId>, snvs_ma_in: &Vec<SnvId>) -> (Vec<bool>, usize) {
    if snvs_ma_use.len() == 0 {
        return (vec![false; snvs_ma_in.len()], 0);
    }

    let m_ma_in: usize = snvs_ma_in.len();

    let use_ida = if snvs_ma_use[0].is_alleles_registered() {
        true
    } else {
        false
    };
    // println!("use_ida: {}", use_ida);

    //  check all index; not [0] only
    for snv in snvs_ma_use.iter() {
        if snv.is_alleles_registered() != use_ida {
            panic!("Alleles are not registered consistently.");
        }
    }

    // map: sid_in -> plink index
    // &str should be fine
    let mut rs_in_to_index: HashMap<String, usize> = HashMap::with_capacity(m_ma_in);
    for (si, s) in snvs_ma_in.iter().enumerate() {
        if use_ida {
            rs_in_to_index.insert(s.ida().to_string(), si);
        } else {
            rs_in_to_index.insert(s.id().to_string(), si);
        }
    }
    // println!("rs_in_to_index: {:?}", rs_in_to_index);

    let mut use_snvs = vec![false; m_ma_in];
    let mut m: usize = 0;

    // dedup
    let snvs_rs: HashSet<String> = HashSet::from_iter(snvs_ma_use.iter().map(|x| {
        if use_ida {
            x.ida().to_string()
        } else {
            x.id().to_string()
        }
    }));

    // for snv in snvs_ma_use.iter() {
    for rs in snvs_rs.into_iter() {
        // let rs = if use_ida { snv.ida() } else { snv.id() };
        //let rs = snv.id();

        // println!("rs: {}", rs);

        match rs_in_to_index.get(&rs) {
            Some(v) => {
                use_snvs[*v] = true;
                m += 1;
            }

            // ignore unfound SNVs
            None => {
                log::info!("SNV in fin_snv was not found in plink: {}.", &rs);
            }
        }
    }

    assert_eq!(vec::count_true(&use_snvs), m);
    (use_snvs, m)

    // duplicated
    // for snv in snvs_ma_use.iter() {
    //     let rs = if use_ida { snv.ida() } else { snv.id() };
    //     //let rs = snv.id();

    //     // println!("rs: {}", rs);

    //     match rs_in_to_index.get(rs) {
    //         Some(v) => {
    //             use_snvs[*v] = true;
    //             m += 1;
    //         }

    //         // ignore unfound SNVs
    //         None => {
    //             log::info!("SNV in fin_snv was not found in plink: {}.", rs);
    //         }
    //     }
    // }
    // (use_snvs, m)
}

// judge by rs
//fn make_use_snvs_buf_vec(snvs_use: &Vec<SnvId>, snvs_in: &Vec<SnvId>) -> (Vec<bool>, usize) {
//    let m_in: usize = snvs_in.len();
//
//    // map: sid_in -> plink index
//    // &str should be fine
//    let mut rs_in_to_index = HashMap::with_capacity(m_in);
//    for (si, s) in snvs_in.iter().enumerate() {
//        rs_in_to_index.insert(s.id(), si);
//    }
//
//    let mut use_snvs = vec![false; m_in];
//    let mut m: usize = 0;
//    for snv in snvs_use.iter() {
//        let rs = snv.id();
//
//        match rs_in_to_index.get(rs) {
//            Some(v) => {
//                use_snvs[*v] = true;
//                m += 1;
//            }
//
//            // ignore unfound SNVs
//            None => {
//                log::info!("SNV in fin_snv was not found in plink: {}.", rs);
//            }
//        }
//
//        // if you want to panic
//        //use_snvs[sida_to_index[&sida]] = true;
//    }
//
//    assert_eq!(vec::count_true(&use_snvs), m);
//
//    (use_snvs, m)
//}

pub fn load_snvs_use(fin_snv: &Path) -> Vec<SnvId> {
    let buf = textfile::read_file_to_end(fin_snv, None).unwrap();

    load_snvs_use_buf(&buf[..])
}

/// SnvId with rs only
pub fn load_snvs_use_buf(snv_buf: &[u8]) -> Vec<SnvId> {
    let mut snvs: Vec<SnvId> = vec![];

    let cols = [0usize];
    let vss: Vec<Vec<String>> = textfile::load_table_cols_buf(snv_buf, &cols, false);

    if vss.len() == 0 {
        return snvs;
    }

    for vi in 0..vss[0].len() {
        snvs.push(SnvId::new_id_ma(&vss[0][vi]));
    }

    // TOFIX: check if ids are unique

    snvs
}

pub fn load_agg_snvs_buf(group_snv_buf: Option<&[u8]>) -> Vec<AggId> {
    let mut groups: Vec<AggId> = vec![];

    if group_snv_buf.is_none() {
        return groups;
    }

    // TODO: introduce header
    let cols = [0usize, 1, 2, 3, 4, 5, 6, 7];
    let vss: Vec<Vec<String>> = textfile::load_table_cols_buf(group_snv_buf.unwrap(), &cols, true);
    // let vss: Vec<Vec<String>> = textfile::load_table_cols_buf(group_snv_buf.unwrap(), &cols, false);

    for vi in 0..vss[0].len() {
        let snv_ids_str: Vec<&str> = if vss[7][vi] == "None" {
            vec![]
        } else {
            vss[7][vi].split(",").collect()
        };

        groups.push(AggId::new_ma(
            vss[0][vi].clone(),
            vss[1][vi].clone(),
            &vss[2][vi],
            &vss[3][vi],
            &snv_ids_str,
        ));
    }
    groups
}

// output should be SnvId
pub fn load_group_snvs_buf(group_snv_buf: Option<&[u8]>) -> Vec<SnvId> {
    let mut groups: Vec<SnvId> = vec![];

    if group_snv_buf.is_none() {
        return groups;
    }

    // TODO: use read_csv
    // now same as agg_snvs
    // now 6=num_snvs
    let cols = [0usize, 1, 2, 3, 4, 5, 6, 7];
    let vss: Vec<Vec<String>> = textfile::load_table_cols_buf(group_snv_buf.unwrap(), &cols, true);
    // let vss: Vec<Vec<String>> = textfile::load_table_cols_buf(group_snv_buf.unwrap(), &cols, false);

    for vi in 0..vss[0].len() {
        let group_ids = if vss[7][vi] == "None" {
            vec![]
        } else {
            vss[7][vi].split(",").collect()
            // vss[7][vi]
            //     .clone()
            //     .split(",")
            //     .map(|x| x.to_string())
            //     .collect::<Vec<String>>(),
        };
        groups.push(SnvId::new_group(
            vss[0][vi].clone(),
            vss[1][vi].clone(),
            &vss[2][vi],
            &vss[3][vi],
            group_ids,
        ));
    }

    groups
}

//fn list_snvs_in_group(group_snvs: &Vec<SnvId>) -> Vec<SnvId> {
fn list_snvs_in_agg(agg_snvs: &Vec<AggId>) -> Vec<SnvId> {
    agg_snvs.iter().flat_map(|x| x.snv_ids().clone()).collect()
    // TODO: dedup?
}

fn list_snvs_in_group(group_snvs: &Vec<SnvId>) -> Vec<SnvId> {
    group_snvs
        .iter()
        .flat_map(|x| x.group_ids().unwrap().clone())
        .collect()
    // TODO: dedup?
}

//fn load_group_to_m_in(group_snv_buf: &[u8]) -> usize {
//    let group_snvs = load_group_snvs_buf(group_snv_buf);
//    let snvs_in_group = list_snvs_in_group(&group_snvs);
//
//    snvs_in_group.len()
//}

// see make_agg_to_m_vec()
pub fn make_group_to_m_in_buf(
    group_snv_buf: Option<&[u8]>,
    snvs_in: &[SnvId],
    // snvs_in: &Vec<SnvId>,
) -> (Option<Vec<Vec<usize>>>, usize) {
    // ) -> (Option<HashMap<usize, Vec<usize>>>, usize) {
    //let m_in: usize = snvs_in.len();

    if group_snv_buf.is_none() {
        return (None, 0);
    }

    let snvs_group = load_group_snvs_buf(group_snv_buf);
    //let group_snvs_use = load_group_snvs_buf(group_snv_buf.unwrap());
    //let snvs_in_group = list_snvs_in_group(&group_snvs_use);
    let m_group = snvs_group.len();

    let group_to_m_in = make_group_to_m_in_vec(&snvs_group, snvs_in);
    // let (group_to_m_in, m_group) = make_group_to_m_in_vec(&snvs_group, snvs_in);
    (Some(group_to_m_in), m_group)
}

pub fn make_group_to_m_in_vec(snvs_group: &[SnvId], snvs_in: &[SnvId]) -> Vec<Vec<usize>> {
    // ) -> (HashMap<usize, Vec<usize>>, usize) {
    let m_in: usize = snvs_in.len();

    // map: sid_in -> plink index
    // &str should be fine
    let mut rs_in_to_index = HashMap::with_capacity(m_in);
    for (si, s) in snvs_in.iter().enumerate() {
        rs_in_to_index.insert(s.ida(), si);
    }

    let mut group_to_m_in: Vec<Vec<usize>> = Vec::with_capacity(snvs_group.len());
    // let mut group_to_m_in: Vec<Vec<usize>> = vec![vec![]; snvs_group.len()];
    // let mut group_to_m_in: HashMap<usize, Vec<usize>> = HashMap::new();
    // let mut m_group_i: usize = 0;
    //let mut m_in_i: usize = 0;

    // for (m_group_i, snv_group) in snvs_group.iter().enumerate() {
    for snv_group in snvs_group.iter() {
        let mut group_to_m_in_i: Vec<usize> = vec![];

        let group_ids_i = snv_group.group_ids().unwrap();
        for snv in group_ids_i.iter() {
            // log::debug!("snv: {:?}", snv);
            // println!("snv: {:?}", snv);
            let ida = snv.ida();
            if ida == "" {
                panic!("ida of SNV is empty. Indicate allele in .group file.")
            }
            match rs_in_to_index.get(ida) {
                Some(&vi) => {
                    group_to_m_in_i.push(vi);

                    //let vi = *v;
                    //let set_id = snv.set_id().unwrap();
                    //match group_to_m_in.get_mut(&set_id) {
                    // match group_to_m_in.get_mut(&m_group_i) {
                    //     Some(v) => {
                    //         v.push(vi);
                    //     }
                    //     None => {
                    //         group_to_m_in.insert(m_group_i, vec![vi]);
                    //     }
                    // }
                }

                // ignore unfound SNVs
                None => {
                    log::info!("SNV in fin_group_snv was not found in plink: {}.", ida);
                }
            }
        }
        group_to_m_in.push(group_to_m_in_i);
        // m_group_i += 1;
    }

    // assert_eq!(group_to_m_in.len(), m_group_i);

    // let m_group = group_to_m_in.len();

    // (group_to_m_in, m_group)
    group_to_m_in
}

pub fn make_agg_to_m_buf(
    agg_snv_buf: Option<&[u8]>,
    snvs: &[SnvId],
) -> (Vec<AggId>, Vec<Vec<usize>>) {
    if agg_snv_buf.is_none() {
        return (vec![], vec![]);
    }

    let mut agg_snvs = load_agg_snvs_buf(agg_snv_buf);
    // log::info!("agg_snvs: {:?}", agg_snvs);

    let agg_to_m = make_agg_to_m_vec(&mut agg_snvs, snvs);

    (agg_snvs, agg_to_m)
}

/// return updated AggId. update snv_ids to SnvId in snvs
pub fn make_agg_to_m_vec(aggs: &mut [AggId], snvs: &[SnvId]) -> Vec<Vec<usize>> {
    let m: usize = snvs.len();

    // map: sid -> snv index
    // &str should be fine
    let mut rs_in_to_index = HashMap::with_capacity(m);
    for (si, s) in snvs.iter().enumerate() {
        rs_in_to_index.insert(s.ida(), (si, s));
    }

    let mut agg_to_m: Vec<Vec<usize>> = vec![];

    // let mut group_to_m_in: HashMap<usize, Vec<usize>> = HashMap::new();
    // let mut m_group_i: usize = 0;

    for agg in aggs.iter_mut() {
        let mut agg_to_m_i: Vec<usize> = vec![];
        let agg_ids = agg.snv_ids();

        let mut snv_ids_update: Vec<SnvId> = vec![];

        for snv in agg_ids.iter() {
            let ida = snv.ida();
            // log::info!("snv: {:?}, {}", snv, ida);
            if ida == "" {
                panic!("ida of SNV is empty. Indicate allele in .agg.snv file.")
            }
            match rs_in_to_index.get(ida) {
                Some(&(vi, snvid)) => {
                    agg_to_m_i.push(vi);
                    snv_ids_update.push(snvid.clone());
                }

                // ignore unfound SNVs
                None => {
                    log::info!("SNV in fin_group_snv was not found in plink: {}.", ida);
                }
            }
        }
        // m_group_i += 1;

        agg.set_snv_ids(snv_ids_update);
        agg_to_m.push(agg_to_m_i);
    }

    // assert_eq!(group_to_m_in.len(), m_group_i);

    agg_to_m
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_make_use_snvs_buf_vec() {
        let snvs_in = vec![
            SnvId::new(
                "rs1".to_owned(),
                "1".to_string(),
                "123",
                "A".to_owned(),
                "C".to_owned(),
            ),
            SnvId::new(
                "rs2".to_owned(),
                "2".to_string(),
                "124",
                "A".to_owned(),
                "C".to_owned(),
            ),
            SnvId::new(
                "rs3".to_owned(),
                "3".to_string(),
                "125",
                "A".to_owned(),
                "C".to_owned(),
            ),
        ];

        let snvs_use = vec![
            SnvId::new(
                "rs3".to_owned(),
                "3".to_string(),
                "125",
                "A".to_owned(),
                "C".to_owned(),
            ),
            SnvId::new(
                "rs1".to_owned(),
                "1".to_string(),
                "123",
                "A".to_owned(),
                "C".to_owned(),
            ),
        ];

        let (use_snvs, m) = make_use_snvs_buf_vec(&snvs_use, &snvs_in);

        assert_eq!(m, 2);
        assert_eq!(use_snvs, vec![true, false, true]);
    }

    #[test]
    fn test_load_snvs_use_buf() {
        let buf: Vec<u8> = "rs1\t1\t123\tA\tC\nrs2\t2\t124\tA\tC\nrs3\t2\t125\tA\tC"
            .to_string()
            .into_bytes();

        let snvs = load_snvs_use_buf(&buf);

        let snvs_ans = vec![
            SnvId::new_id_ma("rs1"),
            SnvId::new_id_ma("rs2"),
            SnvId::new_id_ma("rs3"),
        ];

        assert_eq!(snvs, snvs_ans);
    }

    #[test]
    fn test_load_snvs_use_buf2() {
        let buf: Vec<u8> = "rs1:A>C\nrs2:A>C\nchr1:10:A>C".to_string().into_bytes();

        let snvs = load_snvs_use_buf(&buf);

        let snvs_ans = vec![
            SnvId::new_id_ma("rs1:A>C"),
            SnvId::new_id_ma("rs2:A>C"),
            SnvId::new_id_ma("chr1:10:A>C"),
        ];

        assert_eq!(snvs, snvs_ans);
    }

    #[test]
    fn test_make_group_snvs_buf() {
        let snvs_in = vec![
            SnvId::new(
                "rs1".to_owned(),
                "1".to_string(),
                "123",
                "A".to_owned(),
                "C".to_owned(),
            ),
            SnvId::new(
                "rs2".to_owned(),
                "2".to_string(),
                "124",
                "A".to_owned(),
                "C".to_owned(),
            ),
            SnvId::new(
                "rs3".to_owned(),
                "3".to_string(),
                "125",
                "A".to_owned(),
                "C".to_owned(),
            ),
        ];

        //let set_snvs = vec![
        //    SnvId::new_set_id(
        //        "set3".to_owned(),
        //        vec!["rs3".to_owned(), "rs4".to_owned(), "rs1".to_owned()],
        //    ),
        //    SnvId::new_set_id("set1".to_owned(), vec!["rs1".to_owned()]),
        //];

        let (use_snvs, m) = make_use_snvs_group_buf(None, &snvs_in);
        assert_eq!(m, 0);
        assert_eq!(use_snvs, vec![false, false, false]);

        let (use_snvs, m) = make_use_snvs_group_buf(
            Some(
                "aggid\tchrom\tpos_start\tpos_end\tfilter\tmaf\tnum_snvs\tsnvs\nagg3\t1\t123\t124\tX\tX\t2\trs1:C>A,rs2:C>A\nagg5\t1\t124\t124\tX\tX\t1\trs2:C>A"
                    .as_bytes(),
            ),
            &snvs_in,
        );
        assert_eq!(m, 2);
        assert_eq!(use_snvs, vec![true, true, false]);
    }
}
