//! module for io genotype-related module
// many func. in boosting::predict should be moved here.

use std::collections::HashMap;

pub fn create_m_to_m_in(use_snvs: &[bool]) -> (HashMap<usize, usize>, usize) {
    let mut m_to_m_in = HashMap::new();
    let mut mi = 0;
    for (m_in_i, v) in use_snvs.iter().enumerate() {
        if *v {
            m_to_m_in.insert(mi, m_in_i);
            //m_to_m_in.insert(m_in_i, mi);
            mi += 1;
        }
    }
    (m_to_m_in, mi)
}

pub fn create_m_in_to_m(use_snvs: &[bool]) -> (HashMap<usize, usize>, usize) {
    let mut m_in_to_m = HashMap::new();
    let mut mi = 0;
    for (m_in_i, v) in use_snvs.iter().enumerate() {
        if *v {
            m_in_to_m.insert(m_in_i, mi);
            //m_to_m_in.insert(m_in_i, mi);
            mi += 1;
        }
    }
    (m_in_to_m, mi)
}

// buf_i is the index in use_snvs_buf
// m_i is the index in use_snvs_part
// use_snvs_part; use_snvs_snv or use_snvs_group
pub fn create_m_to_buf(
    use_snvs_buf: &[bool],
    use_snvs_snv: &[bool],
) -> (HashMap<usize, usize>, usize) {
    let mut m_to_buf = HashMap::new();
    let mut buf_i = 0;
    let mut mi = 0;
    for (v, v_snv) in use_snvs_buf.iter().zip(use_snvs_snv.iter()) {
        if *v {
            if *v_snv {
                m_to_buf.insert(mi, buf_i);
                mi += 1;
            }
            buf_i += 1;
        }
    }
    (m_to_buf, mi)
}

// buf_i is the index in use_snvs_buf
// m is the index in use_snvs_part
// use_snvs_part; use_snvs or use_snvs_group
// pub fn create_buf_to_m(
//     use_snvs_buf: &[bool],
//     use_snvs_part: &[bool],
// ) -> (HashMap<usize, usize>, usize) {
//     let mut buf_to_m = HashMap::new();
//     let mut buf_i = 0;
//     let mut mi = 0;
//     for (v, v_part) in use_snvs_buf.iter().zip(use_snvs_part.iter()) {
//         if *v {
//             if *v_part {
//                 buf_to_m.insert(buf_i, mi);
//                 mi += 1;
//             }
//             buf_i += 1;
//         }
//     }
//     (buf_to_m, mi)
// }

pub fn create_m_in_to_buf(use_snvs_buf: &[bool]) -> (HashMap<usize, usize>, usize) {
    create_m_in_to_m(use_snvs_buf)
}

/*
fn create_wgt_to_genotype_index(fin: &str, wgts: &[Wgt]) -> HashMap<usize, Option<usize>> {
    let m_in_wgts = wgts.len();
    let wgt_to_genotype_index = HashMap::with_capacity(m_in_wgts);

    let m_in: usize = plink::compute_num_snv(fin);
    let snvs_in = plink::load_snvs(fin, m_in);

    wgt_to_genotype_index
}
*/

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_create_m_to_m_in() {
        let use_snvs = vec![true, false, true, false, true];
        let (m_to_m_in, m) = create_m_to_m_in(&use_snvs);
        assert_eq!(m, 3);
        assert_eq!(m_to_m_in[&0], 0);
        assert_eq!(m_to_m_in[&1], 2);
        assert_eq!(m_to_m_in[&2], 4);
    }
}
