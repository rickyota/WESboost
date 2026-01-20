//! Group of SNVs
//! Candidate of aggregation of SNVs; input of aggreagation_res
//! Use SnvId for fixed group of SNVs; output of aggregation_res or input of boosting_res

use std::{convert::TryFrom, str::FromStr};

use crate::genot::prelude::*;
use crate::{genot_calc, LdCriteria, SnvId};

// use super::Chrom;

// Now, avoid identifying flipped alelles for multi-length allele
/// reverse: A1 <-> A2
/// flip : in A1, A2, A <-> T
/// TODO: register ref/alt and minor/major alleles
/// all programs uses minor/major but input plink might have different one
#[derive(Clone, Hash, Debug, Default)]
pub struct AggId {
    id: String,
    chrom: String,
    pos_start: usize,
    pos_end: usize,
    // filter: Option<String>,
    // maf: Option<String>,
    // (a1, a2)
    // alleles: Alleles,
    // only for one letter
    //alleles_flip: Alleles,
    // assume (chrom, pos, a1, a2) does not include ":" and sida is unique
    // (chrom):(pos):(a1):(a2)
    // sida: String,
    // (chrom):(pos)
    sid: String,
    // ida: String,
    //group_ids: Option<Vec<String>>,
    // ok?
    //set_id: Option<Vec<SnvId>>,
    //snv_ids: Option<Vec<SnvId>>,
    snv_ids: Vec<SnvId>,
}

// how to implement?
// 1. hash["A"]="T" -> cannot create const hash
// 2. match "A" => "T"
/// assume str is A,C,G,T
//fn flip_allele(a: &str) -> String {
//    match a {
//        "A" => "T".to_owned(),
//        "C" => "G".to_owned(),
//        "G" => "C".to_owned(),
//        "T" => "A".to_owned(),
//        _ => panic!("allele is not one of A,C,G,T"),
//    }
//}

impl AggId {
    // a1, a2 cannot include other than "ACGT".
    //pub fn new_snv_index(
    pub fn new(
        id: String,
        chrom: String,
        pos_start: &str,
        pos_end: &str,
        // a1: String,
        // a2: String,
        // this is ok but need changes all over the codes; instead use in vid()
        //use_snv_pos: bool,
        snv_ids: Vec<SnvId>,
    ) -> Self {
        let mut snv = Self {
            id,
            chrom,
            pos_start: pos_start.parse::<usize>().unwrap(),
            pos_end: pos_end.parse::<usize>().unwrap(),
            // alleles: (a1, a2),
            //alleles_flip: ("".to_owned(), "".to_owned()),
            // sida: "".to_string(),
            sid: "".to_string(),
            // ida: "".to_string(),
            //vid: "".to_string(),
            // group_ids: None,
            snv_ids,
        };
        // TODO: use
        //snv.check_alleles();
        //snv.set_alleles_flip();
        snv.set_id_all();
        //snv.set_sida();
        //snv.set_sid();
        snv
    }

    pub fn new_ma(
        id: String,
        chrom: String,
        pos_start: &str,
        pos_end: &str,
        snv_ids_str: &Vec<&str>,
        // snv_ids_str: &Vec<String>,
    ) -> Self {
        let snv_ids = snv_ids_str.iter().map(|x| SnvId::new_id_ma(x)).collect();
        Self::new(id, chrom, pos_start, pos_end, snv_ids)
    }

    // for use_snvs
    //pub fn new_id(id: String) -> SnvId {
    //    let snv = SnvId {
    //        id,
    //        // dummy for chrom1
    //        chrom: Chrom::try_from(1).unwrap(),
    //        pos: 0,
    //        alleles: ("".to_owned(), "".to_owned()),
    //        //alleles_flip: ("".to_owned(), "".to_owned()),
    //        sida: "".to_string(),
    //        sid: "".to_string(),
    //        group_ids: None,
    //    };
    //    //snv.check_alleles();
    //    //snv.set_alleles_revcomp();
    //    snv.set_sida();
    //    snv
    //}

    // pub fn new_set_ids(id: String, set_ids: Vec<String>) -> SnvId {
    //     let snv = SnvId {
    //         id,
    //         // dummy for chrom1
    //         chrom: Chrom::try_from(1).unwrap(),
    //         pos: 0,
    //         alleles: ("".to_owned(), "".to_owned()),
    //         //alleles_flip: ("".to_owned(), "".to_owned()),
    //         sida: "".to_string(),
    //         sid: "".to_string(),
    //         ida: "".to_string(),
    //         group_ids: Some(set_ids),
    //     };
    //     //snv.check_alleles();
    //     //snv.set_alleles_revcomp();
    //     //snv.set_sida();
    //     snv
    // }

    // TODO: add N etc.
    // TODO: 1kg contained '<DEL>'
    /// should consist of A,C,G,T
    // #[allow(dead_code)]
    // fn check_alleles(&self) {
    //     fn check_allele(a: &str) {
    //         if !(a
    //             .chars()
    //             .all(|v| (v == 'A') || (v == 'C') || (v == 'G') || (v == 'T')))
    //         {
    //             panic!("Alleles should be one of A, C, G, or T: {}.", a);
    //         }
    //     }
    //     check_allele(self.a1());
    //     check_allele(self.a2());
    // }

    //fn set_alleles_flip(&mut self) {
    //    if self.is_one_letter() {
    //        self.alleles_flip = (flip_allele(self.a1()), flip_allele(self.a2()))
    //        //self.alleles_rev = (complement_allele(self.a2()), complement_allele(self.a1()))
    //    }
    //    // else stay ("","")
    //}

    fn set_id_all(&mut self) {
        // self.set_sida();
        self.set_sid();
        // self.set_ida();
    }

    // fn set_sida(&mut self) {
    //     self.sida = self.chrom.to_string()
    //         + ":"
    //         + &self.pos.to_string()
    //         + ":"
    //         + &self.a1()
    //         + ":"
    //         + &self.a2();
    // }

    fn set_sid(&mut self) {
        self.sid = self.chrom.to_string()
            + ":"
            + &self.pos_start.to_string()
            + "-"
            + &self.pos_end.to_string();
    }

    // fn set_ida(&mut self) {
    //     self.ida = self.id.clone() + ":" + &self.a1() + ":" + &self.a2();
    // }

    // pub fn reverse_alleles(&mut self) {
    //     // update alleles, alleles_flip, sida
    //     self.alleles = (self.alleles.1.clone(), self.alleles.0.clone());
    //     //self.set_alleles_flip();
    //     self.set_sida();
    // }

    pub fn id(&self) -> &str {
        &self.id
    }

    // pub fn chrom(&self) -> &Chrom {
    pub fn chrom(&self) -> &str {
        &self.chrom
    }

    pub fn pos_start(&self) -> usize {
        self.pos_start
    }
    pub fn pos_end(&self) -> usize {
        self.pos_end
    }

    // pub fn sida(&self) -> &str {
    //     &self.sida
    // }

    // pub fn ida(&self) -> &str {
    //     &self.ida
    // }

    pub fn sid(&self) -> &str {
        &self.sid
    }

    // /// To compare snvs with or without alelles
    // pub fn vid(&self, use_snv_pos: bool) -> &str {
    //     if use_snv_pos {
    //         self.sid()
    //         //&self.sid
    //     } else {
    //         self.id()
    //         //&self.id
    //     }
    // }

    // fn alleles(&self) -> (&str, &str) {
    //     (&self.alleles.0, &self.alleles.1)
    // }

    // fn alleles_rev(&self) -> (&str, &str) {
    //     (&self.alleles.1, &self.alleles.0)
    // }

    //fn alleles_flip(&self) -> (&str, &str) {
    //    (&self.alleles_flip.0, &self.alleles_flip.1)
    //}

    //fn alleles_rev_flip(&self) -> (&str, &str) {
    //    (&self.alleles_flip.1, &self.alleles_flip.0)
    //}

    // pub fn a1(&self) -> &str {
    //     &self.alleles.0
    // }
    // pub fn a2(&self) -> &str {
    //     &self.alleles.1
    // }
    // flip
    //fn _a1f(&self) -> &str {
    //    &self.alleles_flip.0
    //}
    // flip
    //fn _a2f(&self) -> &str {
    //    &self.alleles_flip.1
    //}

    //pub fn snv_ids(&self) -> Option<&Vec<String>> {
    pub fn snv_ids(&self) -> &Vec<SnvId> {
        &self.snv_ids
    }

    pub fn set_snv_ids(&mut self, snv_ids: Vec<SnvId>) {
        self.snv_ids = snv_ids;
    }

    // pub fn is_alleles_registed(&self) -> bool {
    //     (self.a1() != "") && (self.a2() != "")
    // }

    //pub fn to_sid(&self) -> String {
    //    self.chrom.to_string() + ":" + &self.pos.to_string()
    //}

    //pub fn to_sida_rev(&self) -> String {
    //    self.chrom.to_string() + ":" + &self.pos.to_string() + ":" + &self.a2() + ":" + &self.a1()
    //}
    //// Only for one letter
    //pub fn to_sida_flip(&self) -> Option<String> {
    //    if self.is_one_letter() {
    //        Some(
    //            self.chrom.to_string()
    //                + ":"
    //                + &self.pos.to_string()
    //                + ":"
    //                + &self.a1f()
    //                + ":"
    //                + &self.a2f(),
    //        )
    //    } else {
    //        None
    //    }
    //}
    //// Only for one letter
    //pub fn to_sida_rev_flip(&self) -> Option<String> {
    //    if self.is_one_letter() {
    //        Some(
    //            self.chrom.to_string()
    //                + ":"
    //                + &self.pos.to_string()
    //                + ":"
    //                + &self.a2f()
    //                + ":"
    //                + &self.a1f(),
    //        )
    //    } else {
    //        None
    //    }
    //}

    //fn is_one_letter(&self) -> bool {
    //    (self.a1().len() == 1) && (self.a2().len() == 1)
    //}

    // list all candidates
    //fn flip_or_rev(&self) -> Option<((&str, &str), (&str, &str), (&str, &str), (&str, &str))> {
    // fn flip_or_rev(
    //     &self,
    // ) -> (
    //     (&str, &str),
    //     (&str, &str),
    //     Option<((&str, &str), (&str, &str))>,
    // ) {
    //     let a = self.alleles();
    //     let a_rev = (a.1, a.0);

    //     if !self.is_one_letter() {
    //         return (a, a_rev, None);
    //         //return None;
    //     }
    //     let a_flip = self.alleles_flip();
    //     let a_rev_flip = (a_flip.1, a_flip.0);

    //     (a, a_rev, Some((a_flip, a_rev_flip)))
    //     //Some((a, a_rev, a_flip, a_rev_flip))
    // }

    // use PartialEq trait
    // pub fn is_match()

    // to reverse genotype
    // pub fn is_rev(&self, snv: &SnvId, use_snv_pos: bool) -> bool {
    //     //if self.sid() != snv.sid() {
    //     if self.vid(use_snv_pos) != snv.vid(use_snv_pos) {
    //         return false;
    //     }
    //     //if self.is_one_letter() {
    //     //    return (self.alleles() == snv.alleles_rev())
    //     //        | (self.alleles() == snv.alleles_rev_flip());
    //     //    // otherwise not rev or alleles do not match
    //     //} else {
    //     return self.alleles() == snv.alleles_rev();
    //     //}
    // }

    // pub fn is_in_region(&self, snv_start: &SnvId, snv_end: &SnvId) -> bool {
    //     if snv_start.chrom() != snv_end.chrom() {
    //         panic!("chroms are different");
    //     }

    //     if snv_start.pos() > snv_end.pos() {
    //         panic!("start pos is larger than end pos");
    //     }

    //     if self.chrom() != snv_start.chrom() {
    //         return false;
    //     }

    //     (self.pos() >= snv_start.pos()) && (self.pos() <= snv_end.pos())

    //     //if self.pos() <= snv_start.pos() {
    //     //if self.pos() < snv_start.pos() {
    //     //    return false;
    //     //}
    //     //if self.pos() >= snv_end.pos() {
    //     //if self.pos() > snv_end.pos() {
    //     //    return false;
    //     //}
    //     //return true;
    // }

    //     // TODO: GenotSnv should not be here?
    //     pub fn is_in_ld_criteria(
    //         &self,
    //         snv: &SnvId,
    //         ld_criteria: LdCriteria,
    //         gsnv_self: &GenotSnvRef,
    //         gsnv: &GenotSnvRef,
    //     ) -> bool {
    //         match ld_criteria {
    //             LdCriteria::R2(ld_r2) => self.is_in_ld_r2(snv, ld_r2, gsnv_self, gsnv),
    //             LdCriteria::Radius(ld_radius) => self.is_in_ld_radius(snv, ld_radius),
    //         }
    //     }

    //     fn is_in_ld_radius(&self, snv: &SnvId, ld_radius: usize) -> bool {
    //         if self.chrom() != snv.chrom() {
    //             return false;
    //         }
    //         // usize cannot be minus
    //         if self.pos() > snv.pos() {
    //             return (self.pos() - snv.pos()) <= ld_radius;
    //         } else {
    //             return (snv.pos() - self.pos()) <= ld_radius;
    //         }
    //     }

    //     /// If not on the same Chrom, do not calculate
    //     fn is_in_ld_r2(
    //         &self,
    //         snv: &SnvId,
    //         ld_r2: f64,
    //         gsnv_self: &GenotSnvRef,
    //         gsnv: &GenotSnvRef,
    //     ) -> bool {
    //         if self.chrom() != snv.chrom() {
    //             return false;
    //         }

    //         let r2 = genot_calc::calc_r2(gsnv_self, gsnv);

    //         r2 >= ld_r2

    //         //// usize cannot be minus
    //         //if self.pos() > snv.pos() {
    //         //    return (self.pos() - snv.pos()) <= ld_radius;
    //         //} else {
    //         //    return (snv.pos() - self.pos()) <= ld_radius;
    //         //}
    //     }
}

/*
impl Default for SnvIndex {
    fn default() -> Self {
        Self {
            rs: "".to_owned(),
            chrom: Chrom::Auto(1),
            pos: 0,
            alleles: ("".to_owned(), "".to_owned()),
            alleles_revcomp: ("".to_owned(), "".to_owned()),
            sida: "".to_owned(),
        }
    }
}
*/

// This auto-implement to_string()
impl std::fmt::Display for AggId {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{}", self.id())
        //write!(f, "{}", self.sida())
    }
}

impl PartialEq for AggId {
    fn eq(&self, other: &Self) -> bool {
        self.id() == other.id()
            && self.chrom() == other.chrom()
            && self.pos_start() == other.pos_start()
            && self.pos_end() == other.pos_end()

        // if self.chrom() != other.chrom() || self.pos() != other.pos() {
        //     return false;
        // }

        // check rev for one letter
        //if let Some((a_rev, a_rev_flip)) = alleles_rev {
        //    if (a_rev == a_other) || (a_rev_flip == a_other) {
        //        return true;
        //    }
        //}
        // return false;

        //self.chrom() == other.chrom()
        //    && self.pos() == other.pos()
        //    && self.a1() == other.a1()
        //    && self.a2() == other.a2()

        // if self.chrom() != other.chrom() || self.pos() != other.pos() {
        //     return false;
        // }
        // // same pos
        // let a_other = other.alleles();

        // let (a, a_flip, alleles_rev) = self.flip_or_rev();

        // if (a == a_other) || (a_flip == a_other) {
        //     return true;
        // }

        // // check rev for one letter
        // if let Some((a_rev, a_rev_flip)) = alleles_rev {
        //     if (a_rev == a_other) || (a_rev_flip == a_other) {
        //         return true;
        //     }
        // }
        // return false;

        //match self.flip_or_rev() {
        //    None => self.alleles() == a_other,
        //    Some((a, a_flip, a_rev, a_rev_flip)) => {
        //        (a == a_other)
        //            || (a_flip == a_other)
        //            || (a_rev == a_other)
        //            || (a_rev_flip == a_other)
        //    }
        //}
    }
}

impl Eq for AggId {}

// partial order is the same as order
impl PartialOrd for AggId {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for AggId {
    // TODO: order and eq do not match; is this all right?
    // -> fixed
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        // if eq, return Order.Equal first
        if self.eq(other) {
            return std::cmp::Ordering::Equal;
        }

        //self.id == other.id
        //    && self.chrom == other.chrom
        //    && self.pos_start == other.pos_start
        //    && self.pos_end == other.pos_end

        let ord = self.chrom().cmp(other.chrom());
        if ord.is_ne() {
            return ord;
        }
        let ord = self.pos_start().cmp(&other.pos_end());
        if ord.is_ne() {
            return ord;
        }
        self.id().cmp(other.id())
    }
}

impl AsRef<AggId> for AggId {
    #[inline]
    fn as_ref(&self) -> &AggId {
        self
    }
}

impl AsMut<AggId> for AggId {
    #[inline]
    fn as_mut(&mut self) -> &mut AggId {
        self
    }
}

// later implement FromStr, TryFrom
// from 1:123:A:C or 1_123_A_C
// but ambiguous on A1 and A2

#[cfg(test)]
mod tests {
    use super::*;

    /// test of TryFrom
    #[test]
    fn test_construct_snv_id_string() {
        let snv_id_1 = SnvId::new(
            "rs1".to_owned(),
            "1".to_string(),
            "123",
            "A".to_owned(),
            "C".to_owned(),
        );
        let snv_id_2 = SnvId::new(
            "rs2".to_owned(),
            "1".to_string(),
            "124",
            "A".to_owned(),
            "C".to_owned(),
        );

        let group_id = AggId::new(
            "group1".to_owned(),
            "1".to_string(),
            "123",
            "124",
            vec![snv_id_1, snv_id_2],
        );

        assert_eq!(group_id.id(), "group1");
        assert_eq!(group_id.chrom(), "1");
        assert_eq!(group_id.pos_start(), 123);
        assert_eq!(group_id.pos_end(), 124);
        assert_eq!(group_id.sid(), "1:123-124");
        assert_eq!(group_id.snv_ids.len(), 2);
        assert_eq!(group_id.snv_ids[0].id(), "rs1");
        assert_eq!(group_id.snv_ids[1].id(), "rs2");
    }

    /// test of Display
    #[test]
    fn test_snv_id_to_string() {
        let group_id = AggId::new("group1".to_owned(), "1".to_string(), "123", "124", vec![]);
        assert_eq!(group_id.to_string(), "group1".to_owned());
    }

    #[test]
    fn test_eq() {
        let group_id_1 = AggId::new("group1".to_owned(), "1".to_string(), "123", "124", vec![]);
        let group_id_2 = AggId::new("group1".to_owned(), "1".to_string(), "123", "124", vec![]);

        assert_eq!(group_id_1, group_id_2);

        let group_id_3 = AggId::new("group2".to_owned(), "1".to_string(), "123", "124", vec![]);

        assert_ne!(group_id_1, group_id_3);
    }

    #[test]
    fn test_ord() {
        let group_id_1 = AggId::new("group1".to_owned(), "1".to_string(), "123", "124", vec![]);
        let group_id_2 = AggId::new("group1".to_owned(), "1".to_string(), "124", "125", vec![]);

        assert!(group_id_1 < group_id_2);

        let group_id_3 = AggId::new("group2".to_owned(), "1".to_string(), "123", "124", vec![]);
        assert!(group_id_1 < group_id_3);

        // now avoid identifying flipped alleles
        //// should be ==
        //assert!(!(snv_id_1 < snv_id_3));
        //assert!((snv_id_1 <= snv_id_3));
        //assert!(!(snv_id_1 > snv_id_3));
    }
}
