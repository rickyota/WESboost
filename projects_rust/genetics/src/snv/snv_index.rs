//! SNV index
//! SNV or group of SNVs
//!
//! Now, does not allow flipped alleles for multi-length alleles.
//! id could be duplicated for multi-allelic SNVs.
//!
//! TOFIX: create enum SnvId{Snv, Group}

use crate::genot::prelude::*;
use crate::{genot_calc, LdCriteria};

// use super::Chrom;

type Alleles = (String, String);

// Now, avoid identifying flipped alelles for multi-length allele
/// reverse: A1 <-> A2
/// flip : in A1, A2, A <-> T
/// TODO: register ref/alt and minor/major alleles
/// all programs uses minor/major but input plink might have different one
#[derive(Clone, Hash, Default)]
pub struct SnvId {
    // TOFIX: add 'is_snv' and 'is_agg' to SnvId
    id: String,
    chrom: String,
    pos: usize,
    // (a1, a2)
    // ("", "") only when
    alleles: Alleles,
    // only for one letter
    //alleles_flip: Alleles,
    // assume (chrom, pos, a1, a2) does not include ":" and sida is unique
    // (chrom):(pos):(a1):(a2)
    sida: String,
    // (chrom):(pos):(a2):(a1)
    sidar: String,
    // Better not to use since this confuse multi-allele
    // (chrom):(pos)
    sid: String,
    // (id):(a1):(a2)
    ida: String,
    // (id):(a2):(a1)
    idar: String,
    pos_end: Option<usize>,
    //group_ids: Option<Vec<String>>,
    // ok?
    // could be Some(vec![]) for temporary group when loading wgt on boosting_score
    // later add group_ids using --file-group-snv
    group_ids: Option<Vec<SnvId>>,
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

impl SnvId {
    // a1, a2 cannot include other than "ACGT".
    //pub fn new_snv_index(
    pub fn new(
        id: String,
        chrom: String,
        // chrom: &str,
        pos: &str,
        a1: String,
        a2: String,
        // this is ok but need changes all over the codes; instead use in vid()
        //use_snv_pos: bool,
    ) -> Self {
        let mut snv = Self {
            id,
            chrom,
            // chrom: Chrom::from_str(&chrom).unwrap(),
            pos: pos.parse::<usize>().unwrap(),
            alleles: (a1, a2),
            //alleles_flip: ("".to_owned(), "".to_owned()),
            sida: "".to_string(),
            sidar: "".to_string(),
            sid: "".to_string(),
            ida: "".to_string(),
            idar: "".to_string(),
            pos_end: None,
            group_ids: None,
        };
        // TODO: use
        //snv.check_alleles();
        //snv.set_alleles_flip();
        snv.set_id_all();
        //snv.set_sida();
        //snv.set_sid();
        snv
    }

    // for use_snvs
    pub fn new_id_ma(ida_in: &str) -> Self {
        if ida_in == "None" {
            panic!("wrong implementation.")
        }
        if ida_in.contains(">") {
            // ida = (id):(a2)>(a1)
            let (ida2, a1) = ida_in.rsplit_once(">").unwrap();
            let (id, a2) = ida2.rsplit_once(":").unwrap();

            let mut snv = Self {
                id: id.to_string(),
                chrom: "".to_string(),
                // dummy for chrom1
                // chrom: Chrom::try_from(1).unwrap(),
                pos: 0,
                alleles: (a1.to_string(), a2.to_string()),
                //alleles: ("".to_owned(), "".to_owned()),
                sida: "".to_string(),
                sidar: "".to_string(),
                sid: "".to_string(),
                ida: "".to_string(),
                idar: "".to_string(),
                // group_ids: None,
                pos_end: None,
                group_ids: None,
            };
            snv.set_id_all();
            //snv.set_sida();
            snv
        } else {
            let snv = Self {
                id: ida_in.to_string(),
                chrom: "".to_string(),
                // dummy for chrom1
                // chrom: Chrom::try_from(1).unwrap(),
                pos: 0,
                alleles: ("".to_owned(), "".to_owned()),
                sida: "".to_string(),
                sidar: "".to_string(),
                sid: "".to_string(),
                ida: "".to_string(),
                idar: "".to_string(),
                // group_ids: None,
                pos_end: None,
                group_ids: None,
            };
            //snv.set_sida();
            snv
        }
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

    pub fn new_group(
        id: String,
        chrom: String,
        // chrom: &str,
        pos_start: &str,
        pos_end: &str,
        group_ids: Vec<&str>,
        // group_ids: Vec<String>,
    ) -> Self {
        // println!("pos_start: {}", pos_start);
        let group_ids = if group_ids.len() == 1 && group_ids[0] == "None" {
            vec![]
        } else {
            group_ids.iter().map(|x| SnvId::new_id_ma(x)).collect()
        };
        let mut snv = Self {
            id,
            chrom,
            // chrom: Chrom::from_str(&chrom).unwrap(),
            pos: pos_start.parse::<usize>().unwrap(),
            alleles: ("&".to_owned(), "&".to_owned()),
            sida: "".to_string(),
            sidar: "".to_string(),
            sid: "".to_string(),
            ida: "".to_string(),
            idar: "".to_string(),
            group_ids: Some(group_ids),
            // group_ids: Some(group_ids.iter().map(|x| SnvId::new_id_ma(x)).collect()),
            // group_ids: Some(group_ids),
            pos_end: Some(pos_end.parse::<usize>().unwrap()),
        };
        // log::debug!("snv {:?}", snv);
        snv.set_id_all();
        //snv.check_alleles();
        //snv.set_alleles_revcomp();
        //snv.set_sida();
        snv
    }

    // could be new() or new_group()
    // SnvId::group_ids() = Some(vec![])
    // add group_ids later
    pub fn new_score(
        id: String,
        chrom: String,
        pos: &str,
        // TODO: add pos_end in .wgt
        a1: String,
        a2: String,
    ) -> Self {
        if a1 == "&" && a2 == "&" {
            // add group_ids later
            Self::new_group(id, chrom, pos, pos, vec![])
        } else {
            Self::new(id, chrom, pos, a1, a2)
        }
    }

    // TODO: add N etc.
    // TODO: 1kg contained '<DEL>'
    /// should consist of A,C,G,T
    #[allow(dead_code)]
    fn check_alleles(&self) {
        fn check_allele(a: &str) {
            if !(a
                .chars()
                .all(|v| (v == 'A') || (v == 'C') || (v == 'G') || (v == 'T')))
            {
                panic!("Alleles should be one of A, C, G, or T: {}.", a);
            }
        }
        check_allele(self.a1());
        check_allele(self.a2());
    }

    //fn set_alleles_flip(&mut self) {
    //    if self.is_one_letter() {
    //        self.alleles_flip = (flip_allele(self.a1()), flip_allele(self.a2()))
    //        //self.alleles_rev = (complement_allele(self.a2()), complement_allele(self.a1()))
    //    }
    //    // else stay ("","")
    //}

    fn set_id_all(&mut self) {
        self.set_sida();
        self.set_sidar();
        self.set_sid();
        self.set_ida();
        self.set_idar();
    }

    fn set_sida(&mut self) {
        self.sida = self.chrom.to_string()
            + ":"
            + &self.pos.to_string()
            + ":"
            + &self.a1()
            + ":"
            + &self.a2();
    }

    fn set_sidar(&mut self) {
        self.sidar = self.chrom.to_string()
            + ":"
            + &self.pos.to_string()
            + ":"
            + &self.a2()
            + ":"
            + &self.a1();
    }

    fn set_sid(&mut self) {
        self.sid = self.chrom.to_string() + ":" + &self.pos.to_string();
    }

    fn set_ida(&mut self) {
        self.ida = self.id.clone() + ":" + &self.a1() + ":" + &self.a2();
    }

    fn set_idar(&mut self) {
        self.idar = self.id.clone() + ":" + &self.a2() + ":" + &self.a1();
    }

    pub fn reverse_alleles(&mut self) {
        // update alleles, alleles_flip, sida
        self.alleles = (self.alleles.1.clone(), self.alleles.0.clone());
        //self.set_alleles_flip();
        // self.set_sida();
        // assume sida~idar is already set.
        self.set_id_all();
    }

    pub fn id(&self) -> &str {
        &self.id
    }

    // pub fn chrom(&self) -> &Chrom {
    pub fn chrom(&self) -> &str {
        &self.chrom
    }

    pub fn pos(&self) -> usize {
        self.pos
    }

    pub fn sida(&self) -> &str {
        &self.sida
    }

    pub fn sidar(&self) -> &str {
        &self.sidar
    }

    pub fn sid(&self) -> &str {
        &self.sid
    }

    pub fn ida(&self) -> &str {
        &self.ida
    }

    pub fn idar(&self) -> &str {
        &self.idar
    }

    /// To compare snvs with or without alelles
    // pub fn vid(&self, use_snv_pos: bool) -> &str {
    //     panic!("deprecated: use vida() and vidar()");
    //     if use_snv_pos {
    //         self.sid()
    //         //&self.sid
    //     } else {
    //         self.id()
    //         //&self.id
    //     }
    // }

    pub fn vida(&self, use_snv_pos: bool) -> &str {
        if use_snv_pos {
            self.sida()
        } else {
            self.ida()
        }
    }

    pub fn vidar(&self, use_snv_pos: bool) -> &str {
        if use_snv_pos {
            self.sidar()
        } else {
            self.idar()
        }
    }

    // TODO: add to struct
    // (id):(a2)>(a1)
    pub fn idma(&self) -> String {
        self.id.clone() + ":" + &self.a2() + ">" + &self.a1()
    }

    fn alleles(&self) -> (&str, &str) {
        (&self.alleles.0, &self.alleles.1)
    }

    fn alleles_rev(&self) -> (&str, &str) {
        (&self.alleles.1, &self.alleles.0)
    }

    //fn alleles_flip(&self) -> (&str, &str) {
    //    (&self.alleles_flip.0, &self.alleles_flip.1)
    //}

    //fn alleles_rev_flip(&self) -> (&str, &str) {
    //    (&self.alleles_flip.1, &self.alleles_flip.0)
    //}

    pub fn a1(&self) -> &str {
        &self.alleles.0
    }
    pub fn a2(&self) -> &str {
        &self.alleles.1
    }
    // flip
    //fn _a1f(&self) -> &str {
    //    &self.alleles_flip.0
    //}
    // flip
    //fn _a2f(&self) -> &str {
    //    &self.alleles_flip.1
    //}

    //pub fn group_ids(&self) -> Option<&Vec<String>> {
    pub fn group_ids(&self) -> Option<&Vec<SnvId>> {
        self.group_ids.as_ref()
    }

    pub fn is_group(&self) -> bool {
        self.group_ids().is_some()
    }

    pub fn is_alleles_registered(&self) -> bool {
        (self.a1() != "") && (self.a2() != "")
    }

    pub fn eq_vid(&self, snv: &SnvId, use_snv_pos: bool) -> bool {
        self.vida(use_snv_pos) == snv.vida(use_snv_pos)
            || self.vida(use_snv_pos) == snv.vidar(use_snv_pos)
    }

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
    pub fn is_rev(&self, snv: &SnvId, use_snv_pos: bool) -> Option<bool> {
        //if self.sid() != snv.sid() {
        // if self.vid(use_snv_pos) != snv.vid(use_snv_pos) {
        if !self.eq_vid(snv, use_snv_pos) {
            return None;
        }
        return Some(self.alleles() == snv.alleles_rev());
    }

    pub fn is_in_region(&self, snv_start: &SnvId, snv_end: &SnvId) -> bool {
        if snv_start.chrom() != snv_end.chrom() {
            panic!("chroms are different");
        }

        if snv_start.pos() > snv_end.pos() {
            panic!("start pos is larger than end pos");
        }

        if self.chrom() != snv_start.chrom() {
            return false;
        }

        (self.pos() >= snv_start.pos()) && (self.pos() <= snv_end.pos())

        //if self.pos() <= snv_start.pos() {
        //if self.pos() < snv_start.pos() {
        //    return false;
        //}
        //if self.pos() >= snv_end.pos() {
        //if self.pos() > snv_end.pos() {
        //    return false;
        //}
        //return true;
    }

    // TODO: GenotSnv should not be here?
    pub fn is_in_ld_criteria(
        &self,
        snv: &SnvId,
        ld_criteria: LdCriteria,
        gsnv_self: &GenotSnvRef,
        gsnv: &GenotSnvRef,
    ) -> bool {
        match ld_criteria {
            LdCriteria::R2(ld_r2) => self.is_in_ld_r2(snv, ld_r2, gsnv_self, gsnv),
            LdCriteria::Radius(ld_radius) => self.is_in_ld_radius(snv, ld_radius),
        }
    }

    fn is_in_ld_radius(&self, snv: &SnvId, ld_radius: usize) -> bool {
        if self.chrom() != snv.chrom() {
            return false;
        }
        // usize cannot be minus
        if self.pos() > snv.pos() {
            return (self.pos() - snv.pos()) <= ld_radius;
        } else {
            return (snv.pos() - self.pos()) <= ld_radius;
        }
    }

    /// If not on the same Chrom, do not calculate
    fn is_in_ld_r2(
        &self,
        snv: &SnvId,
        ld_r2: f64,
        gsnv_self: &GenotSnvRef,
        gsnv: &GenotSnvRef,
    ) -> bool {
        if self.chrom() != snv.chrom() {
            return false;
        }

        let r2 = genot_calc::calc_r2(gsnv_self, gsnv);

        r2 >= ld_r2

        //// usize cannot be minus
        //if self.pos() > snv.pos() {
        //    return (self.pos() - snv.pos()) <= ld_radius;
        //} else {
        //    return (snv.pos() - self.pos()) <= ld_radius;
        //}
    }
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
impl std::fmt::Display for SnvId {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{}", self.sida())
    }
}

// Printing all group_ids may be too long
impl std::fmt::Debug for SnvId {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(
            f,
            "SnvId {{ id: {}, chrom: {}, pos: {}, a1: {}, a2: {}, sida: {}, sid: {}, pos_end: {:?}, group_ids: {} }}",
            self.id, self.chrom, self.pos, self.a1(), self.a2(), self.sida, self.sid,self.pos_end, self.group_ids.as_ref().map_or("None".to_string(), |g| format!("{:?} snvs", g.len() ) )
        )

        // write!(
        //     f,
        //     "SnvId {{ id: {}, chrom: {}, pos: {}, a1: {}, a2: {}, sida: {}, sidar: {}, sid: {}, ida: {}, idar: {} }}",
        //     self.id, self.chrom, self.pos, self.a1(), self.a2(), self.sida, self.sidar, self.sid, self.ida, self.idar
        // )
    }
}

impl PartialEq for SnvId {
    /// For snv comparison including reversed allele, use eq_vid()
    /// mainly for sorting
    fn eq(&self, other: &Self) -> bool {
        self.chrom() == other.chrom()
            && self.pos() == other.pos()
            && self.a1() == other.a1()
            && self.a2() == other.a2()

        // if self.chrom() != other.chrom() || self.pos() != other.pos() {
        //     return false;
        // }

        // // same pos
        // let a_other = other.alleles();

        // //let (a, a_flip, alleles_rev) = self.flip_or_rev();
        // let a = self.alleles();
        // let a_rev = self.alleles_rev();

        // if (a == a_other) || (a_rev == a_other) {
        //     return true;
        // }

        // // check rev for one letter
        // //if let Some((a_rev, a_rev_flip)) = alleles_rev {
        // //    if (a_rev == a_other) || (a_rev_flip == a_other) {
        // //        return true;
        // //    }
        // //}
        // return false;

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

impl Eq for SnvId {}

// partial order is the same as order
impl PartialOrd for SnvId {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for SnvId {
    // not considering reversed alleles for sort consistency
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        let ord = self.chrom().cmp(other.chrom());
        if ord.is_ne() {
            return ord;
        }
        let ord = self.pos().cmp(&other.pos());
        if ord.is_ne() {
            return ord;
        }
        let ord = self.a1().cmp(other.a1());
        if ord.is_ne() {
            return ord;
        }
        self.a2().cmp(other.a2())
    }

    // TOFIX: order and eq do not match; is this all right?
    // rs1:A:T = rs1:T:A < rs1:A:C
    // fn cmp(&self, other: &Self) -> std::cmp::Ordering {
    //     // if eq, return Order.Equal first
    //     if self.eq(other) {
    //         return std::cmp::Ordering::Equal;
    //     }

    //     let ord = self.chrom().cmp(other.chrom());
    //     if ord.is_ne() {
    //         return ord;
    //     }
    //     let ord = self.pos().cmp(&other.pos());
    //     if ord.is_ne() {
    //         return ord;
    //     }
    //     let ord = self.a1().cmp(other.a1());
    //     if ord.is_ne() {
    //         return ord;
    //     }
    //     self.a2().cmp(other.a2())
    // }
}

impl AsRef<SnvId> for SnvId {
    #[inline]
    fn as_ref(&self) -> &SnvId {
        self
    }
}

impl AsMut<SnvId> for SnvId {
    #[inline]
    fn as_mut(&mut self) -> &mut SnvId {
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
        let snv_id = SnvId::new(
            "rs1".to_owned(),
            "1".to_string(),
            "123",
            "A".to_owned(),
            "C".to_owned(),
        );
        assert_eq!(snv_id.id(), "rs1");
        assert_eq!(snv_id.chrom(), "1");
        // assert_eq!(snv_id.chrom(), &Chrom::Auto(1));
        assert_eq!(snv_id.pos(), 123);
        assert_eq!(snv_id.a1(), "A");
        assert_eq!(snv_id.a2(), "C");
        assert_eq!(snv_id.sida(), "1:123:A:C");
        assert_eq!(snv_id.sid(), "1:123");
        assert_eq!(snv_id.ida(), "rs1:A:C");
    }

    #[test]
    fn test_new_id_ma() {
        let snv_id = SnvId::new_id_ma("rs1:1:1:A>C");
        assert_eq!(snv_id.id(), "rs1:1:1");
        // assert_eq!(snv_id.chrom(), "1");
        // assert_eq!(snv_id.chrom(), &Chrom::Auto(1));
        // assert_eq!(snv_id.pos(), 123);
        assert_eq!(snv_id.a1(), "C");
        assert_eq!(snv_id.a2(), "A");
        // assert_eq!(snv_id.sida(), "1:123:A:C");
        // assert_eq!(snv_id.sid(), "1:123");
        assert_eq!(snv_id.ida(), "rs1:1:1:C:A");
        assert_eq!(snv_id.idar(), "rs1:1:1:A:C");
    }

    #[test]
    #[should_panic]
    fn test_new_id_ma_panic() {
        let snv_id = SnvId::new_id_ma("None");
    }

    // #[test]
    // #[should_panic]
    // fn test_construct_snv_id_string_panic() {
    //     let _ = SnvId::new("rs1".to_owned(), "1", "123", "N".to_owned(), "N".to_owned());
    // }

    // #[test]
    // fn test_rev_flip() {
    //     let snv_id = SnvId::new("rs1".to_owned(), "1", "123", "A".to_owned(), "C".to_owned());
    //     //assert_eq!(snv_id.alleles_rev(), ("C", "A"));
    //     //assert_eq!(snv_id.alleles_flip(), ("T", "G"));
    //     //assert_eq!(snv_id.alleles_rev_flip(), ("G", "T"));
    // }

    // #[test]
    // fn test_rev_flip_long() {
    //     let snv_id = SnvId::new(
    //         "rs1".to_owned(),
    //         "1",
    //         "123",
    //         "AAT".to_owned(),
    //         "C".to_owned(),
    //     );
    //     assert_eq!(snv_id.alleles_rev(), ("C", "AAT"));
    //     assert_eq!(snv_id.alleles_flip(), ("", ""));
    //     assert_eq!(snv_id.alleles_rev_flip(), ("", ""));
    // }

    // #[test]
    // fn test_flip_or_rev() {
    //     let snv_id = SnvId::new("rs1".to_owned(), "1", "123", "A".to_owned(), "C".to_owned());
    //     assert_eq!(snv_id.id(), "rs1");
    //     assert_eq!(snv_id.chrom(), &Chrom::Auto(1));
    //     assert_eq!(snv_id.pos(), 123);
    //     assert_eq!(snv_id.a1(), "A");
    //     assert_eq!(snv_id.a2(), "C");

    //     let alleles = snv_id.flip_or_rev();
    //     assert_eq!(alleles.0, ("A", "C"));
    //     assert_eq!(alleles.1, ("C", "A"));
    //     assert_eq!(alleles.2.unwrap().0, ("T", "G"));
    //     assert_eq!(alleles.2.unwrap().1, ("G", "T"));
    // }

    // #[test]
    // fn test_flip_or_rev_long() {
    //     let snv_id = SnvId::new(
    //         "rs1".to_owned(),
    //         "1",
    //         "123",
    //         "AAT".to_owned(),
    //         "C".to_owned(),
    //     );

    //     let alleles = snv_id.flip_or_rev();
    //     assert_eq!(alleles.0, ("AAT", "C"));
    //     assert_eq!(alleles.1, ("C", "AAT"));
    //     assert_eq!(alleles.2, None);
    // }

    // #[test]
    // fn test_is_rev() {
    //     let snv_id_1 = SnvId::new("rs1".to_owned(), "1", "123", "A".to_owned(), "C".to_owned());
    //     let snv_id_2 = SnvId::new("rs1".to_owned(), "1", "123", "A".to_owned(), "C".to_owned());
    //     assert!(!snv_id_1.is_rev(&snv_id_2, false));

    //     let snv_id_2 = SnvId::new("rs1".to_owned(), "1", "123", "T".to_owned(), "G".to_owned());
    //     assert!(!snv_id_1.is_rev(&snv_id_2, false));

    //     let snv_id_2 = SnvId::new("rs1".to_owned(), "1", "123", "C".to_owned(), "A".to_owned());
    //     assert!(snv_id_1.is_rev(&snv_id_2, false));

    //     let snv_id_2 = SnvId::new("rs1".to_owned(), "1", "123", "G".to_owned(), "T".to_owned());
    //     assert!(snv_id_1.is_rev(&snv_id_2, false));

    //     // alleles do not match
    //     let snv_id_2 = SnvId::new("rs1".to_owned(), "1", "123", "A".to_owned(), "G".to_owned());
    //     assert!(!snv_id_1.is_rev(&snv_id_2, false));

    //     // snv does not match
    //     let snv_id_2 = SnvId::new("rs2".to_owned(), "1", "124", "C".to_owned(), "A".to_owned());
    //     assert!(!snv_id_1.is_rev(&snv_id_2, false));
    // }

    #[test]
    fn test_is_in_region() {
        let snv_id = SnvId::new(
            "rs1".to_owned(),
            "1".to_string(),
            "100",
            "A".to_owned(),
            "C".to_owned(),
        );

        let inputs = [
            ((1, 50, 150), true),   // in
            ((1, 100, 150), true),  // on the border
            ((1, 120, 150), false), // out
            ((1, 50, 80), false),   // out
            ((2, 100, 150), false), // chrom is different
        ];

        for ((chrom, pos_start, pos_end), exp) in inputs {
            let snv_start = SnvId::new(
                "rs1".to_owned(),
                chrom.to_string(),
                &pos_start.to_string(),
                "A".to_owned(),
                "C".to_owned(),
            );
            let snv_end = SnvId::new(
                "rs1".to_owned(),
                chrom.to_string(),
                &pos_end.to_string(),
                "A".to_owned(),
                "C".to_owned(),
            );
            assert_eq!(snv_id.is_in_region(&snv_start, &snv_end), exp);
        }
    }

    //#[test]
    //fn test_is_in_region() {
    //    let snv_id = SnvId::new("rs1".to_owned(), "1", "100", "A".to_owned(), "C".to_owned());

    //    let snv_start = SnvId::new("rs1".to_owned(), "1", "50", "A".to_owned(), "C".to_owned());
    //    let snv_end = SnvId::new("rs1".to_owned(), "1", "150", "A".to_owned(), "C".to_owned());

    //    assert!(snv_id.is_in_region(&snv_start, &snv_end));
    //}

    //#[test]
    //fn test_is_in_region_2() {
    //    // when snv_start == snv_end
    //    let snv_id = SnvId::new("rs1".to_owned(), "1", "100", "A".to_owned(), "C".to_owned());

    //    let snv_start = SnvId::new("rs1".to_owned(), "1", "100", "A".to_owned(), "C".to_owned());
    //    let snv_end = SnvId::new("rs1".to_owned(), "1", "150", "A".to_owned(), "C".to_owned());

    //    assert!(snv_id.is_in_region(&snv_start, &snv_end));
    //}

    //#[test]
    //fn test_is_in_region_3() {
    //    // outside
    //    let snv_id = SnvId::new("rs1".to_owned(), "1", "100", "A".to_owned(), "C".to_owned());

    //    let snv_start = SnvId::new("rs1".to_owned(), "1", "120", "A".to_owned(), "C".to_owned());
    //    let snv_end = SnvId::new("rs1".to_owned(), "1", "150", "A".to_owned(), "C".to_owned());

    //    assert!(!snv_id.is_in_region(&snv_start, &snv_end));
    //}

    //#[test]
    //fn test_is_in_region_4() {
    //    // outside
    //    let snv_id = SnvId::new("rs1".to_owned(), "1", "100", "A".to_owned(), "C".to_owned());

    //    let snv_start = SnvId::new("rs1".to_owned(), "1", "50", "A".to_owned(), "C".to_owned());
    //    let snv_end = SnvId::new("rs1".to_owned(), "1", "80", "A".to_owned(), "C".to_owned());

    //    assert!(!snv_id.is_in_region(&snv_start, &snv_end));
    //}

    //#[test]
    //fn test_is_in_region_5() {
    //    // chrom is different
    //    let snv_id = SnvId::new("rs1".to_owned(), "1", "100", "A".to_owned(), "C".to_owned());

    //    let snv_start = SnvId::new("rs1".to_owned(), "2", "100", "A".to_owned(), "C".to_owned());
    //    let snv_end = SnvId::new("rs1".to_owned(), "2", "150", "A".to_owned(), "C".to_owned());

    //    assert!(!snv_id.is_in_region(&snv_start, &snv_end));
    //}

    #[test]
    #[should_panic]
    fn test_is_in_region_panic_1() {
        // when chrom of snv_start != snv_end
        let snv_id = SnvId::new(
            "rs1".to_owned(),
            "1".to_string(),
            "100",
            "A".to_owned(),
            "C".to_owned(),
        );

        let snv_start = SnvId::new(
            "rs1".to_owned(),
            "1".to_string(),
            "50",
            "A".to_owned(),
            "C".to_owned(),
        );
        let snv_end = SnvId::new(
            "rs1".to_owned(),
            "2".to_string(),
            "150",
            "A".to_owned(),
            "C".to_owned(),
        );

        assert!(snv_id.is_in_region(&snv_start, &snv_end));
    }

    #[test]
    #[should_panic]
    fn test_is_in_region_panic_2() {
        // when pos of start < end
        let snv_id = SnvId::new(
            "rs1".to_owned(),
            "1".to_string(),
            "100",
            "A".to_owned(),
            "C".to_owned(),
        );

        let snv_start = SnvId::new(
            "rs1".to_owned(),
            "1".to_string(),
            "150",
            "A".to_owned(),
            "C".to_owned(),
        );
        let snv_end = SnvId::new(
            "rs1".to_owned(),
            "1".to_string(),
            "50",
            "A".to_owned(),
            "C".to_owned(),
        );

        assert!(snv_id.is_in_region(&snv_start, &snv_end));
    }

    #[test]
    fn test_is_in_ld_radius() {
        let snv_id = SnvId::new(
            "rs1".to_owned(),
            "1".to_string(),
            "100",
            "A".to_owned(),
            "C".to_owned(),
        );

        let ldr = 50;

        let inputs = [
            ((1, 70), true),   // in
            ((1, 120), true),  // in
            ((1, 50), true),   // on the border
            ((1, 20), false),  // out
            ((1, 180), false), // out
            ((2, 100), false), // chrom is different
        ];

        for ((chrom, pos_start), exp) in inputs {
            let snv_start = SnvId::new(
                "rs1".to_owned(),
                chrom.to_string(),
                &pos_start.to_string(),
                "A".to_owned(),
                "C".to_owned(),
            );
            assert_eq!(snv_id.is_in_ld_radius(&snv_start, ldr), exp);
        }
    }

    /// test of Display
    #[test]
    fn test_snv_id_to_string() {
        let snv_id = SnvId::new(
            "rs1".to_owned(),
            "1".to_string(),
            "123",
            "A".to_owned(),
            "C".to_owned(),
        );
        assert_eq!(snv_id.to_string(), "1:123:A:C".to_owned());
    }

    #[test]
    fn test_check_alleles() {
        let snv_id = SnvId::new(
            "rs1".to_owned(),
            "1".to_string(),
            "123",
            "A".to_owned(),
            "C".to_owned(),
        );
        snv_id.check_alleles()
    }

    #[test]
    fn test_is_rev() {
        let snv_id_1 = SnvId::new(
            "rs1".to_owned(),
            "1".to_string(),
            "123",
            "A".to_owned(),
            "C".to_owned(),
        );
        let snv_id_2 = SnvId::new(
            "rs1".to_owned(),
            "1".to_string(),
            "123",
            "A".to_owned(),
            "C".to_owned(),
        );
        let snv_id_2_2 = SnvId::new(
            "rs1".to_owned(),
            "1".to_string(),
            "123",
            "C".to_owned(),
            "A".to_owned(),
        );
        let snv_id_2_3 = SnvId::new(
            "rs1".to_owned(),
            "1".to_string(),
            "123",
            "A".to_owned(),
            "T".to_owned(),
        );

        assert_eq!(snv_id_1.is_rev(&snv_id_2, false), Some(false));
        assert_eq!(snv_id_1.is_rev(&snv_id_2_2, false), Some(true));
        assert_eq!(snv_id_1.is_rev(&snv_id_2_3, false), None);

        let snv_id_3 = SnvId::new(
            "1:123".to_owned(),
            "1".to_string(),
            "123",
            "A".to_owned(),
            "C".to_owned(),
        );
        let snv_id_3_2 = SnvId::new(
            "1:123".to_owned(),
            "1".to_string(),
            "123",
            "C".to_owned(),
            "A".to_owned(),
        );
        let snv_id_3_3 = SnvId::new(
            "1:123".to_owned(),
            "1".to_string(),
            "123",
            "A".to_owned(),
            "T".to_owned(),
        );

        assert_eq!(snv_id_1.is_rev(&snv_id_3, true), Some(false));
        assert_eq!(snv_id_1.is_rev(&snv_id_3_2, true), Some(true));
        assert_eq!(snv_id_1.is_rev(&snv_id_3_3, true), None);
    }

    #[test]
    fn test_eq_vid() {
        let snv_id_1 = SnvId::new(
            "rs1".to_owned(),
            "1".to_string(),
            "123",
            "A".to_owned(),
            "C".to_owned(),
        );
        let snv_id_2 = SnvId::new(
            "rs1".to_owned(),
            "1".to_string(),
            "123",
            "A".to_owned(),
            "C".to_owned(),
        );
        let snv_id_2_2 = SnvId::new(
            "rs1".to_owned(),
            "1".to_string(),
            "123",
            "C".to_owned(),
            "A".to_owned(),
        );
        let snv_id_2_3 = SnvId::new(
            "rs1".to_owned(),
            "1".to_string(),
            "123",
            "A".to_owned(),
            "T".to_owned(),
        );

        assert!(snv_id_1.eq_vid(&snv_id_2, false));
        assert!(snv_id_1.eq_vid(&snv_id_2_2, false));
        assert!(!snv_id_1.eq_vid(&snv_id_2_3, false));

        let snv_id_3 = SnvId::new(
            "1:123".to_owned(),
            "1".to_string(),
            "123",
            "A".to_owned(),
            "C".to_owned(),
        );
        let snv_id_3_2 = SnvId::new(
            "1:123".to_owned(),
            "1".to_string(),
            "123",
            "C".to_owned(),
            "A".to_owned(),
        );
        let snv_id_3_3 = SnvId::new(
            "1:123".to_owned(),
            "1".to_string(),
            "123",
            "A".to_owned(),
            "T".to_owned(),
        );

        assert!(snv_id_1.eq_vid(&snv_id_3, true));
        assert!(snv_id_1.eq_vid(&snv_id_3_2, true));
        assert!(!snv_id_1.eq_vid(&snv_id_3_3, true));
    }

    #[test]
    fn test_eq() {
        let snv_id_1 = SnvId::new(
            "rs1".to_owned(),
            "1".to_string(),
            "123",
            "A".to_owned(),
            "C".to_owned(),
        );
        let snv_id_2 = SnvId::new(
            "1:123:G:T".to_owned(),
            "1".to_string(),
            "123",
            "A".to_owned(),
            "C".to_owned(),
        );

        assert_eq!(snv_id_1, snv_id_2);

        let snv_id_3 = SnvId::new(
            "rs1".to_owned(),
            "1".to_string(),
            "123",
            "C".to_owned(),
            "A".to_owned(),
        );

        assert_ne!(snv_id_1, snv_id_3);

        let snv_id_4 = SnvId::new(
            "rs1".to_owned(),
            "1".to_string(),
            "123",
            "A".to_owned(),
            "G".to_owned(),
        );

        assert_ne!(snv_id_1, snv_id_4);
    }

    #[test]
    fn test_ord() {
        let snv_id_1 = SnvId::new(
            "rs1".to_owned(),
            "1".to_string(),
            "123",
            "A".to_owned(),
            "C".to_owned(),
        );
        let snv_id_2 = SnvId::new(
            "rs1".to_owned(),
            "1".to_string(),
            "123",
            "A".to_owned(),
            "G".to_owned(),
        );

        assert!(snv_id_1 < snv_id_2);

        let snv_id_3 = SnvId::new(
            "rs1".to_owned(),
            "1".to_string(),
            "123",
            "G".to_owned(),
            "T".to_owned(),
        );
        assert!(snv_id_1 < snv_id_3);

        // now avoid identifying flipped alleles
        //// should be ==
        //assert!(!(snv_id_1 < snv_id_3));
        //assert!((snv_id_1 <= snv_id_3));
        //assert!(!(snv_id_1 > snv_id_3));
    }

    // TODO: test Debug for SnvId

    #[test]
    fn test_debug() {
        let snv_id = SnvId::new(
            "rs1".to_owned(),
            "1".to_string(),
            "123",
            "A".to_owned(),
            "C".to_owned(),
        );
        let out: String = format!("{:?}", snv_id);
        assert_eq!(
            out,
            "SnvId { id: rs1, chrom: 1, pos: 123, a1: A, a2: C, sida: 1:123:A:C, sid: 1:123, pos_end: None, group_ids: None }"
        );
    }

    #[test]
    fn test_debug_group() {
        let snv_id = SnvId::new_group(
            "rs1".to_owned(),
            "1".to_string(),
            "123",
            "456",
            vec!["rs2", "rs3"],
        );

        let out: String = format!("{:?}", snv_id);
        assert_eq!(
            out,
            "SnvId { id: rs1, chrom: 1, pos: 123, a1: &, a2: &, sida: 1:123:&:&, sid: 1:123, pos_end: Some(456), group_ids: 2 snvs }"
        );
    }
}
