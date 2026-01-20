//! Base trait for CBits family
//!
//! To avoid name collision to base_cvec, use _b
//!

use super::iterator::BoolIter;
use super::{calc, CBits};
use super::{BaseCMatrix, CBitsRef, B8};
use super::{BaseCVec, BaseCVecMut};

/// for CBits, CBitsRef, CBitsMut
pub trait BaseCBits: BaseCVec {
    fn clone_to_cbits(&self) -> CBits {
        CBits::new_from_cmatrix(self.clone_to_cmatrix())
    }

    #[inline]
    fn get_unchecked_b(&self, col_i: usize) -> bool {
        self.get_unchecked(0, col_i, 0)
    }

    #[inline]
    fn get_b(&self, col_i: usize) -> bool {
        self.get(0, col_i, 0)
    }

    #[inline]
    fn inner_vals_col_digit_b(&self) -> Vec<bool> {
        self.inner_vals_col_digit(0, 0)
    }

    #[allow(unreachable_code)]
    fn count(&self) -> usize {
        #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
        {
            return count_simd(self.inner(), self.col_n());
        }
        count_nosimd(self.inner(), self.col_n())
    }

    fn count_false(&self) -> usize {
        self.col_n() - self.count()
    }

    fn as_cbits_ref_b(&self) -> CBitsRef {
        CBitsRef::new(self.inner(), self.col_n())
    }

    fn iter(&self) -> BoolIter {
        BoolIter::new(self.as_cbits_ref_b())
    }

    // move to BaseGenotBits
    fn stat_contingency_table(&self, ys: &[B8], case_n: usize) -> (usize, usize, usize, usize) {
        //use std::time::Instant;
        //let start_time = Instant::now();

        let n = self.col_n();

        // assert_eq!(self.n(), phe.n());

        //println!("in table afr init: {} sec",  start_time.elapsed().as_micros());

        let pred_s0m = self.inner();

        // TODO: counting both and substract case could be faster
        fn sum_byte(y: u32, p0: u32) -> (usize, usize) {
            // do not use d0/n0 since could be mixed with padding
            let d1 = y & p0;
            let n1 = (!y) & p0;
            (crate::popcnt(d1), crate::popcnt(n1))
        }

        fn add_tuple2(sums: (usize, usize), sums_: (usize, usize)) -> (usize, usize) {
            (sums.0 + sums_.0, sums.1 + sums_.1)
        }

        let mut sums: (usize, usize) = (0usize, 0);
        for ni in 0..(n / 32 + 1) {
            let pred_s0_b32 =
                u32::from_le_bytes(pred_s0m[4 * ni..4 * (ni + 1)].try_into().unwrap());
            let ys_b32 = u32::from_le_bytes(ys[4 * ni..4 * (ni + 1)].try_into().unwrap());

            let sums_32 = sum_byte(ys_b32, pred_s0_b32);

            sums = add_tuple2(sums, sums_32);
        }

        let d1 = sums.0;
        let n1 = sums.1;

        let d0 = case_n - d1;
        let n0 = n - case_n - n1;

        // let dall = phe.count();
        // let nall = n - dall;
        // let d0 = dall - (d2 + d1 + dm);
        // let n0 = nall - (n2 + n1 + nm);
        //println!("in table afr last: {} sec",  start_time.elapsed().as_micros());

        // (d2, n2, d1, n1, d0, n0, dm, nm)
        (d1, n1, d0, n0)
    }

    // move to BaseGenotBits
    fn maf_group(&self) -> f64 {
        self.count() as f64 / self.col_n() as f64
    }
}

#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
fn count_simd(ys: &[u8], n: usize) -> usize {
    let mut count = 0;

    // padding should be false
    for ni in 0..(n / 32 + 1) {
        let ys_b32 = u32::from_le_bytes(ys[4 * ni..4 * (ni + 1)].try_into().unwrap());
        count = count + calc::popcnt_simd(ys_b32);
    }

    count
}

fn count_nosimd(ys: &[u8], n: usize) -> usize {
    let mut count = 0;

    // padding should be false
    for ni in 0..(n / 32 + 1) {
        let ys_b32 = u32::from_le_bytes(ys[4 * ni..4 * (ni + 1)].try_into().unwrap());
        count = count + calc::popcnt_nosimd(ys_b32);
    }

    count
}

/// for CBits, CBitsMut
pub trait BaseCBitsMut: BaseCVecMut + BaseCBits {
    fn set_bool_unchecked_b(&mut self, b: bool, col_i: usize) {
        self.set_bool_unchecked(b, 0, col_i, 0);
    }

    fn set_bool_b(&mut self, b: bool, col_i: usize) {
        self.set_bool(b, 0, col_i, 0);
    }

    fn or_bitwise(&mut self, other: CBitsRef) {
        let v = self.inner_mut();
        let v_other = other.inner();
        v.iter_mut()
            .zip(v_other.iter())
            .for_each(|(b1, b2)| *b1 = *b1 | *b2);
    }
}
