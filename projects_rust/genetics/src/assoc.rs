use cmatrix::dense::{BaseCBits, CBitsRef};

use crate::sample::{BasePhe, Phe};

pub fn judge_sign_odds_ratio(gbits: &CBitsRef, phe: &Phe, case_n: usize) -> bool {
    let (d1, n1, d0, n0) = gbits.stat_contingency_table(phe.inner(), case_n);

    let odds_ratio =
        (d1 as f64 + 0.5) * (n0 as f64 + 0.5) / ((d0 as f64 + 0.5) * (n1 as f64 + 0.5));
    // println!(
    //     "odds_ratio: {}, d1: {}, n1: {}, d0: {}, n0: {}",
    //     odds_ratio, d1, n1, d0, n0
    // );

    odds_ratio >= 1.0
}

#[cfg(test)]
mod tests {
    use cmatrix::dense::CBits;

    use super::*;

    #[test]
    fn test_judge_sign_odds_ratio() {
        // contingency table:
        //           | high | low |
        // case |  1  |  1  |
        // control |  1   |  2  |
        // or=2.5/1.5>1

        let gbits = CBits::new(&vec![true, false, true, false, false]);
        let phe = Phe::new(&vec![true, true, false, false, false]);
        let case_n = phe.count();

        assert!(judge_sign_odds_ratio(&gbits.as_cbits_ref_b(), &phe, case_n));

        // contingency table:
        //           | high | low |
        // case |  1  |  1  |
        // control |  1   |  1  |
        // or=1

        let gbits = CBits::new(&vec![true, false, true, false]);
        let phe = Phe::new(&vec![true, true, false, false]);
        let case_n = phe.count();

        assert!(judge_sign_odds_ratio(&gbits.as_cbits_ref_b(), &phe, case_n));

        // contingency table:
        //           | high | low |
        // case |  0  |  1  |
        // control |  1   |  5  |
        // or=0.5*5.5/(1.5*1.5)>1
        let gbits = CBits::new(&vec![false, true, false, false, false, false, false]);
        let phe = Phe::new(&vec![true, false, false, false, false, false, false]);
        let case_n = phe.count();

        assert!(judge_sign_odds_ratio(&gbits.as_cbits_ref_b(), &phe, case_n));

        // contingency table:
        //           | high | low |
        // case |  0  |  1  |
        // control |  1   |  2  |
        // or=0.5*2.5/(1.5*1.5)<1

        let gbits = CBits::new(&vec![false, true, false, false]);
        let phe = Phe::new(&vec![true, false, false, false]);
        let case_n = phe.count();

        assert!(!judge_sign_odds_ratio(
            &gbits.as_cbits_ref_b(),
            &phe,
            case_n
        ));

        // Sign should be positive for too rare vars (no samples have rare var)
        // contingency table:
        //           | high | low |
        // case |  0  |  1  |
        // control |  0   |  5  |
        // or=0.5*5.5/(0.5*1.5)>1
        let gbits = CBits::new(&vec![false, false, false, false, false, false]);
        let phe = Phe::new(&vec![true, false, false, false, false, false]);
        let case_n = phe.count();

        assert!(judge_sign_odds_ratio(&gbits.as_cbits_ref_b(), &phe, case_n));
    }
}
