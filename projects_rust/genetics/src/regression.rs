use rstat::{normal::UvNormal, Distribution};
use statrs::distribution::{ChiSquared, ContinuousCDF, Normal};
use std::f64::consts::{LN_10, LN_2};

use crate::samples::prelude::*;
use crate::{wgt::Coef, Wgt};
use crate::{CovId, Covs, CovsTrait, Samples};

//mod logistic;
//mod smartcore;
mod linfa;
#[cfg(feature = "pyo3")]
mod sklearn_py;

use crate::vec;

// linfa
// [ref](https://github.com/zupzup/rust-ml-example/blob/main/src/main.rs)
// [linfa_logistic ref](https://crates.io/crates/linfa-logistic)
// https://docs.rs/linfa-logistic/0.7.0/linfa_logistic/type.LogisticRegression.html
pub fn logistic_regression_covs(samples: &Samples) -> Vec<Wgt> {
    if let Some(covs) = samples.covs() {
        let (coefs, intercept) = linfa::logistic_regression_covs_samples(samples);
        //let (coefs, intercept) = smartcore::logistic_regression_covs_samples(samples);

        assert_eq!(coefs.len(), covs.covs_n());

        wgts_cov_from_coef(&coefs, intercept, covs)

        //let mut wgts_cov: Vec<Wgt> = Vec::new();

        //let cov_id = CovId::new_const();
        //let wgt_const = Wgt::new_cov(cov_id, Coef::Linear(intercept));
        //wgts_cov.push(wgt_const);

        ////let coefs = lr.params().iter().map(|x| *x).collect::<Vec<f64>>();
        //////let coefs = lr.coefficients().clone().to_row_vector();
        ////assert_eq!(coefs.len(), covs_val.covs_n());
        //let cov_indexs = covs.cov_indexs().unwrap();
        //for wgt_i in 0..covs.covs_n() {
        //    let coef = coefs[wgt_i];
        //    if coef.is_nan() {
        //        panic!("coef is NaN.");
        //    }
        //    let cov_id = CovId::new_cov(cov_indexs[wgt_i].name().to_owned());
        //    let wgt_cov = Wgt::new_cov(cov_id, Coef::Linear(coef));
        //    wgts_cov.push(wgt_cov);
        //}

        //wgts_cov
    } else {
        // regression on const

        log::debug!("Regression on const value only.");

        log::debug!("ys true count {}", samples.phe_unwrap().count());
        log::debug!("ys false count {}", samples.phe_unwrap().count_false());

        let intercept = logistic_regression_const_samples(samples);
        wgts_cov_from_intercept(intercept)

        //let y1_n = samples.phe_unwrap().count();
        //let y0_n = samples.phe_unwrap().count_false();

        //let intercept = logistic_regression_const(y0_n, y1_n);

        //let mut wgts_cov: Vec<Wgt> = Vec::new();
        //if intercept.is_nan() {
        //    panic!("Intercept is NaN.");
        //}
        //let cov_id = CovId::new_const();
        //let wgt_const = Wgt::new_cov(cov_id, Coef::Linear(intercept));
        //wgts_cov.push(wgt_const);

        //wgts_cov
    }
}

fn logistic_regression_const_samples(samples: &Samples) -> f64 {
    let y1_n = samples.phe_unwrap().count();
    let y0_n = samples.phe_unwrap().count_false();

    let intercept = logistic_regression_const(y0_n, y1_n);

    if intercept.is_nan() {
        panic!("Intercept is NaN.");
    }
    intercept
}

fn logistic_regression_const(y0_n: usize, y1_n: usize) -> f64 {
    // regression on const
    let y1_n = y1_n as f64;
    let y0_n = y0_n as f64;

    // TODO: test
    let intercept = (y1_n / y0_n).ln();
    intercept
}

fn wgts_cov_from_coef(coefs: &[f64], intercept: f64, covs: &Covs) -> Vec<Wgt> {
    let mut wgts_cov: Vec<Wgt> = Vec::new();

    let cov_id = CovId::new_const();
    let wgt_const = Wgt::new_cov(cov_id, Coef::Linear(intercept));
    wgts_cov.push(wgt_const);

    let cov_indexs = covs.cov_indexs().unwrap();
    for wgt_i in 0..covs.covs_n() {
        let coef = coefs[wgt_i];
        // already checked
        //if coef.is_nan() {
        //    panic!("coef is NaN.");
        //}
        let cov_id = CovId::new_cov(cov_indexs[wgt_i].name().to_owned());
        let wgt_cov = Wgt::new_cov(cov_id, Coef::Linear(coef));
        wgts_cov.push(wgt_cov);
    }

    wgts_cov
}

fn wgts_cov_from_intercept(intercept: f64) -> Vec<Wgt> {
    let cov_id = CovId::new_const();
    let wgt_const = Wgt::new_cov(cov_id, Coef::Linear(intercept));
    vec![wgt_const]
}

#[cfg(not(feature = "pyo3"))]
pub fn logistic_regression_covs_sampleweights(_: &Samples, _: &[f64]) -> Vec<Wgt> {
    panic!("Use feature=pyo3");
}

#[cfg(feature = "pyo3")]
pub fn logistic_regression_covs_sampleweights(samples: &Samples, ps: &[f64]) -> Vec<Wgt> {
    if let Some(covs) = samples.covs() {
        let (coefs, intercept) =
            sklearn_py::logistic_regression_covs_sample_sampleweights(samples, ps);

        assert_eq!(coefs.len(), covs.covs_n());

        wgts_cov_from_coef(&coefs, intercept, covs)
    } else {
        // regression on const

        unimplemented!("Not implemented yet.")

        //log::debug!("Regression on const value only.");

        //log::debug!("ys true count {}", samples.phe_unwrap().count());
        //log::debug!("ys false count {}", samples.phe_unwrap().count_false());

        //let y1_n = samples.phe_unwrap().count();
        //let y0_n = samples.phe_unwrap().count_false();

        //let intercept = logistic_regression_const(y0_n, y1_n);

        //let mut wgts_cov: Vec<Wgt> = Vec::new();
        //if intercept.is_nan() {
        //    panic!("Intercept is NaN.");
        //}
        //let cov_id = CovId::new_const();
        //let wgt_const = Wgt::new_cov(cov_id, Coef::Linear(intercept));
        //wgts_cov.push(wgt_const);

        //wgts_cov
    }
}

// for nagelkerke
pub fn logistic_regression_vec2(
    v: Vec<Vec<f64>>, // col major; unnormed 2d
    ys: Vec<bool>,
) -> (Vec<f64>, f64) {
    let means = vec::mean_v2(&v);
    let stds = vec::std_v2(&v);

    // convert to row major
    let v_norm = vec::convert_vec2d_to_row_major(&vec::norm_vec2d(v));
    let (v_norm_1d, (row_n, col_n)) = vec::convert_vec2d_to_vec1(v_norm);
    //let (v_norm_1d, (col_n, row_n)) = vec::convert_vec2d_to_vec1(v_norm);

    assert_eq!(col_n, 1);

    log::debug!("means {:?}", means);
    log::debug!("stds {:?}", stds);

    log::debug!("v_norm_1d[0],[1] {:?}, {:?}", v_norm_1d[0], v_norm_1d[1]);

    linfa::logistic_regression_vec_linfa_normalize(v_norm_1d, col_n, row_n, ys, &means, &stds)
}

// for nagelkerke with small score variance
// 1d only for now
pub fn logistic_regression_vec_no_intercept_1col(
    //v: Vec<Vec<f64>>, // col major; unnormed 2d
    v: Vec<f64>, // col major; unnormed 2d
    ys: Vec<bool>,
) -> f64 {
    //let means = vec::mean_v2(&v);
    //let stds = vec::std_v2(&v);

    // convert to row major
    //let v_norm = vec::convert_vec2d_to_row_major(&vec::norm_vec2d(v));
    //let (v_norm_1d, (row_n, col_n)) = vec::convert_vec2d_to_vec1(v_norm);
    //let (v_norm_1d, (col_n, row_n)) = vec::convert_vec2d_to_vec1(v_norm);

    //assert_eq!(col_n, 1);

    //log::debug!("means {:?}", means);
    //log::debug!("stds {:?}", stds);

    //log::debug!("v_norm_1d[0],[1] {:?}, {:?}", v_norm_1d[0], v_norm_1d[1]);

    linfa::logistic_regression_vec_linfa_no_intercept_1col(v, ys)
}

/// error: digit loss
/// -> same as p_value_From_z
fn mlog10_pval_from_z(z: f64) -> f64 {
    let dist = UvNormal::new(0.0, 1.0).unwrap();
    let x = z.abs();
    let log_sf = dist.log_ccdf(&x);
    // 2* for two-tailed test
    let ln_p = LN_2 + log_sf;
    -ln_p / LN_10
}

/// error: digit loss
/// Computes the two-sided p-value from a z-statistic using the standard normal CDF.
fn p_value_from_z(z: f64) -> f64 {
    let normal = Normal::new(0.0, 1.0).unwrap();
    2.0 * (1.0 - normal.cdf(z.abs()))
}

/// Given a 2x2 contingency table with counts:
///
/// |              | Outcome = 1 | Outcome = 0 |
/// |--------------|-------------|-------------|
/// | Exposed      | a           | b           |
/// | Not Exposed  | c           | d           |
///
/// This function computes the logistic regression p-value (using a Wald test)
/// for testing whether the binary predictor is associated with the outcome.
/// It does so by calculating the log odds ratio and its standard error, then
/// computing the two-sided p-value.
///
/// TODO: If any cell is zero, a continuity correction (adding 0.5) is applied.
fn wald_test_mlog10_pval_zscore(a: u64, b: u64, c: u64, d: u64) -> Option<(f64, f64)> {
    // Convert counts to f64 and apply continuity correction if any cell is zero.
    let (a, b, c, d) = if a == 0 || b == 0 || c == 0 || d == 0 {
        return None;
        // (
        //     a as f64 + 0.5,
        //     b as f64 + 0.5,
        //     c as f64 + 0.5,
        //     d as f64 + 0.5,
        // )
    } else {
        (a as f64, b as f64, c as f64, d as f64)
    };

    // Calculate the odds ratio: OR = (a*d)/(b*c)
    let or = (a * d) / (b * c);
    // Compute the log odds ratio (this is the logistic regression coefficient)
    let log_or = or.ln();

    // Standard error of the log odds ratio:
    let se = (1.0 / a + 1.0 / b + 1.0 / c + 1.0 / d).sqrt();

    // Wald z-statistic:
    let z = log_or / se;

    // Two-sided p-value from z:
    Some((mlog10_pval_from_z(z), z))
}

/// Given a 2x2 contingency table with counts:
///
/// |              | Outcome = 1 | Outcome = 0 |
/// |--------------|-------------|-------------|
/// | Exposed      | a           | b           |
/// | Not Exposed  | c           | d           |
///
/// This function computes the logistic regression p-value (using a Wald test)
/// for testing whether the binary predictor is associated with the outcome.
/// It does so by calculating the log odds ratio and its standard error, then
/// computing the two-sided p-value.
///
/// TODO: If any cell is zero, a continuity correction (adding 0.5) is applied.
// fn wald_test_p_value(a: u64, b: u64, c: u64, d: u64) -> Option<f64> {
//     // Convert counts to f64 and apply continuity correction if any cell is zero.
//     let (a, b, c, d) = if a == 0 || b == 0 || c == 0 || d == 0 {
//         return None;
//         // (
//         //     a as f64 + 0.5,
//         //     b as f64 + 0.5,
//         //     c as f64 + 0.5,
//         //     d as f64 + 0.5,
//         // )
//     } else {
//         (a as f64, b as f64, c as f64, d as f64)
//     };

//     // Calculate the odds ratio: OR = (a*d)/(b*c)
//     let or = (a * d) / (b * c);
//     // Compute the log odds ratio (this is the logistic regression coefficient)
//     let log_or = or.ln();

//     // Standard error of the log odds ratio:
//     let se = (1.0 / a + 1.0 / b + 1.0 / c + 1.0 / d).sqrt();

//     // Wald z-statistic:
//     let z = log_or / se;

//     // Two-sided p-value from z:
//     Some(p_value_from_z(z))
// }

pub fn logistic_regression_cont_table_2_2_mlog10_pval_zscore(
    cont_table: &[Vec<f64>],
) -> Option<(f64, f64)> {
    wald_test_mlog10_pval_zscore(
        cont_table[0][0] as u64,
        cont_table[0][1] as u64,
        cont_table[1][0] as u64,
        cont_table[1][1] as u64,
    )
}

// error: digit loss
// pub fn logistic_regression_cont_table_2_2_p_value(cont_table: &[Vec<f64>]) -> Option<f64> {
//     wald_test_p_value(
//         cont_table[0][0] as u64,
//         cont_table[0][1] as u64,
//         cont_table[1][0] as u64,
//         cont_table[1][1] as u64,
//     )
// }

// TOFIX: digit loss
pub fn chi_square_cont_table_p_value(cont_table: &[Vec<f64>]) -> f64 {
    let rows = cont_table.len();
    let cols = cont_table[0].len();

    // Compute row sums, column sums, and total sum.
    let mut row_sums = vec![0.0; rows];
    let mut col_sums = vec![0.0; cols];
    let mut total = 0.0;
    for i in 0..rows {
        for j in 0..cols {
            row_sums[i] += cont_table[i][j];
            col_sums[j] += cont_table[i][j];
            total += cont_table[i][j];
        }
    }

    // Calculate the chi-square statistic.
    let mut chi_square_stat = 0.0;
    for i in 0..rows {
        for j in 0..cols {
            let expected = row_sums[i] * col_sums[j] / total;
            // Avoid division by zero if expected is 0.
            if expected > 0.0 {
                chi_square_stat += (cont_table[i][j] - expected).powi(2) / expected;
            }
        }
    }

    // Degrees of freedom = (number of rows - 1) * (number of columns - 1)
    let df = (rows - 1) * (cols - 1);
    let chi2_dist = ChiSquared::new(df as f64).unwrap();

    // p-value = 1 - CDF(chi_square_stat)
    1.0 - chi2_dist.cdf(chi_square_stat)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_logisticregression_const() {
        let y1_n = 100;
        let y0_n = 10;
        let intercept = logistic_regression_const(y0_n, y1_n);
        assert_float_absolute_eq!(intercept, (100.0f64 / 10.0).ln());
    }

    #[test]
    fn test_logistic_regression_vec2() {
        // col major
        let v = vec![vec![-1.0f64, -0.01, 0.01, 1.0]];
        let y = vec![false, true, false, true];

        //let v_norm = vec::convert_vec2d_to_row_major(&vec::norm_vec2d(v));
        //let (v_norm_1d, row_n, col_n) = vec::convert_vec2d_to_vec1_row_major(v_norm);
        let (coefs, intercept) = logistic_regression_vec2(v, y);
        //logistic_regression_covs_vec_linfa(v_norm_1d, row_n, col_n, y, &means, &stds);
        assert_eq!(coefs.len(), 1);
        assert_float_absolute_eq!(intercept, 0.0);

        // expected answer was calculated using https://stats.blue/Stats_Suite/logistic_regression_calculator.html
        assert_float_absolute_eq!(coefs[0], 5.2672, 1e-3);
    }

    #[test]
    fn test_mlog10_pval_from_z() {
        let z = 1.96; // Example z-statistic
        let mlog10_pval = mlog10_pval_from_z(z);
        // Expected value calculated using a statistical calculator
        // Two-tailed test: https://www.socscistatistics.com/pvalues/normaldistribution.aspx
        assert_float_absolute_eq!(mlog10_pval, 1.3010, 1e-4);

        let mlog10_pval2 = p_value_from_z(z).log10() * -1.0;
        assert_float_absolute_eq!(mlog10_pval, mlog10_pval2, 1e-4);
    }

    // TOFIX: remove [should_panic] when fixed
    #[test]
    #[should_panic]
    fn test_mlog10_pval_from_z_2() {
        let z = 10.0; // Example z-statistic
        let mlog10_pval = mlog10_pval_from_z(z);
        // Expected value calculated using a statistical calculator
        // Two-tailed test: https://www.socscistatistics.com/pvalues/normaldistribution.aspx
        assert_float_absolute_eq!(mlog10_pval, 1.3010, 1e-4);
    }

    #[test]
    fn test_logistic_regression_cont_table_2_2_p_value() {
        // TODO
    }

    #[test]
    fn test_chi_square_p_value() {
        let observed = vec![vec![10.0, 20.0], vec![30.0, 40.0]];
        let p_value = chi_square_cont_table_p_value(&observed);
        // https://www.socscistatistics.com/tests/chisquare/default2.aspx
        assert_float_absolute_eq!(p_value, 0.372998, 1e-6);
    }
}
