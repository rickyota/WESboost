//  When adding one parameter, modify the following 4 places
// 1. BoostParamCommon.set_batch_way()
// 2. BatchWay::from_str()
// 3. BoostParamCommonTrait.batch_way()
// 4. BoostParamCommon.batch_way()

//use core::panic;
use std::{collections::HashSet, str::FromStr};

use genetics::{LdCriteria, SnvId};

pub type AggParam = AggParamCommon;

#[derive(PartialEq, Clone, Debug, Default)]
pub struct AggParamCommon {
    loss_func: AggLossFunc,
    group_freq_thre: Option<f64>,
    skip_neg: bool,
}

pub trait AggParamCommonTrait {
    fn boost_param_common(&self) -> &AggParamCommon;

    //fn boost_param_common_mut(&mut self) -> &mut BoostParamCommon;

    fn loss_func(&self) -> AggLossFunc {
        self.boost_param_common().loss_func()
    }

    fn group_freq_thre(&self) -> Option<f64> {
        self.boost_param_common().group_freq_thre
    }
    fn skip_neg(&self) -> bool {
        self.boost_param_common().skip_neg
    }
}

impl AggParamCommonTrait for AggParamCommon {
    fn boost_param_common(&self) -> &AggParamCommon {
        self
    }
}

impl AggParamCommon {
    pub fn loss_func(&self) -> AggLossFunc {
        self.loss_func
    }

    // pub fn new_loss_func(loss_func: &str) -> Self {
    //     Self {
    //         loss_func: AggLossFunc::from_str(loss_func).unwrap(),
    //     }
    // }

    pub fn set_loss_func(self, loss_func: &str) -> Self {
        Self {
            loss_func: AggLossFunc::from_str(loss_func).unwrap(),
            ..self
        }
    }

    pub fn set_group_freq_thre(self, group_freq_thre: Option<f64>) -> Self {
        Self {
            group_freq_thre,
            ..self
        }
    }

    pub fn set_skip_neg(self, skip_neg: bool) -> Self {
        Self { skip_neg, ..self }
    }
}

#[derive(Eq, PartialEq, Copy, Clone, Hash, Debug)]
pub enum AggLossFunc {
    // Exp,
    ErrorRate,
    // logistic regression
    Logistic,
    // logitboost
    // Logit,
    Pval,
}

impl FromStr for AggLossFunc {
    type Err = String;
    fn from_str(str: &str) -> Result<Self, Self::Err> {
        match str {
            // "exp" => Ok(AggLossFunc::Exp),
            "logistic" => Ok(AggLossFunc::Logistic),
            "error" => Ok(AggLossFunc::ErrorRate),
            "pval" => Ok(AggLossFunc::Pval),
            _ => Err(format!("Unknown LossFunc: {}", str)),
        }
    }
}

impl Default for AggLossFunc {
    fn default() -> Self {
        Self::Logistic
        // Self::Logit
        // Self::Exp
    }
}

#[cfg(test)]
mod tests {
    //use super::*;

    //#[test]
    //fn test_boosting_params() {
    //    let lrs = vec![1.0, 0.1];
    //    let boost_params = BoostParamLrs::default().set_learning_rates(lrs);

    //    let mut boost_iter = boost_params.into_iter();

    //    assert_eq!(boost_iter.next().unwrap().learning_rate(), 1.0);
    //    assert_eq!(boost_iter.next().unwrap().learning_rate(), 0.1);
    //}
}
