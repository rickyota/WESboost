# WESBoost 


## <a name="started"></a>Getting Started

Run `wesboost.sh` for a full example of training and scoring.


```bash
$ wesboost train \
    --dir ./result \
    --file-genot ./example/genot \
    --file-phe ./example/genot.cov \
    --cov age,sex \
    --major-a2-train \
    --verbose
```



## Table of Contents
- [WESBoost](#wesboost)
  - [Getting Started](#getting-started)
  - [Table of Contents](#table-of-contents)
  - [News](#news)
  - [Introduction](#introduction)
  - [Users' Guide](#users-guide)
    - [Installation](#installation)
      - [Plink1 Input](#plink1-input)
      - [Plink2 Input (compile)](#plink2-input-compile)
    - [Train WESBoost Model](#train-wesboost-model)
      - [Simplest Usage](#simplest-usage)
      - [Without Validation](#without-validation)
      - [Input Plink2](#input-plink2)
    - [Calculate Sample Scores](#calculate-sample-scores)
      - [Simplest Usage](#simplest-usage-1)
      - [Without Validation](#without-validation-1)
      - [Input Plink2](#input-plink2-1)


## <a name="news"></a>News

- v0.1 (Jan 20, 2026)
    - Initial version.


## <a name="introduction"></a>Introduction

WESBoost is a polygenic score method to capture rare variant effects. It is based on the GenoBoost framework, which uses gradient boosting machines to model complex genetic architectures. WESBoost extends GenoBoost by incorporating functional annotations and optimizing for whole-exome sequencing (WES) data.



## <a name="user-guide"></a>Users' Guide

For now, the input genotype format is allowed for plink1 or plink2 only.

### <a name="install"></a>Installation

#### <a name="install-plink1"></a>Plink1 Input

If you want to input plink1, download a compiled program for Linux (tested on Rocky Linux<=8.9), macOS (tested on <=14.3.1), and Windows (not tested) from [here][release]. This should take less than 1 minute.

#### <a name="install-plink2-compile"></a>Plink2 Input (compile)

If you want to input plink2 genotype file, you can compile program by yourself as below or [use docker or singularity](#advanced-installation). You can use plink1 format as well.

First, install `rust` as instructed [here][rust-install] if not installed. Then,

```bash
git clone https://github.com/rickyota/wesboost.git
cd wesboost
cargo build --manifest-path ./projects_rust/Cargo.toml --release --bin wesboost
cp ./projects_rust/target/release/wesboost ./wesboost
```

and you can use `wesboost` program. This should take less than 5 minutes.



### <a name="train"></a>Train WESBoost Model

#### <a name="train-simple"></a>Simplest Usage

You can run WESBoost at least with plink1 genotype files and, in most cases, a covariates file.

```bash
$ ./wesboost train \
    --dir ./result \
    --file-genot ./example/genot \
    --file-phe ./example/genot.cov \
    --cov age,sex \
    --major-a2-train \
    --seed 55
```

This test code should take less than 2 minutes.


#### <a name="train-train-only"></a>Without Validation

If you want to treat all samples as a training dataset, use `--train-only` option. WESBoost produces SNV weights each for learning rate. Use `--iter-snv` or `--iter` to control the maximum number of SNVs or iterations for training.

```bash
$ ./wesboost train \
    --dir ./result \
    --file-genot ./example/genot \
    --file-phe ./example/genot.cov \
    --cov age,sex \
    --major-a2-train \
    --train-only \
    --iter-snv 10000
```

#### <a name="train-plink2"></a>Input Plink2

If you use plink2 genotype file (`.pgen`, `.psam` and `.pvar` or `.pvar.zst`), use `--genot-format plink2` or `--genot-format plink2-vzs`.

If the phenotype is accompanied by covariates in the phenotype file, use `--phe` for the phenotype name. If phenotypes and covariates are in plink2 psam file, do not use `--file-phe`.

Control/case format should be `0/1` or `1/2`.

```bash
$ ./genoboost train \
    --dir ./result \
    --file-genot ./example/genot2 \
    --genot-format plink2-vzs \
    --file-phe ./example/genot2.phe \
    --phe PHENO1 \
    --cov age,sex \
    --major-a2-train \
    --seed 55
```



### <a name="score"></a>Calculate Sample Scores

WESBoost returns a polygenic score for each sample. WESBoost outputs scores without covariates (`score.tsv`) and with covariates (`score.withcov.tsv`).



#### <a name="score-simple"></a>Simplest Usage

With the minimum options, WESBoost will calculate sample scores from SNV weights with the best parameters determined in the validation dataset.

```bash
$ ./wesboost score \
    --dir-score ./result_score \
    --dir-wgt ./result \
    --file-genot ./example/genot \
    --file-phe ./example/genot.cov \
    --cov age,sex
```


#### <a name="score-train-only"></a>Without Validation

If you did not use the validation dataset in the training phase, WESBoost will output sample scores for all parameters. You have to specify the number of SNVs in `--iters`.

```bash
$ ./wesboost score \
    --dir-score ./result_score \
    --dir-wgt ./result \
    --file-genot ./example/genot \
    --file-phe ./example/genot.cov \
    --cov age,sex \
    --train-only \
    --iters "10 20 50"
```


#### <a name="score-plink2"></a>Input Plink2

Use `--genot-format`, `--file-phe` etc. for plink2 as shown in [training phase](#train-plink2).

```bash
$ ./wesboost score \
    --dir ./result \
    --file-genot ./example/genot2 \
    --genot-format plink2-vzs \
    --file-phe ./example/genot2.phe \
    --cov age,sex
```
