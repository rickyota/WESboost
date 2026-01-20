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


## <a name="news"></a>News

- v0.1 (Jan 20, 2026)
    - Initial version.


## <a name="introduction"></a>Introduction

WESboost is a polygenic score method to capture rare variant effects. It is based on the GenoBoost framework, which uses gradient boosting machines to model complex genetic architectures. WESboost extends GenoBoost by incorporating functional annotations and optimizing for whole-exome sequencing (WES) data.


