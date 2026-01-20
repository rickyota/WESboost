
/*
 *
 * File obtained from pgenlibr R library:
 * https://github.com/chrchang/plink-ng/tree/master/2.0/pgenlibr
 * Version:
 * https://github.com/chrchang/plink-ng/releases/tag/v2.0.0-a.6.8
 *
 * License info obtained from DESCRIPTION file:
 * https://github.com/chrchang/plink-ng/blob/master/2.0/pgenlibr/DESCRIPTION
 * -----------------------------------------------------
    Package: pgenlibr
    Type: Package
    Title: PLINK 2 Binary (.pgen) Reader
    Version: 0.4.0
    Date: 2025-01-15
    Author: Christopher Chang
    Maintainer: Christopher Chang <chrchang@alumni.caltech.edu>
    Description: A thin wrapper over PLINK 2's core libraries which provides an R
    interface for reading .pgen files.  A minimal .pvar loader is also included.
    License: LGPL (>= 3)
    Imports: Rcpp (>= 1.0.1)
    LinkingTo: Rcpp
 * -----------------------------------------------------

 *  Modified by ricky
 *
 * This file remains under LGPL v3 license (license is in same directory as this file)
 */

#pragma once
#include "include/pgenlib_ffi_support.h"
#include "include/pgenlib_read.h"
#include "include/pvar_ffi_support.h"
#include "pvar.h" // includes Rcpp
#include <memory>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector>

int get_int_abc();

class RPgenReader {
  public:
    // imitates Python/pgenlib.pyx
    RPgenReader();

    int abc();
    int get_a();
    std::string def();

#if __cplusplus >= 201103L
    RPgenReader(const RPgenReader &) = delete;
    RPgenReader &operator=(const RPgenReader &) = delete;
#endif

    // pvar is necessary for multi-allelic
    void Load(std::string filename, RPvar *pvar, uint32_t cur_sample_ct, std::vector<int> sample_subset_1based, int nthr);
    // void Load(String filename, Nullable<List> pvar, Nullable<int> raw_sample_ct, Nullable<IntegerVector> sample_subset_1based);

    // uint32_t GetRawSampleCt() const;

    // uint32_t GetSubsetSize() const;

    // uint32_t GetVariantCt() const;

    // uint32_t GetAlleleCt(uint32_t variant_idx) const;

    // uint32_t GetMaxAlleleCt() const;

    // bool HardcallPhasePresent() const;

    void ReadIntHardcalls(int *buf, size_t const &n, int const &thr, int variant_idx, int allele_idx);
    // void ReadIntHardcalls(IntegerVector buf, int variant_idx, int allele_idx);

    void ReadHardcalls(double *buf, size_t const &n, int const &thr, int variant_idx, int allele_idx);
    // void ReadHardcalls(double *buf, size_t const &n, int const &thr, int variant_idx, int allele_idx);
    //  void ReadHardcalls(NumericVector buf, int variant_idx, int allele_idx);

    void Read(double *buf, size_t const &n, int const &thr, int variant_idx, int allele_idx);
    // void Read(NumericVector buf, int variant_idx, int allele_idx);

    // void ReadAlleles(IntegerMatrix acbuf, Nullable<LogicalVector> phasepresent_buf, int variant_idx);

    // void ReadAllelesNumeric(NumericMatrix acbuf, Nullable<LogicalVector> phasepresent_buf, int variant_idx);

    // void ReadIntList(IntegerMatrix buf, IntegerVector variant_subset);

    // void ReadList(NumericMatrix buf, IntegerVector variant_subset, bool meanimpute);

    // void FillVariantScores(NumericVector result, NumericVector weights, Nullable<IntegerVector> variant_subset);

    void Close();

    ~RPgenReader();

  private:
    plink2::PgenFileInfo *_info_ptr;
    plink2::RefcountedWptr *_allele_idx_offsetsp;
    plink2::RefcountedWptr *_nonref_flagsp;

    // have all below be threads specific
    std::vector<plink2::PgenReader *> _state_ptr;
    std::vector<uintptr_t *> _subset_include_vec;
    std::vector<uintptr_t *> _subset_include_interleaved_vec;
    std::vector<uint32_t *> _subset_cumulative_popcounts;
    std::vector<plink2::PgrSampleSubsetIndex> _subset_index;
    std::vector<uint32_t> _subset_size;
    // plink2::PgenReader *_state_ptr;
    // uintptr_t *_subset_include_vec;
    // uintptr_t *_subset_include_interleaved_vec;
    // uint32_t *_subset_cumulative_popcounts;
    // plink2::PgrSampleSubsetIndex _subset_index;
    // uint32_t _subset_size;

    std::vector<std::shared_ptr<plink2::PgenVariant>> _pgv;
    // std::vector<plink2::PgenVariant> _pgv;
    //  plink2::PgenVariant _pgv;

    std::vector<plink2::VecW *> _transpose_batch_buf;
    // plink2::VecW *_transpose_batch_buf;

    // kPglNypTransposeBatch (= 256) variants at a time, and then transpose
    std::vector<uintptr_t *> _multivar_vmaj_geno_buf;
    std::vector<uintptr_t *> _multivar_vmaj_phasepresent_buf;
    std::vector<uintptr_t *> _multivar_vmaj_phaseinfo_buf;
    std::vector<uintptr_t *> _multivar_smaj_geno_batch_buf;
    std::vector<uintptr_t *> _multivar_smaj_phaseinfo_batch_buf;
    std::vector<uintptr_t *> _multivar_smaj_phasepresent_batch_buf;

    void SetSampleSubsetInternal(std::vector<int> &sample_subset_1based, int const &thr);
    // void SetSampleSubsetInternal(IntegerVector sample_subset_1based);

    // void ReadAllelesPhasedInternal(int variant_idx);
};
