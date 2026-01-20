// C Wrapper for C++
// https://akitsu-sanae.hatenablog.com/entry/2016/12/21/010321
// https://nachtimwald.com/2017/08/18/wrapping-c-objects-in-()

#pragma once
#include "pgenlibr.h"
#include "pvar.h"

// for int8_t
#include <stdint.h>
// error on mc
// #include <cstdint>

void foo();
int gcd(int a, int b);
// std::string ghi();
int ghi(const char *);

// int ghi_wrap(const char *str) { return 89; }

// class MyClass {
//   public:
//     MyClass();
//
//     int get_abc();
//     ~MyClass();
//
//   private:
//     int abc;
// };

// MyClass *myclass_create();
// void myclass_delete(MyClass *m);
// int myclass_get_abc(MyClass *m);

// int pgenreader_get_a();

// only make this since debugging in C is troublesome
// create pgenreader_load_whole() etc. in rust
// signed char is i32 in rust
int pgenreader_load_snvs_extract(int8_t *genot, const char *fgenot, const char *fsnv, int m_in, bool const *const use_snvs, const int *snvs_idx,
                                 const int *alleles_idx, int sample_n_in, int const *const sample_subset, int sample_subset_n, int nthr);
// int pgenreader_load_snvs_extract(int8_t *genot, const char *fgenot,const char *fsnv, int snv_start, int snv_end, bool const *const use_snvs, int
// sample_n_in,
//                                  int const *const sample_subset, int sample_subset_n, int nthr);
