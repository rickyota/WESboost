
#include "pgenlibr_wrapc.hpp"

// #include <stdio.h>
#include <string>

#ifdef _OPENMP
#include <omp.h>
#endif

void foo() { printf("Hello, pgenlib!!\n"); }
int gcd(int a, int b) {
    int t;
    while(b != 0) {
        t = a % b;
        a = b;
        b = t;
    }
    return a;
}

int ghi(const char *str) { return 89; }

// genot double* : #snvs(true in use_snvs) * sample_subset_n
// genot is number of non-reference alleles
// If either of the allele is missing, then 3
int pgenreader_load_snvs_extract(int8_t *genot, const char *fgenot, const char *fsnv, int m_in, bool const *const use_snvs, const int *snvs_idx,
                                 const int *alleles_idx, int sample_n_in, int const *const sample_subset, int sample_subset_n, int nthr) {

    printf("fgenot: %s \n", fgenot);
    printf("fsnv: %s \n", fsnv);
    printf("cpp nthr: %d\n", nthr);

    RPgenReader *pread = new RPgenReader();
    std::string fgenot_str(fgenot);

    // might have better way
    std::vector<int> sample_subset_1based_v(sample_subset_n);
    std::copy(sample_subset, sample_subset + sample_subset_n, sample_subset_1based_v.begin());
    for(int &d : sample_subset_1based_v) {
        d += 1;
    }

    // size_t const snv_len = snv_end - snv_start;

    // TODO: use map
    // std::unordered_map<size_t, size_t> snv_index;
    size_t *snv_index = static_cast<size_t *>(malloc(sizeof(size_t) * m_in));
    size_t mi = 0;
    for(size_t m_in_i = 0; m_in_i < m_in; m_in_i++) {
        if(use_snvs[m_in_i]) {
            snv_index[m_in_i] = mi;
            mi++;
        } else {
            // otherwise not initialized
            snv_index[m_in_i] = SIZE_MAX;
        }
    }
    // mi = m_read

    // printf("use snvs %zu\n", snv_index_use);

    // load pvar
    RPvar *pvar = new RPvar();
    std::string fsnv_str(fsnv);
    pvar->Load(fsnv_str, false, false);

    printf("loaded pvar\n");
    fflush(stdout);

    printf("cpp nthr %d\n", nthr);
    fflush(stdout);
    pread->Load(fgenot_str, pvar, sample_n_in, sample_subset_1based_v, nthr);

    size_t m_in_i;
#ifdef _OPENMP
    //  disable dynamic
    //  https://stackoverflow.com/questions/11095309/openmp-set-num-threads-is-not-working
    // omp_set_dynamic(0);
    omp_set_num_threads(nthr);
#endif

#ifdef _OPENMP
    size_t const thrn = omp_get_max_threads();
#else
    size_t const thrn = 1;
#endif
    int **buf_thrs = static_cast<int **>(malloc(sizeof(int *) * thrn));
    // double **buf_thrs = static_cast<double **>(malloc(sizeof(double *) * thrn));
    for(size_t thri = 0; thri < thrn; thri++) {
        buf_thrs[thri] = static_cast<int *>(malloc(sizeof(int) * sample_subset_n));
        // buf_thrs[thri] = static_cast<double *>(malloc(sizeof(double) * sample_subset_n));
    }

#ifdef _OPENMP
#pragma omp parallel for private(mi)
#endif
    for(m_in_i = 0; m_in_i < m_in; m_in_i++) {
        //  for(mi = 0; mi < snv_end - snv_start; mi++) {

        if(!use_snvs[m_in_i]) {
            continue;
        }

#ifdef _OPENMP
        size_t const thri = omp_get_thread_num();
#else
        size_t const thri = 0;
#endif

        int *buf_thr = buf_thrs[thri];
        // double *buf_thr = buf_thrs[thri];

        int snv_in_i = snvs_idx[m_in_i];
        // 1-based
        int allele_in_i = alleles_idx[m_in_i] + 1;

        pread->ReadIntHardcalls(buf_thr, sample_subset_n, thri, snv_in_i, allele_in_i);

        //  should use ReadHardcalls?
        //  -> different from [manual of pgenlibr](https://cran.r-project.org/web/packages/pgenlibr/pgenlibr.pdf).
        //  1: # of first ALT
        // pread->ReadIntHardcalls(buf_thr, sample_subset_n, thri, m_in_i + snv_start, 1);
        // memory inefficient
        // pread->ReadHardcalls(buf_thr, sample_subset_n, thri, mi + snv_start, 1);
        //  pread->Read(buf_thr, sample_subset_n, thri, mi + snv_start, 1);
        //   2: # of second ALT
        // pread->ReadIntHardcalls(buf_thr, sample_subset_n, thri, mi + snv_start, 2);
        //   pread->Read(buf_thr, sample_subset_n, thri, mi + snv_start, 2);
        //   3: # of third ALT
        //   pread->Read(buf_thr, sample_subset_n, thri, mi + snv_start, 3);

        size_t const snv_index_use = snv_index[m_in_i];
        // printf("mi, use_index %zu, %zu\n", mi, use_index);

        for(size_t samplei = 0; samplei < sample_subset_n; samplei++) {
            genot[snv_index_use * sample_subset_n + samplei] = (int8_t)buf_thr[samplei];
        }

        // pread->Read(genot + mi * sample_subset_n, sample_subset_n, thread_num, mi + snv_start, 1);
    }

    return nthr;
}

// int pgenreader_load_snvs_extract(int8_t *genot, const char *filename, int snv_start, int snv_end, bool const *const use_snvs, int sample_n_in,
//                                  int const *const sample_subset, int sample_subset_n, int nthr) {}

// bi-allelic
// genot double* : #snvs(true in use_snvs) * sample_subset_n
// genot is number of non-reference alleles
// If either of the allele is missing, then 3
// int pgenreader_load_snvs_extract(int8_t *genot, const char *fgenot, const char *fsnv, int snv_start, int snv_end, bool const *const use_snvs,
//                                  int sample_n_in, int const *const sample_subset, int sample_subset_n, int nthr) {

//     printf("fgenot: %s \n", fgenot);
//     printf("fsnv: %s \n", fsnv);
//     printf("cpp nthr: %d\n", nthr);

//     RPgenReader *pread = new RPgenReader();
//     std::string fgenot_str(fgenot);

//     // might have better way
//     std::vector<int> sample_subset_1based_v(sample_subset_n);
//     std::copy(sample_subset, sample_subset + sample_subset_n, sample_subset_1based_v.begin());
//     for(int &d : sample_subset_1based_v) {
//         d += 1;
//     }

//     size_t const snv_len = snv_end - snv_start;

//     // TODO: use map
//     // std::unordered_map<size_t, size_t> snv_index;
//     size_t *snv_index = static_cast<size_t *>(malloc(sizeof(size_t) * snv_len));
//     size_t use_index = 0;
//     for(size_t mi = 0; mi < snv_len; mi++) {
//         if(use_snvs[mi]) {
//             snv_index[mi] = use_index;
//             use_index++;
//         } else {
//             // otherwise not initialized
//             snv_index[mi] = SIZE_MAX;
//         }
//     }

//     printf("use snvs %zu\n", use_index);

//     // load pvar
//     RPvar *pvar = new RPvar();
//     std::string fsnv_str(fsnv);
//     // pvar->Load(fsnv_str,true,true);
//     pvar->Load(fsnv_str, false, false);

//     printf("loaded pvar\n");
//     fflush(stdout);

//     printf("cpp nthr %d\n", nthr);
//     fflush(stdout);
//     pread->Load(fgenot_str, pvar, sample_n_in, sample_subset_1based_v, nthr);

//     size_t mi;
// #ifdef _OPENMP
//     //  disable dynamic
//     //  https://stackoverflow.com/questions/11095309/openmp-set-num-threads-is-not-working
//     // omp_set_dynamic(0);
//     omp_set_num_threads(nthr);
// #endif

// #ifdef _OPENMP
//     size_t const thrn = omp_get_max_threads();
// #else
//     size_t const thrn = 1;
// #endif
//     int **buf_thrs = static_cast<int **>(malloc(sizeof(int *) * thrn));
//     // double **buf_thrs = static_cast<double **>(malloc(sizeof(double *) * thrn));
//     for(size_t thri = 0; thri < thrn; thri++) {
//         buf_thrs[thri] = static_cast<int *>(malloc(sizeof(int) * sample_subset_n));
//         // buf_thrs[thri] = static_cast<double *>(malloc(sizeof(double) * sample_subset_n));
//     }

// #ifdef _OPENMP
// #pragma omp parallel for private(mi)
// #endif
//     for(mi = 0; mi < snv_len; mi++) {
//         // for(mi = 0; mi < snv_end - snv_start; mi++) {

//         if(!use_snvs[mi]) {
//             continue;
//         }

// #ifdef _OPENMP
//         size_t const thri = omp_get_thread_num();
// #else
//         size_t const thri = 0;
// #endif

//         int *buf_thr = buf_thrs[thri];
//         // double *buf_thr = buf_thrs[thri];
//         //  should use ReadHardcalls?
//         //  -> different from [manual of pgenlibr](https://cran.r-project.org/web/packages/pgenlibr/pgenlibr.pdf).
//         //  1: # of first ALT
//         pread->ReadIntHardcalls(buf_thr, sample_subset_n, thri, mi + snv_start, 1);
//         // memory inefficient
//         // pread->ReadHardcalls(buf_thr, sample_subset_n, thri, mi + snv_start, 1);
//         //  pread->Read(buf_thr, sample_subset_n, thri, mi + snv_start, 1);
//         //   2: # of second ALT
//         // pread->ReadIntHardcalls(buf_thr, sample_subset_n, thri, mi + snv_start, 2);
//         //   pread->Read(buf_thr, sample_subset_n, thri, mi + snv_start, 2);
//         //   3: # of third ALT
//         //   pread->Read(buf_thr, sample_subset_n, thri, mi + snv_start, 3);

//         size_t const use_index = snv_index[mi];
//         // printf("mi, use_index %zu, %zu\n", mi, use_index);

//         for(size_t samplei = 0; samplei < sample_subset_n; samplei++) {
//             genot[use_index * sample_subset_n + samplei] = (int8_t)buf_thr[samplei];
//         }

//         // pread->Read(genot + mi * sample_subset_n, sample_subset_n, thread_num, mi + snv_start, 1);
//     }

//     return nthr;
// }
