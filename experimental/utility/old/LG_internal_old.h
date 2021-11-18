
#if 0
int LG_qsort_1b    // sort array A of size 2-by-n, using 1 key (A [0][])
(
    int64_t *LG_RESTRICT A_0,       // size n array
    LG_void *LG_RESTRICT A_1,       // size n array
    const size_t xsize,             // size of entries in A_1
    const int64_t n
) ;

void LG_qsort_1b_size1  // LG_qsort_1b with A1 with sizeof = 1
(
    int64_t *LG_RESTRICT A_0,       // size n array
    uint8_t *LG_RESTRICT A_1,       // size n array
    const int64_t n
) ;

void LG_qsort_1b_size2  // LG_qsort_1b with A1 with sizeof = 2
(
    int64_t *LG_RESTRICT A_0,       // size n array
    uint16_t *LG_RESTRICT A_1,      // size n array
    const int64_t n
) ;

void LG_qsort_1b_size4  // LG_qsort_1b with A1 with sizeof = 4
(
    int64_t *LG_RESTRICT A_0,       // size n array
    uint32_t *LG_RESTRICT A_1,      // size n array
    const int64_t n
) ;

void LG_qsort_1b_size8  // LG_qsort_1b with A_1 with sizeof = 8
(
    int64_t *LG_RESTRICT A_0,       // size n array
    uint64_t *LG_RESTRICT A_1,      // size n array
    const int64_t n
) ;

typedef struct
{
    uint8_t stuff [16] ;            // not accessed directly
}
LG_blob16 ;                         // sizeof (LG_blob16) is 16.

void LG_qsort_1b_size16 // LG_qsort_1b with A_1 with sizeof = 16
(
    int64_t *LG_RESTRICT A_0,       // size n array
    LG_blob16 *LG_RESTRICT A_1,     // size n array
    const int64_t n
) ;
#endif

