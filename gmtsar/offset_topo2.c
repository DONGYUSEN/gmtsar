#include "gmtsar.h"
#include <omp.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char **argv) {
    int i, j, k, i1, j1, k1;
    int is, js, ns;
    int ni, nj, ntot;
    int xshft, yshft, ib = 200;
    int imax = 0, jmax = 0;
    double ra, rt, avea;
    double suma, sumt, sumc, corr, denom;
    double maxcorr = -1e30;

    // 检查OpenMP支持
#ifdef _OPENMP
    printf("=== OpenMP 已启用 ===\n");
    // 设置线程数（使用所有可用核心）
	printf("最大可用线程数: %d\n", omp_get_max_threads());
    int nthreads = omp_get_max_threads()*5/6;
	if (nthreads < 2) nthreads = 2; 
    omp_set_num_threads(nthreads);
    printf("仅使用 %d 个线程进行计算\n", nthreads);
#else
    printf("=== 警告: OpenMP 未启用，程序将串行运行 ===\n");
    printf("编译时请添加 -fopenmp 选项启用并行计算\n");
#endif


    void *API = NULL;
    struct GMT_GRID *A = NULL, *T = NULL, *TS = NULL;

    if (argc < 6) {
        fprintf(stderr,
                "Usage: offset_topo2 amp_master.grd topo_ra.grd rshift ashift ns [topo_shift.grd]\n");
        exit(EXIT_FAILURE);
    }

    API = GMT_Create_Session(argv[0], 0U, 0U, NULL);

    xshft = atoi(argv[3]);
    yshft = atoi(argv[4]);
    ns = atoi(argv[5]);

    A = GMT_Read_Data(API, GMT_IS_GRID, GMT_IS_FILE, GMT_IS_SURFACE,
                      GMT_GRID_HEADER_ONLY, NULL, argv[1], NULL);
    T = GMT_Read_Data(API, GMT_IS_GRID, GMT_IS_FILE, GMT_IS_SURFACE,
                      GMT_GRID_HEADER_ONLY, NULL, argv[2], NULL);

    if (A->header->n_columns != T->header->n_columns) {
        fprintf(stderr, "Grid width mismatch\n");
        exit(EXIT_FAILURE);
    }

    if (argc >= 7) {
        TS = GMT_Create_Data(API, GMT_IS_GRID, GMT_IS_SURFACE, GMT_GRID_ALL,
                             NULL, A->header->wesn, A->header->inc,
                             A->header->registration, GMT_NOTSET, NULL);
    }

    GMT_Read_Data(API, GMT_IS_GRID, GMT_IS_FILE, GMT_IS_SURFACE,
                  GMT_GRID_DATA_ONLY, NULL, argv[1], A);
    GMT_Read_Data(API, GMT_IS_GRID, GMT_IS_FILE, GMT_IS_SURFACE,
                  GMT_GRID_DATA_ONLY, NULL, argv[2], T);

    ni = (A->header->n_rows < T->header->n_rows) ?
         A->header->n_rows : T->header->n_rows;
    nj = T->header->n_columns;

    /* ---------- Mean calculation (OpenMP) ---------- */
    suma = sumt = 0.0;
    ntot = 0;

    #pragma omp parallel for private(k) reduction(+:suma,sumt,ntot)
    for (i = 0; i < ni; i++) {
        for (j = 0; j < nj; j++) {
            k = i * nj + j;
            suma += A->data[k];
            sumt += T->data[k];
            ntot++;
        }
    }
    avea = suma / ntot;

    /* ---------- Cross-correlation search (OpenMP) ---------- */
    #pragma omp parallel
    {
        double maxcorr_t = -1e30;
        int imax_t = 0, jmax_t = 0;

        #pragma omp for collapse(2) schedule(dynamic)
        for (is = -ns + yshft; is <= ns + yshft; is++) {
            for (js = -ns + xshft; js <= ns + xshft; js++) {

                double sumc_l = 0.0, suma_l = 0.0, sumt_l = 0.0;

                for (i = ib; i < ni - ib; i++) {
                    i1 = i - is;
                    if (i1 < 1 || i1 >= ni - 1) continue;

                    for (j = ib; j < nj - ib; j++) {
                        j1 = j - js;
                        if (j1 < 1 || j1 >= nj - 1) continue;

                        k = i * nj + j;
                        k1 = i1 * nj + j1;

                        ra = A->data[k] - avea;
                        rt = T->data[k1 + 1] - T->data[k1 - 1];

                        sumc_l += ra * rt;
                        suma_l += ra * ra;
                        sumt_l += rt * rt;
                    }
                }

                denom = suma_l * sumt_l;
                if (denom > 0.0) {
                    corr = sumc_l / sqrt(denom);
                    if (corr > maxcorr_t) {
                        maxcorr_t = corr;
                        imax_t = is;
                        jmax_t = js;
                    }
                }
            }
        }

        #pragma omp critical
        {
            if (maxcorr_t > maxcorr) {
                maxcorr = maxcorr_t;
                imax = imax_t;
                jmax = jmax_t;
            }
        }
    }

    printf("optimal: rshift=%d ashift=%d maxcorr=%g\n",
           jmax, imax, maxcorr);

    /* ---------- Optional output ---------- */
    if (argc >= 7) {
        #pragma omp parallel for private(i1,j1,k,k1)
        for (i = 0; i < ni; i++) {
            i1 = i - imax;
            for (j = 0; j < nj; j++) {
                j1 = j - jmax;
                k = i * nj + j;
                k1 = i1 * nj + j1;
                TS->data[k] =
                    (i1 >= 0 && i1 < ni && j1 >= 0 && j1 < nj) ? T->data[k1] : 0.0f;
            }
        }

        GMT_Write_Data(API, GMT_IS_GRID, GMT_IS_FILE, GMT_IS_SURFACE,
                       GMT_GRID_ALL, NULL, argv[6], TS);
    }

    GMT_Destroy_Session(API);
    return EXIT_SUCCESS;
}
