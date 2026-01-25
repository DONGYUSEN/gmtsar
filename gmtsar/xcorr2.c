/*	$Id: xcorr.c 73 2013-04-19 17:59:45Z pwessel $	*/
/***************************************************************************/
/* xcorr does a 2-D cross correlation on complex or real images            */
/* either using a time convolution or wavenumber multiplication.           */
/*                                                                         */
/***************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#ifdef _OPENMP
#include <omp.h>
#else
#define omp_get_thread_num() 0
#define omp_get_max_threads() 1
#define omp_get_num_threads() 1
#define omp_set_num_threads(n)
#endif

#include "gmtsar.h"
#include "xcorr.h" /* 确保 struct xcorr/locs/FCOMPLEX 定义可见 */

extern int verbose;
extern int debug;

extern void die(const char *msg, const char *arg); /* gmtsar 中通常提供 */
void fft_interpolate_1d(void *API, struct FCOMPLEX *c, int nx, struct FCOMPLEX *work, int ri);
float Cabs(struct FCOMPLEX z);
void read_xcorr_data(struct xcorr *xc, int iloc);
void print_complex(struct FCOMPLEX *c, int npy, int npx, int flag);
void do_time_corr(struct xcorr *xc, int iloc);
void do_freq_corr(void *API, struct xcorr *xc, int iloc);
void do_highres_corr(void *API, struct xcorr *xc, int iloc);
void print_results(struct xcorr *xc, int iloc);
void set_defaults(struct xcorr *xc);
void parse_command_line(int argc, char **argv, struct xcorr *xc, int *nfiles, int *input_flag, char *USAGE);
void handle_prm(void *API, char **argv, struct xcorr *xc, int nfiles);
void print_params(struct xcorr *xc);
void get_locations(struct xcorr *xc);

/* Usage 保留 */
char *USAGE = "xcorr2 [GMTSAR] - Compute 2-D cross-correlation of two images\n\n"
              /* ... (略) ... */
              ;

/* range 插值包装 */
static inline int do_range_interpolate(void *API, struct FCOMPLEX *c, int nx, int ri, struct FCOMPLEX *work) {
    fft_interpolate_1d(API, c, nx, work, ri);
    for (int i = 0; i < nx; ++i) {
        c[i].r = work[i + nx / 2].r;
        c[i].i = work[i + nx / 2].i;
    }
    return EXIT_SUCCESS;
}

/* 将值拷贝到线程私有缓冲并处理（线程私有版 assign_values） */
void assign_values_thread(void *API, struct xcorr *xc,
                          struct FCOMPLEX *d1_t, struct FCOMPLEX *d2_t,
                          struct FCOMPLEX *c1_t, struct FCOMPLEX *c2_t, struct FCOMPLEX *c3_t,
                          int *i1_t, int *i2_t, short *mask_t, struct FCOMPLEX *ritmp_t,
                          int iloc)
{
    if (!xc || !xc->loc) die("assign_values_thread: null xc or loc", NULL);
    if (iloc < 0 || iloc >= xc->nloc) die("assign_values_thread: iloc out of range", NULL);

    const int npx = xc->npx;
    const int npy = xc->npy;
    const int m_nx = xc->m_nx;
    const int s_nx = xc->s_nx;

    const int mx = xc->loc[iloc].x - npx / 2;
    const int sx = xc->loc[iloc].x + xc->x_offset - npx / 2;

    /* 填充 c1,c2,c3 */
    for (int iy = 0; iy < npy; iy++) {
        int d1_row_base = iy * m_nx + mx;
        int d2_row_base = iy * s_nx + sx;
        int c_row_base = iy * npx;
        for (int j = 0; j < npx; j++) {
            int k = c_row_base + j;
            c3_t[k].r = c3_t[k].i = 0.0f;
            c1_t[k].r = d1_t[d1_row_base + j].r;
            c1_t[k].i = d1_t[d1_row_base + j].i;
            c2_t[k].r = d2_t[d2_row_base + j].r;
            c2_t[k].i = d2_t[d2_row_base + j].i;
        }
    }

    if (xc->ri > 1) {
        int ri = xc->ri;
        for (int iy = 0; iy < npy; iy++) {
            do_range_interpolate(API, &c1_t[iy * npx], npx, ri, ritmp_t);
            do_range_interpolate(API, &c2_t[iy * npx], npx, ri, ritmp_t);
        }
    }

    const int total = npy * npx;
    double mean1 = 0.0, mean2 = 0.0;
    for (int i = 0; i < total; ++i) {
        float mag1 = Cabs(c1_t[i]);
        float mag2 = Cabs(c2_t[i]);
        c1_t[i].r = mag1; c1_t[i].i = 0.0f;
        c2_t[i].r = mag2; c2_t[i].i = 0.0f;
        mean1 += mag1;
        mean2 += mag2;
    }

    double inv_total = 1.0 / (double)total;
    mean1 *= inv_total;
    mean2 *= inv_total;

    for (int i = 0; i < total; ++i) {
        c1_t[i].r -= (float)mean1;
        c2_t[i].r -= (float)mean2;
        c1_t[i].i = c2_t[i].i = 0.0f;
        c2_t[i].r *= (float)mask_t[i];
        i1_t[i] = (int)c1_t[i].r;
        i2_t[i] = (int)c2_t[i].r;
    }
}

/* make_mask 简洁安全实现 */
void make_mask(struct xcorr *xc) {
    if (!xc || !xc->mask) die("make_mask: null", NULL);
    const int npy = xc->npy;
    const int npx = xc->npx;
    const int xsearch = xc->xsearch;
    const int ysearch = xc->ysearch;
    short *mask = xc->mask;

    const int y_min = ysearch;
    const int y_max = npy - ysearch;
    const int x_min = xsearch;
    const int x_max = npx - xsearch;

    for (int i = 0; i < npy; ++i) {
        int in_y = (i >= y_min && i < y_max);
        int base = i * npx;
        for (int j = 0; j < npx; ++j) {
            int idx = base + j;
            mask[idx] = (in_y && j >= x_min && j < x_max) ? 1 : 0;
        }
    }
}

/* allocate_arrays：为共享（旧接口）缓冲分配，保持原行为并检查失败 */
void allocate_arrays(struct xcorr *xc) {
    if (!xc) die("allocate_arrays: null", NULL);
    const size_t fcomplex_size = sizeof(struct FCOMPLEX);
    const size_t int_size = sizeof(int);
    const size_t short_size = sizeof(short);
    const size_t double_size = sizeof(double);

    xc->d1 = (struct FCOMPLEX *)malloc((size_t)xc->m_nx * xc->m_ny * fcomplex_size);
    xc->d2 = (struct FCOMPLEX *)malloc((size_t)xc->s_nx * xc->s_ny * fcomplex_size);

    xc->i1 = (int *)malloc((size_t)xc->npx * xc->npy * int_size);
    xc->i2 = (int *)malloc((size_t)xc->npx * xc->npy * int_size);

    xc->c1 = (struct FCOMPLEX *)malloc((size_t)xc->npx * xc->npy * fcomplex_size);
    xc->c2 = (struct FCOMPLEX *)malloc((size_t)xc->npx * xc->npy * fcomplex_size);
    xc->c3 = (struct FCOMPLEX *)malloc((size_t)xc->npx * xc->npy * fcomplex_size);

    xc->ritmp = (struct FCOMPLEX *)malloc((size_t)xc->ri * xc->npx * fcomplex_size);
    xc->mask = (short *)malloc((size_t)xc->npx * xc->npy * short_size);

    xc->corr = (double *)malloc((size_t)2 * xc->ri * xc->nxc * xc->nyc * double_size);

    if (xc->interp_flag == 1) {
        const int nx = 2 * xc->n2x;
        const int ny = 2 * xc->n2y;
        const int nx_exp = nx * xc->interp_factor;
        const int ny_exp = ny * xc->interp_factor;
        xc->md = (struct FCOMPLEX *)malloc((size_t)nx * ny * fcomplex_size);
        xc->cd_exp = (struct FCOMPLEX *)malloc((size_t)nx_exp * ny_exp * fcomplex_size);
    }

    if (!xc->d1 || !xc->d2 || !xc->i1 || !xc->i2 || !xc->c1 || !xc->c2 || !xc->c3 ||
        !xc->ritmp || !xc->mask || !xc->corr || (xc->interp_flag && (!xc->md || !xc->cd_exp))) {
        die("Memory allocation failed in allocate_arrays", NULL);
    }
}

/* do_correlation：并行实现，线程私有缓冲 + 最小读临界区 */
void do_correlation(void *API, struct xcorr *xc) {
    if (!xc) die("do_correlation: null xc", NULL);

    allocate_arrays(xc);
    make_mask(xc);

    const int nloc = xc->nloc;
    const int npx = xc->npx;
    const int npy = xc->npy;
    const int m_nx = xc->m_nx;
    const int s_nx = xc->s_nx;

    const size_t fsize = sizeof(struct FCOMPLEX);
    const size_t isize = sizeof(int);
    const size_t ssize = sizeof(short);

    int max_threads = omp_get_max_threads();
    if (max_threads < 1) max_threads = 1;

    /* 为每线程分配指针数组 */
    struct FCOMPLEX **d1_thr = calloc((size_t)max_threads, sizeof(struct FCOMPLEX *));
    struct FCOMPLEX **d2_thr = calloc((size_t)max_threads, sizeof(struct FCOMPLEX *));
    struct FCOMPLEX **c1_thr = calloc((size_t)max_threads, sizeof(struct FCOMPLEX *));
    struct FCOMPLEX **c2_thr = calloc((size_t)max_threads, sizeof(struct FCOMPLEX *));
    struct FCOMPLEX **c3_thr = calloc((size_t)max_threads, sizeof(struct FCOMPLEX *));
    struct FCOMPLEX **ritmp_thr = calloc((size_t)max_threads, sizeof(struct FCOMPLEX *));
    int **i1_thr = calloc((size_t)max_threads, sizeof(int *));
    int **i2_thr = calloc((size_t)max_threads, sizeof(int *));
    short **mask_thr = calloc((size_t)max_threads, sizeof(short *));

    if (!d1_thr || !d2_thr || !c1_thr || !c2_thr || !c3_thr || !ritmp_thr || !i1_thr || !i2_thr || !mask_thr) {
        /* 内存指针数组分配失败 -> 回退串行 */
        if (d1_thr) free(d1_thr); if (d2_thr) free(d2_thr);
        if (c1_thr) free(c1_thr); if (c2_thr) free(c2_thr); if (c3_thr) free(c3_thr);
        if (ritmp_thr) free(ritmp_thr);
        if (i1_thr) free(i1_thr); if (i2_thr) free(i2_thr); if (mask_thr) free(mask_thr);

        for (int iloc = 0; iloc < nloc; ++iloc) {
            read_xcorr_data(xc, iloc);
            assign_values_thread(API, xc, xc->d1, xc->d2, xc->c1, xc->c2, xc->c3, xc->i1, xc->i2, xc->mask, xc->ritmp, iloc);
            if (xc->corr_flag < 2) do_time_corr(xc, iloc);
            else if (xc->corr_flag == 2) do_freq_corr(API, xc, iloc);
            if (xc->interp_flag == 1) do_highres_corr(API, xc, iloc);
            print_results(xc, iloc);
        }
        return;
    }

    int allocation_failed = 0;
    for (int t = 0; t < max_threads; ++t) {
        d1_thr[t] = malloc((size_t)m_nx * xc->m_ny * fsize);
        d2_thr[t] = malloc((size_t)s_nx * xc->s_ny * fsize);
        c1_thr[t] = malloc((size_t)npx * npy * fsize);
        c2_thr[t] = malloc((size_t)npx * npy * fsize);
        c3_thr[t] = malloc((size_t)npx * npy * fsize);
        ritmp_thr[t] = malloc((size_t)xc->ri * npx * fsize);
        i1_thr[t] = malloc((size_t)npx * npy * isize);
        i2_thr[t] = malloc((size_t)npx * npy * isize);
        mask_thr[t] = malloc((size_t)npx * npy * ssize);

        if (!d1_thr[t] || !d2_thr[t] || !c1_thr[t] || !c2_thr[t] || !c3_thr[t] ||
            !ritmp_thr[t] || !i1_thr[t] || !i2_thr[t] || !mask_thr[t]) {
            allocation_failed = 1;
            for (int u = 0; u <= t; ++u) {
                if (d1_thr[u]) free(d1_thr[u]);
                if (d2_thr[u]) free(d2_thr[u]);
                if (c1_thr[u]) free(c1_thr[u]);
                if (c2_thr[u]) free(c2_thr[u]);
                if (c3_thr[u]) free(c3_thr[u]);
                if (ritmp_thr[u]) free(ritmp_thr[u]);
                if (i1_thr[u]) free(i1_thr[u]);
                if (i2_thr[u]) free(i2_thr[u]);
                if (mask_thr[u]) free(mask_thr[u]);
            }
            break;
        }
        memcpy(mask_thr[t], xc->mask, (size_t)npx * npy * ssize);
    }

    if (allocation_failed) {
        free(d1_thr); free(d2_thr); free(c1_thr); free(c2_thr); free(c3_thr);
        free(ritmp_thr); free(i1_thr); free(i2_thr); free(mask_thr);
        /* 回退串行 */
        for (int iloc = 0; iloc < nloc; ++iloc) {
            read_xcorr_data(xc, iloc);
            assign_values_thread(API, xc, xc->d1, xc->d2, xc->c1, xc->c2, xc->c3, xc->i1, xc->i2, xc->mask, xc->ritmp, iloc);
            if (xc->corr_flag < 2) do_time_corr(xc, iloc);
            else if (xc->corr_flag == 2) do_freq_corr(API, xc, iloc);
            if (xc->interp_flag == 1) do_highres_corr(API, xc, iloc);
            print_results(xc, iloc);
        }
        return;
    }

    /* 并行处理：读取在短临界区内，后续计算使用线程私有缓冲区 */
    #pragma omp parallel
    {
        int tid = omp_get_thread_num();
        struct FCOMPLEX *d1_t = d1_thr[tid];
        struct FCOMPLEX *d2_t = d2_thr[tid];
        struct FCOMPLEX *c1_t = c1_thr[tid];
        struct FCOMPLEX *c2_t = c2_thr[tid];
        struct FCOMPLEX *c3_t = c3_thr[tid];
        struct FCOMPLEX *ritmp_t = ritmp_thr[tid];
        int *i1_t = i1_thr[tid];
        int *i2_t = i2_thr[tid];
        short *mask_t = mask_thr[tid];

        #pragma omp for schedule(dynamic)
        for (int iloc = 0; iloc < nloc; ++iloc) {
            /* 只在临界区读取共享数据到共享缓冲 xc->d1/xc->d2，然后 memcpy 到线程私有缓冲 */
            #pragma omp critical(read_data)
            {
                read_xcorr_data(xc, iloc);
                memcpy(d1_t, xc->d1, (size_t)m_nx * xc->m_ny * fsize);
                memcpy(d2_t, xc->d2, (size_t)s_nx * xc->s_ny * fsize);
            }

            assign_values_thread(API, xc, d1_t, d2_t, c1_t, c2_t, c3_t, i1_t, i2_t, mask_t, ritmp_t, iloc);

            if (xc->corr_flag < 2) do_time_corr(xc, iloc);
            else if (xc->corr_flag == 2) do_freq_corr(API, xc, iloc);

            if (xc->interp_flag == 1) do_highres_corr(API, xc, iloc);

            #pragma omp critical(write_results)
            {
                print_results(xc, iloc);
            }
        } /* end for iloc */
    } /* end parallel */

    /* 释放线程私有内存 */
    for (int t = 0; t < max_threads; ++t) {
        free(d1_thr[t]); free(d2_thr[t]);
        free(c1_thr[t]); free(c2_thr[t]); free(c3_thr[t]);
        free(ritmp_thr[t]);
        free(i1_thr[t]); free(i2_thr[t]);
        free(mask_thr[t]);
    }
    free(d1_thr); free(d2_thr); free(c1_thr); free(c2_thr); free(c3_thr);
    free(ritmp_thr); free(i1_thr); free(i2_thr); free(mask_thr);
}

/* main：仅作少量保护和资源释放 */
int main(int argc, char **argv) {
    int input_flag = 0, nfiles = 2;
    struct xcorr *xc = NULL;
    clock_t start, end;
    double cpu_time;
    void *API = NULL;

#ifdef _OPENMP
    printf("=== OpenMP enabled ===\n");
    int nthreads = omp_get_max_threads();
    if (nthreads < 1) nthreads = 1;
    printf("max threads available: %d\n", nthreads);
#else
    printf("=== Warning: OpenMP not enabled; running serial ===\n");
#endif

    xc = (struct xcorr *)malloc(sizeof(struct xcorr));
    if (!xc) die("Memory allocation failed for xcorr struct", NULL);

    verbose = 0;
    debug = 0;
    xc->interp_flag = 0;
    xc->corr_flag = 2;

    if ((API = GMT_Create_Session(argv[0], 0U, 0U, NULL)) == NULL) {
        free(xc);
        return EXIT_FAILURE;
    }

    if (argc < 3) die(USAGE, "");

    set_defaults(xc);
    parse_command_line(argc, argv, xc, &nfiles, &input_flag, USAGE);

    if (input_flag == 0) handle_prm(API, argv, xc, nfiles);

    if (xc->corr_flag == 0) strcpy(xc->filename, "time_xcorr.dat");
    else if (xc->corr_flag == 1) strcpy(xc->filename, "time_xcorr_Gatelli.dat");
    else strcpy(xc->filename, "freq_xcorr.dat");

    xc->file = fopen(xc->filename, "w");
    if (!xc->file) die("Can't open output file", xc->filename);

    get_locations(xc);

    start = clock();
    do_correlation(API, xc);
    end = clock();

    cpu_time = ((double)(end - start)) / CLOCKS_PER_SEC;
    fprintf(stdout, " elapsed time: %lf \n", cpu_time);

    if (xc->format == 0 || xc->format == 1) {
        if (xc->data1) fclose(xc->data1);
        if (xc->data2) fclose(xc->data2);
    }
    if (xc->file) fclose(xc->file);

    GMT_Destroy_Session(API);

    if (xc->d1) free(xc->d1);
    if (xc->d2) free(xc->d2);
    if (xc->i1) free(xc->i1);
    if (xc->i2) free(xc->i2);
    if (xc->c1) free(xc->c1);
    if (xc->c2) free(xc->c2);
    if (xc->c3) free(xc->c3);
    if (xc->ritmp) free(xc->ritmp);
    if (xc->mask) free(xc->mask);
    if (xc->corr) free(xc->corr);
    if (xc->interp_flag == 1) {
        if (xc->md) free(xc->md);
        if (xc->cd_exp) free(xc->cd_exp);
    }
    if (xc->loc) free(xc->loc);
    free(xc);

    return EXIT_SUCCESS;
}
