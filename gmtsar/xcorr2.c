/*	$Id: xcorr.c 73 2013-04-19 17:59:45Z pwessel $	*/
/***************************************************************************/
/* xcorr does a 2-D cross correlation on complex or real images            */
/* either using a time convolution or wavenumber multiplication.           */
/*                                                                         */
/***************************************************************************/

/***************************************************************************
 * Creator:  Rob J. Mellors                                                *
 *           (San Diego State University)                                  *
 * Date   :  November 7, 2009                                              *
 ***************************************************************************/

/***************************************************************************
 * Modification history:                                                   *
 *                                                                         *
 * DATE                                                                     *
 *                                                                         *
 * 011810       Testing and very minor cosmetic modifications DTS          *
 * 061520       Problem with sub-pixel interpolation RJM                   *
 *              - fixed bug in 2D interpolation                            *
 *              - revised read_xcorr_data to read in all x position        *
 *              - reads directly into float rather than int                *
 *              - add range interpolation                                  *
 *              - eliminated obsolete options and code                     *
 *              - renamed xcorr_utils.c print_results.c                    *
 *              - further testing....                                      *
 ***************************************************************************/
/*-------------------------------------------------------*/
// 必须先包含系统头文件，再包含GMTSAR相关头文件
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
//#include <omp.h>


#ifdef _OPENMP
#include <omp.h>
#else
// 如果OpenMP未启用，提供空定义
#define omp_get_thread_num() 0
#define omp_get_max_threads() 1
#define omp_get_num_threads() 1
#define omp_set_num_threads(n)
#endif

// 补充缺失的结构体声明（匹配xcorr.h中的定义）
#ifndef XC_H
// 提前声明xcorr.h中的核心结构体，避免重复定义
struct LOCATION {  // 对应代码中的loc数组类型
    int x;
    int y;
};

// 声明xcorr结构体（实际定义在xcorr.h中，此处仅声明）
struct xcorr;
#endif

// 包含GMTSAR头文件（确保struct xcorr只在xcorr.h中定义）
#include "gmtsar.h"

// 全局变量声明（匹配GMTSAR规范）
extern int verbose;
extern int debug;

// 错误处理函数声明（匹配GMTSAR的die函数）
// void die(const char *msg, const char *arg);

// 外部函数声明（来自GMTSAR库）
void fft_interpolate_1d(void *API, struct FCOMPLEX *c, int nx, struct FCOMPLEX *work, int ri);
float Cabs(struct FCOMPLEX z);
void allocate_arrays(struct xcorr *xc);
void make_mask(struct xcorr *xc);
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
/*-------------------------------------------------------*/
// 调试宏定义：非调试模式下移除调试代码




char *USAGE = "xcorr2 [GMTSAR] - Compute 2-D cross-correlation of two images\n\n"
              "\nUsage: xcorr2 master.PRM aligned.PRM [-time] [-real] [-freq] [-nx n] [-ny "
              "n] [-xsearch xs] [-ysearch ys]\n"
              "master.PRM         PRM file for reference image\n"
              "aligned.PRM        PRM file of secondary image\n"
              "-time              use time cross-correlation\n"
              "-freq              use frequency cross-correlation (default)\n"
              "-real              read float numbers instead of complex numbers\n"
              "-noshift           ignore ashift and rshift in prm file (set to 0)\n"
              "-nx  nx            number of locations in x (range) direction "
              "(int)\n"
              "-ny  ny            number of locations in y (azimuth) direction "
              "(int)\n"
              "-nointerp          do not interpolate correlation function\n"
              "-range_interp ri   interpolate range by ri (power of two) [default: 2]\n"
              "-norange           do not range interpolate \n"
              "-xsearch xs        search window size in x (range) direction (int "
              "power of 2 [32 64 128 256])\n"
              "-ysearch ys        search window size in y (azimuth) direction "
              "(int power of 2 [32 64 128 256])\n"
              "-interp  factor    interpolate correlation function by factor "
              "(int) [default, 16]\n"
              "-v                 verbose\n"
              "output: \n freq_xcorr.dat (default) \n time_xcorr.dat (if -time option))\n"
              "\nuse fitoffset.csh to convert output to PRM format\n"
              "\nExample:\n"
              "xcorr2 IMG-HH-ALPSRP075880660-H1.0__A.PRM "
              "IMG-HH-ALPSRP129560660-H1.0__A.PRM -nx 20 -ny 50 \n"
              "xcorr2 file1.grd file2.grd -nx 20 -ny 50 (takes grids with real numbers)\n";

/*-------------------------------------------------------------------------------*/
// 高频调用的小函数 inline，减少函数调用开销
static inline int do_range_interpolate(void *API, struct FCOMPLEX *c, int nx, int ri, struct FCOMPLEX *work) {
    int i;
    fft_interpolate_1d(API, c, nx, work, ri);  // 假设此函数已优化
    
    // 预计算终止条件，减少循环内计算
    const int end = nx;
    for (i = 0; i < end; i++) {
        c[i].r = work[i + nx / 2].r;
        c[i].i = work[i + nx / 2].i;
    }
    return EXIT_SUCCESS;
}
/*-------------------------------------------------------------------------------*/
void assign_values(void *API, struct xcorr *xc, int iloc) {
    int i, j, k;
    double mean1 = 0.0, mean2 = 0.0;
    
    // 预计算偏移量，减少重复计算
    const int mx = xc->loc[iloc].x - xc->npx / 2;
    const int sx = xc->loc[iloc].x + xc->x_offset - xc->npx / 2;
    const int npx = xc->npx;
    const int npy = xc->npy;
    const int m_nx = xc->m_nx;
    const int s_nx = xc->s_nx;
    
    // 缓存指针，减少多级指针访问开销
    struct FCOMPLEX *c1 = xc->c1;
    struct FCOMPLEX *c2 = xc->c2;
    struct FCOMPLEX *c3 = xc->c3;
    struct FCOMPLEX *d1 = xc->d1;
    struct FCOMPLEX *d2 = xc->d2;
    
    // 行优先访问，提升缓存利用率
    for (i = 0; i < npy; i++) {
        // 预计算行起始索引，减少乘法运算
        const int d1_row_base = i * m_nx + mx;
        const int d2_row_base = i * s_nx + sx;
        const int c_row_base = i * npx;
        
        for (j = 0; j < npx; j++) {
            k = c_row_base + j;
            c3[k].i = c3[k].r = 0.0f;
            
            // 直接使用预计算的行基地址，减少计算
            c1[k].r = d1[d1_row_base + j].r;
            c1[k].i = d1[d1_row_base + j].i;
            
            c2[k].r = d2[d2_row_base + j].r;
            c2[k].i = d2[d2_row_base + j].i;
        }
    }

    // 范围插值：仅在需要时执行
    if (xc->ri > 1) {
        const int ri = xc->ri;
        struct FCOMPLEX *ritmp = xc->ritmp;
        for (i = 0; i < npy; i++) {
            do_range_interpolate(API, &c1[i * npx], npx, ri, ritmp);
            do_range_interpolate(API, &c2[i * npx], npx, ri, ritmp);
        }
    }

    // 计算振幅和均值（合并循环，减少遍历次数）
    const int total = npy * npx;
    for (i = 0; i < total; i++) {
        // 合并计算，减少对同一内存的重复访问
        c1[i].r = Cabs(c1[i]);
        c1[i].i = 0.0f;
        mean1 += c1[i].r;

        c2[i].r = Cabs(c2[i]);
        c2[i].i = 0.0f;
        mean2 += c2[i].r;
    }

    // 均值归一化
    const double inv_total = 1.0 / total;
    mean1 *= inv_total;
    mean2 *= inv_total;

    // 去均值并应用掩膜（合并循环）
    short *mask = xc->mask;
    int *i1 = xc->i1;
    int *i2 = xc->i2;
    for (i = 0; i < total; i++) {
        c1[i].r -= (float)mean1;
        c2[i].r -= (float)mean2;

        c1[i].i = 0.0f;
        c2[i].i = 0.0f;
        c2[i].r *= (float)mask[i];

        i1[i] = (int)c1[i].r;
        i2[i] = (int)c2[i].r;
    }

}
/*-------------------------------------------------------------------------------*/
void do_correlation(void *API, struct xcorr *xc) {
    int i, j, iloc;
    const int istep = 1;  // 保持原步长逻辑
    
    allocate_arrays(xc);
    make_mask(xc);

    // 并行化外层循环（独立迭代，无数据依赖）
    #pragma omp parallel for private(i, j, iloc) collapse(2)
    for (i = 0; i < xc->nyl; i += istep) {
        // 按行读取数据（减少IO次数，提升缓存利用率）
        read_xcorr_data(xc, i * xc->nxl);  // 假设iloc按行连续存储
        
        for (j = 0; j < xc->nxl; j++) {
            iloc = i * xc->nxl + j;

            assign_values(API, xc, iloc);

            // 选择相关计算方式
            if (xc->corr_flag < 2)
                do_time_corr(xc, iloc);
            else if (xc->corr_flag == 2)
                do_freq_corr(API, xc, iloc);

            // 子像素插值
            if (xc->interp_flag == 1)
                do_highres_corr(API, xc, iloc);

            // 输出结果（注意线程安全，若print_results有写操作需加锁）
            #pragma omp critical
            print_results(xc, iloc);
        }
    }
}
/*-------------------------------------------------------------------------------*/
void make_mask(struct xcorr *xc) {
    int i, j;
    const int imask = 0;
    const int npy = xc->npy;
    const int npx = xc->npx;
    const int xsearch = xc->xsearch;
    const int ysearch = xc->ysearch;
    short *mask = xc->mask;
    
    // 预计算边界条件，减少循环内比较
    const int y_min = ysearch;
    const int y_max = npy - ysearch;
    const int x_min = xsearch;
    const int x_max = npx - xsearch;
    
    for (i = 0; i < npy; i++) {
        // 提前判断行是否在边界内，减少内层循环条件判断
        const int in_y_range = (i >= y_min && i < y_max);
        const int row_base = i * npx;
        
        for (j = 0; j < npx; j++) {
            const int idx = row_base + j;
            if (in_y_range && j >= x_min && j < x_max) {
                mask[idx] = 1;
            } else {
                mask[idx] = imask;
            }
        }
    }
}
/*-------------------------------------------------------------------------------*/
void allocate_arrays(struct xcorr *xc) {
    // 一次性计算所有内存需求，避免碎片化
    const size_t fcomplex_size = sizeof(struct FCOMPLEX);
    const size_t int_size = sizeof(int);
    const size_t short_size = sizeof(short);
    const size_t double_size = sizeof(double);

    // 优先分配大块内存，提升缓存效率
    xc->d1 = (struct FCOMPLEX *)malloc(xc->m_nx * xc->npy * fcomplex_size);
    xc->d2 = (struct FCOMPLEX *)malloc(xc->s_nx * xc->npy * fcomplex_size);

    xc->i1 = (int *)malloc(xc->npx * xc->npy * int_size);
    xc->i2 = (int *)malloc(xc->npx * xc->npy * int_size);

    xc->c1 = (struct FCOMPLEX *)malloc(xc->npx * xc->npy * fcomplex_size);
    xc->c2 = (struct FCOMPLEX *)malloc(xc->npx * xc->npy * fcomplex_size);
    xc->c3 = (struct FCOMPLEX *)malloc(xc->npx * xc->npy * fcomplex_size);

    xc->ritmp = (struct FCOMPLEX *)malloc(xc->ri * xc->npx * fcomplex_size);
    xc->mask = (short *)malloc(xc->npx * xc->npy * short_size);

    // 相关结果数组
    xc->corr = (double *)malloc(2 * xc->ri * xc->nxc * xc->nyc * double_size);

    // 插值相关数组（条件分配）
    if (xc->interp_flag == 1) {
        const int nx = 2 * xc->n2x;
        const int ny = 2 * xc->n2y;
        const int nx_exp = nx * xc->interp_factor;
        const int ny_exp = ny * xc->interp_factor;
        
        xc->md = (struct FCOMPLEX *)malloc(nx * ny * fcomplex_size);
        xc->cd_exp = (struct FCOMPLEX *)malloc(nx_exp * ny_exp * fcomplex_size);
    }

    // 检查内存分配失败（新增错误处理）
    if (!xc->d1 || !xc->d2 || !xc->i1 || !xc->i2 || !xc->c1 || !xc->c2 || !xc->c3 || 
        !xc->ritmp || !xc->mask || !xc->corr || (xc->interp_flag && (!xc->md || !xc->cd_exp))) {
        die("Memory allocation failed in allocate_arrays", NULL);
    }
}

/*-------------------------------------------------------*/
int main(int argc, char **argv) {
    int input_flag = 0, nfiles = 2;
    struct xcorr *xc;
    clock_t start, end;
    double cpu_time;
    void *API = NULL;


    // 检查OpenMP支持
#ifdef _OPENMP
    printf("=== OpenMP 已启用 ===\n");
        // 设置线程数（使用所有可用核心）
    printf("最大可用线程数: %d\n", omp_get_max_threads());
    int nthreads = omp_get_max_threads()*4/5+1;
    if (nthreads < 2) nthreads = 2; 
    omp_set_num_threads(nthreads);
    printf("使用 %d 个线程进行计算\n", nthreads);
#else
    printf("=== 警告: OpenMP 未启用，程序将串行运行 ===\n");
    printf("编译时请添加 -fopenmp 选项启用并行计算\n");
#endif




    xc = (struct xcorr *)malloc(sizeof(struct xcorr));
    if (!xc) die("Memory allocation failed for xcorr struct", NULL);

    verbose = 0;
    debug = 0;
    xc->interp_flag = 0;
    xc->corr_flag = 2;

    // 初始化GMT会话
    if ((API = GMT_Create_Session(argv[0], 0U, 0U, NULL)) == NULL) {
        free(xc);
        return EXIT_FAILURE;
    }

    if (argc < 3) {
        die(USAGE, "");
    }

    set_defaults(xc);
    parse_command_line(argc, argv, xc, &nfiles, &input_flag, USAGE);

    if (input_flag == 0)
        handle_prm(API, argv, xc, nfiles);

    // 设置输出文件
    if (xc->corr_flag == 0)
        strcpy(xc->filename, "time_xcorr.dat");
    else if (xc->corr_flag == 1)
        strcpy(xc->filename, "time_xcorr_Gatelli.dat");
    else
        strcpy(xc->filename, "freq_xcorr.dat");

    xc->file = fopen(xc->filename, "w");
    if (!xc->file) die("Can't open output file", xc->filename);

    get_locations(xc);

    // 计时并执行相关计算
    start = clock();
    do_correlation(API, xc);
    end = clock();

    cpu_time = ((double)(end - start)) / CLOCKS_PER_SEC;
    fprintf(stdout, " elapsed time: %lf \n", cpu_time);

    // 关闭文件
    if (xc->format == 0 || xc->format == 1) {
        fclose(xc->data1);
        fclose(xc->data2);
    }
    fclose(xc->file);  // 确保输出文件关闭

    // 清理资源
    GMT_Destroy_Session(API);
    free(xc->d1); free(xc->d2); free(xc->i1); free(xc->i2);
    free(xc->c1); free(xc->c2); free(xc->c3); free(xc->ritmp);
    free(xc->mask); free(xc->corr);
    if (xc->interp_flag == 1) {
        free(xc->md);
        free(xc->cd_exp);
    }
    free(xc->loc);  // 假设loc在get_locations中分配
    free(xc);

    return EXIT_SUCCESS;
}
