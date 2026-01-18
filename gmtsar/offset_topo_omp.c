#include "gmtsar.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>

#ifdef _OPENMP
#include <omp.h>
#else
// 如果OpenMP未启用，提供空定义
#define omp_get_thread_num() 0
#define omp_get_max_threads() 1
#define omp_get_num_threads() 1
#define omp_set_num_threads(n)
#endif

// 辅助宏函数
#ifndef max
#define max(a,b) ((a) > (b) ? (a) : (b))
#endif
#ifndef min
#define min(a,b) ((a) < (b) ? (a) : (b))
#endif

// 计时函数
double get_time() {
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return tv.tv_sec + tv.tv_usec * 1e-6;
}

int main(int argc, char **argv) {
    double total_start = get_time();
    
    int i, j, is, js;
    int ni, nj, ntot;
    int xshft, yshft, ib = 200;
    int ns;
    int imax = 0, jmax = 0;
    double avea;
    double suma, sumt, maxcorr = -1e30;

    // 检查OpenMP支持
#ifdef _OPENMP
    printf("=== OpenMP 已启用 ===\n");
    // 设置线程数（使用所有可用核心）
	printf("最大可用线程数: %d\n", omp_get_max_threads());
    int nthreads = omp_get_max_threads()*5/6;
	if (nthreads < 2) nthreads = 2; 
    omp_set_num_threads(nthreads);
    printf("使用 %d 个线程进行计算\n", nthreads);
#else
    printf("=== 警告: OpenMP 未启用，程序将串行运行 ===\n");
    printf("编译时请添加 -fopenmp 选项启用并行计算\n");
#endif

    void *API = NULL;
    struct GMT_GRID *A = NULL, *T = NULL, *TS = NULL;

    if (argc < 6) {
        fprintf(stderr,
                "用法: offset_topo2 amp_master.grd topo_ra.grd rshift ashift ns [topo_shift.grd]\n");
        exit(EXIT_FAILURE);
    }

    API = GMT_Create_Session(argv[0], 0U, 0U, NULL);

    xshft = atoi(argv[3]);
    yshft = atoi(argv[4]);
    ns = atoi(argv[5]);

    // 读取网格数据
    double read_start = get_time();
    A = GMT_Read_Data(API, GMT_IS_GRID, GMT_IS_FILE, GMT_IS_SURFACE,
                      GMT_GRID_HEADER_ONLY, NULL, argv[1], NULL);
    T = GMT_Read_Data(API, GMT_IS_GRID, GMT_IS_FILE, GMT_IS_SURFACE,
                      GMT_GRID_HEADER_ONLY, NULL, argv[2], NULL);

    if (A->header->n_columns != T->header->n_columns) {
        fprintf(stderr, "错误: 网格宽度不匹配\n");
        exit(EXIT_FAILURE);
    }

    // 读取实际数据
    GMT_Read_Data(API, GMT_IS_GRID, GMT_IS_FILE, GMT_IS_SURFACE,
                  GMT_GRID_DATA_ONLY, NULL, argv[1], A);
    GMT_Read_Data(API, GMT_IS_GRID, GMT_IS_FILE, GMT_IS_SURFACE,
                  GMT_GRID_DATA_ONLY, NULL, argv[2], T);
    printf("数据读取完成: %.3f秒\n", get_time() - read_start);

    ni = min(A->header->n_rows, T->header->n_rows);
    nj = T->header->n_columns;
    printf("网格尺寸: %d x %d, 总像素数: %d\n", ni, nj, ni * nj);

    // 为输出网格分配内存（如果需要）
    if (argc >= 7) {
        TS = GMT_Create_Data(API, GMT_IS_GRID, GMT_IS_SURFACE, GMT_GRID_ALL,
                             NULL, A->header->wesn, A->header->inc,
                             A->header->registration, GMT_NOTSET, NULL);
    }

    /* ---------- 1. 预计算阶段 ---------- */
    double precompute_start = get_time();
    
    // 计算A的均值
    suma = 0.0;
    ntot = ni * nj;
    
    // 并行计算均值
#ifdef _OPENMP
    #pragma omp parallel for reduction(+:suma)
#endif
    for (i = 0; i < ni; i++) {
        int row_start = i * nj;
        double row_sum = 0.0;
        for (j = 0; j < nj; j++) {
            row_sum += A->data[row_start + j];
        }
        suma += row_sum;
    }
    avea = suma / ntot;
    printf("均值计算完成: avea = %g\n", avea);

    // 分配预计算数组
    float *A_centered = (float*)malloc(ni * nj * sizeof(float));
    float *T_grad_x = (float*)malloc(ni * nj * sizeof(float));
    
    if (!A_centered || !T_grad_x) {
        fprintf(stderr, "内存分配失败\n");
        exit(EXIT_FAILURE);
    }

    // 预计算A的中心化值
#ifdef _OPENMP
    #pragma omp parallel for collapse(2)
#endif
    for (i = 0; i < ni; i++) {
        for (j = 0; j < nj; j++) {
            int idx = i * nj + j;
            A_centered[idx] = (float)(A->data[idx] - avea);
        }
    }

    // 预计算T的x方向梯度
#ifdef _OPENMP
    #pragma omp parallel for collapse(2)
#endif
    for (i = 0; i < ni; i++) {
        for (j = 1; j < nj - 1; j++) {
            int idx = i * nj + j;
            T_grad_x[idx] = (float)(T->data[idx + 1] - T->data[idx - 1]);
        }
    }
    
    // 边界处理
#ifdef _OPENMP
    #pragma omp parallel for
#endif
    for (i = 0; i < ni; i++) {
        T_grad_x[i * nj] = 0.0f;
        T_grad_x[i * nj + nj - 1] = 0.0f;
    }
    
    printf("预计算完成: %.3f秒\n", get_time() - precompute_start);

    /* ---------- 2. 互相关搜索 ---------- */
    double corr_start = get_time();
    
    // 计算有效搜索区域
    int valid_start_i = max(ib, ns);
    int valid_end_i = ni - max(ib, ns);
    int valid_start_j = max(ib, ns);
    int valid_end_j = nj - max(ib, ns);
    
    int is_start = -ns + yshft;
    int is_end = ns + yshft;
    int js_start = -ns + xshft;
    int js_end = ns + xshft;
    
    printf("搜索范围: is=[%d, %d], js=[%d, %d]\n", 
           is_start, is_end, js_start, js_end);
    printf("有效计算区域: i=[%d, %d], j=[%d, %d]\n",
           valid_start_i, valid_end_i, valid_start_j, valid_end_j);

    // 主搜索循环 - 使用更简单的并行策略
#ifdef _OPENMP
    #pragma omp parallel
#endif
    {
#ifdef _OPENMP
        double local_maxcorr = -1e30;
        int local_imax = 0, local_jmax = 0;
        
        // 使用static调度，每个线程处理不同的位移
        #pragma omp for collapse(2) schedule(static)
#endif
        for (is = is_start; is <= is_end; is++) {
            for (js = js_start; js <= js_end; js++) {
                
                double sumc = 0.0, sumaa = 0.0, sumtt = 0.0;
                
                // 计算当前位移的有效区域
                int start_i = max(valid_start_i, abs(is));
                int end_i = valid_end_i;
                if (is < 0) end_i = min(end_i, ni + is);
                
                int start_j = max(valid_start_j, abs(js));
                int end_j = valid_end_j;
                if (js < 0) end_j = min(end_j, nj + js);
                
                if (start_i >= end_i || start_j >= end_j) continue;
                
                // 计算相关系数
                for (i = start_i; i < end_i; i++) {
                    int i1 = i - is;
                    long base_idx = i * nj;
                    long base_idx1 = i1 * nj;
                    
                    for (j = start_j; j < end_j; j++) {
                        int j1 = j - js;
                        int idx = base_idx + j;
                        int idx1 = base_idx1 + j1;
                        
                        float ra = A_centered[idx];
                        float rt = T_grad_x[idx1];
                        
                        sumc += ra * rt;
                        sumaa += ra * ra;
                        sumtt += rt * rt;
                    }
                }
                
                // 计算相关系数
                double denom = sumaa * sumtt;
                if (denom > 0.0) {
                    double corr = sumc / sqrt(denom);
#ifdef _OPENMP
                    if (corr > local_maxcorr) {
                        local_maxcorr = corr;
                        local_imax = is;
                        local_jmax = js;
                    }
#else
                    if (corr > maxcorr) {
                        maxcorr = corr;
                        imax = is;
                        jmax = js;
                    }
#endif
                }
            }
        }
        
#ifdef _OPENMP
        // 更新全局最大值
        #pragma omp critical
        {
            if (local_maxcorr > maxcorr) {
                maxcorr = local_maxcorr;
                imax = local_imax;
                jmax = local_jmax;
            }
        }
#endif
    }
    
    printf("互相关搜索完成: %.3f秒\n", get_time() - corr_start);
    printf("最优结果: rshift=%d ashift=%d maxcorr=%g\n", jmax, imax, maxcorr);

    /* ---------- 3. 输出移位后的DEM ---------- */
    if (argc >= 7) {
        double output_start = get_time();
        
#ifdef _OPENMP
        #pragma omp parallel for collapse(2)
#endif
        for (i = 0; i < ni; i++) {
            for (j = 0; j < nj; j++) {
                int i1 = i - imax;
                int j1 = j - jmax;
                int idx = i * nj + j;
                
                if (i1 >= 0 && i1 < ni && j1 >= 0 && j1 < nj) {
                    int idx1 = i1 * nj + j1;
                    TS->data[idx] = T->data[idx1];
                } else {
                    TS->data[idx] = 0.0f;
                }
            }
        }
        
        GMT_Write_Data(API, GMT_IS_GRID, GMT_IS_FILE, GMT_IS_SURFACE,
                       GMT_GRID_ALL, NULL, argv[6], TS);
        
        printf("输出文件已保存: %s (%.3f秒)\n", argv[6], get_time() - output_start);
    }

    /* ---------- 4. 清理内存 ---------- */
    free(A_centered);
    free(T_grad_x);
    
    GMT_Destroy_Session(API);
    
    printf("总运行时间: %.3f秒\n", get_time() - total_start);
    
    return EXIT_SUCCESS;
}