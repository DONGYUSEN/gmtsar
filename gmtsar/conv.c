/*	$Id: conv.c 109 2015-01-19 23:01:24Z sandwell $	*/
/***************************************************************************/
/* conv convolves a 2-D filter with an array and outputs the results       */
/***************************************************************************/

#include "gmtsar.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#ifndef min
#define min(a,b) ((a) < (b) ? (a) : (b))
#endif


char *USAGE = "conv [GMTSAR] - 2-D image convolution\n\n"
              "Usage: conv idec jdec filter_file input output \n";

int input_file_type, format_flag;

/*-------------------------------------------------------------*/
int determine_file_type(char *name, int *input_file_type) {
    int n, m;
    char tail[8];

    *input_file_type = 1;

    n = (int)strlen(name);
    m = n - 3;
    strncpy(&tail[0], &name[m], 4);

    if ((strncmp(tail, "PRM", 3) == 0) || (strncmp(tail, "prm", 3) == 0)) {
        *input_file_type = 2;
    }

    return (EXIT_SUCCESS);
}

/*-------------------------------------------------------------*/
FILE *read_PRM_file(char *prmfilename, char *input_file_name, struct PRM p, int *xdim, int *ydim) {
    FILE *f_input_prm, *f_input;

    if ((f_input_prm = fopen(prmfilename, "r")) == NULL)
        die("Can't open input header", prmfilename);
    null_sio_struct(&p);
    get_sio_struct(f_input_prm, &p);
    strcpy(input_file_name, p.SLC_file);
    format_flag = 2;
    if (strncmp(p.dtype, "c", 1) == 0)
        format_flag = 3;
    if ((f_input = fopen(input_file_name, "r")) == NULL)
        die("Can't open input data ", input_file_name);
    *xdim = p.num_rng_bins;
    *ydim = p.num_valid_az * p.num_patches;

    return (f_input);
}

/*-------------------------------------------------------------*/
int read_float(float *indat, int xdim, FILE *f_input, int yarr, float *buffer, int ibuff) {
    int i, j;
    for (i = 0; i < ibuff; i++) {
        fread(indat, sizeof(float), xdim, f_input);
        for (j = 0; j < xdim; j++)
            buffer[j + xdim * (i + yarr)] = indat[j];
    }
    return (EXIT_SUCCESS);
}

int read_SLC_int(short *ci2, int xdim, FILE *f_input, int yarr, float *buffer, double dfact, int ibuff) {
    int i, j;
    double df2 = dfact * dfact;
    for (i = 0; i < ibuff; i++) {
        fread(ci2, 2 * sizeof(short), xdim, f_input);
        for (j = 0; j < xdim; j++)
            buffer[j + xdim * (i + yarr)] =
                (float)(df2 * ci2[2 * j] * ci2[2 * j] + df2 * ci2[2 * j + 1] * ci2[2 * j + 1]);
    }
    return (EXIT_SUCCESS);
}

int read_SLC_float(float *cf2, int xdim, FILE *f_input, int yarr, float *buffer, double dfact, int ibuff) {
    int i, j;
    double df2 = dfact * dfact;
    for (i = 0; i < ibuff; i++) {
        fread(cf2, 2 * sizeof(float), xdim, f_input);
        for (j = 0; j < xdim; j++)
            buffer[j + xdim * (i + yarr)] =
                (float)(df2 * cf2[2 * j] * cf2[2 * j] + df2 * cf2[2 * j + 1] * cf2[2 * j + 1]);
    }
    return (EXIT_SUCCESS);
}

/*-------------------------------------------------------------*/
int main(int argc, char **argv) {
    int idec, jdec;
    int iout, jout;
    int i, j, ic, jc, norm, ic0, ic1;
    int ydim = 0, xdim = 0;
    int xarr, yarr, narr, yarr2;
    int nbuff, ibuff, imove;
    int iend, ylen, iread;
    uint64_t left_node;
    unsigned int row;
    char input_name[128], output_name[128], prmfilename[128], *c = NULL;
    short *cindat = NULL;
    float *cfdat = NULL;
    double inc[2], wesn[4], xmax = 0.0, ymax = 0.0;
    float *filter = NULL, *buffer = NULL, *indat = NULL;
    float filtin, rnormax, anormax;
    FILE *f_filter = NULL, *f_input = NULL;
    struct PRM p;
    void *API = NULL;
    struct GMT_GRID *Out = NULL;
    struct GMT_GRID *In = NULL;
    
    // 新增：输出网格的行列数
    int output_n_columns, output_n_rows;

    if (argc < 6)
        die("\n", USAGE);

    if ((API = GMT_Create_Session(argv[0], 0U, 0U, NULL)) == NULL)
        return EXIT_FAILURE;

    // 检查OpenMP支持
#ifdef _OPENMP
    printf("=== OpenMP 已启用 ===\n");
    printf("最大可用线程数: %d\n", omp_get_max_threads());
    int nthreads = omp_get_max_threads() * 4 / 5 + 1;
    if (nthreads < 2) nthreads = 2; 
    omp_set_num_threads(nthreads);
    printf("使用 %d 个线程进行计算\n", nthreads);
#else
    printf("=== 警告: OpenMP 未启用，程序将串行运行 ===\n");
    printf("编译时请添加 -fopenmp 选项启用并行计算\n");
#endif

    ibuff = 512;
    verbose = 0;

    null_sio_struct(&p);
    input_file_type = 1;
    format_flag = 1;

    idec = atoi(argv[1]);
    jdec = atoi(argv[2]);

    if ((f_filter = fopen(argv[3], "r")) == NULL)
        die("Can't open filter", "");

    strcpy(input_name, argv[4]);
    strcpy(output_name, argv[5]);

    determine_file_type(input_name, &input_file_type);

    switch (input_file_type) {
    case 1:
        if ((In = GMT_Read_Data(API, GMT_IS_GRID, GMT_IS_FILE, GMT_IS_SURFACE,
                                GMT_GRID_HEADER_ONLY, NULL, input_name, NULL)) == NULL)
            die("Can't open ", input_name);
        if ((c = strstr(input_name, "=bf")))
            c[0] = '\0';
        if ((f_input = fopen(input_name, "r")) == NULL)
            die("Can't open ", input_name);
        fseek(f_input, 892L, SEEK_SET);
        xdim = In->header->n_columns;
        ydim = In->header->n_rows;
        xmax = In->header->wesn[GMT_XHI];
        ymax = In->header->wesn[GMT_YHI];
        format_flag = 1;
        break;
    case 2:
        strcpy(prmfilename, input_name);
        f_input = read_PRM_file(prmfilename, input_name, p, &xdim, &ydim);
        xmax = xdim;
        ymax = ydim;
        break;
    default:
        die("confused about input file type", "quitting");
    }

    if (fscanf(f_filter, "%d%d", &xarr, &yarr) != 2)
        die("filter incomplete", "");

    narr = xarr * yarr;
    yarr2 = yarr / 2;

    if ((filter = (float *)malloc(sizeof(float) * narr)) == NULL)
        die("memory allocation", "");
    if ((buffer = (float *)malloc(2 * sizeof(float) * xdim * ibuff)) == NULL)
        die("memory allocation", "");

    anormax = rnormax = 0.0f;
    for (i = 0; i < narr; i++) {
        if (fscanf(f_filter, "%f", &filtin) == EOF)
            die("filter incomplete", "");
        filter[i] = filtin;
        anormax += fabs(filter[i]);
        rnormax += filter[i];
    }

    norm = 0;
    if (fabs(rnormax) > 0.05 * anormax)
        norm = 1;

    // 计算输出网格的维度
    output_n_columns = (int)ceil((double)xdim / (double)jdec);
    output_n_rows = (int)ceil((double)ydim / (double)idec);
    
    // 创建输出网格 - 修复：正确创建输出网格
    wesn[0] = 0.0; wesn[1] = output_n_columns;
    wesn[2] = 0.0; wesn[3] = output_n_rows;
    inc[0] = 1.0; inc[1] = 1.0;
    
    Out = GMT_Create_Data(API, GMT_IS_GRID, GMT_IS_SURFACE, 
                         GMT_GRID_ALL, NULL, wesn, inc, 
                         GMT_GRID_PIXEL_REG, 0, NULL);
    if (Out == NULL) {
        die("Failed to create output grid", "");
    }
    
    // 设置输出网格维度
    Out->header->n_columns = output_n_columns;
    Out->header->n_rows = output_n_rows;

    if (format_flag == 1) {
        indat = (float *)malloc(4 * xdim);
        if (indat == NULL) die("memory allocation", "");
        read_float(indat, xdim, f_input, 0, buffer, ibuff);
    }

    ic0 = 0;
    iend = ylen = ibuff;

    for (ic = 0, row = 0; ic < ydim; ic += idec, row++) {
        // 确保行索引不越界
        if (row >= output_n_rows) break;
        
        left_node = GMT_Get_Index(API, Out->header, row, 0);
        ic1 = ic - ic0;

        // 计算当前行需要处理的列数
        int ncol = output_n_columns;
        
        // 并行化卷积处理
#ifdef _OPENMP
#pragma omp parallel for schedule(static) private(jc)
#endif
        for (int jidx = 0; jidx < ncol; jidx++) {
            jc = jidx * jdec;
            
            // 确保列索引不越界
            if (jc >= xdim) continue;
            
            float filtdat_loc = 0.0f;
            float rnorm_loc = 0.0f;
            
            // 调用卷积函数
            conv2d(buffer, &ylen, &xdim,
                   filter, &yarr, &xarr,
                   &filtdat_loc, &ic1, &jc, &rnorm_loc);
            
            float outv = 0.0f;
            
            if (norm > 0) {
                if (fabs(rnorm_loc) > (0.01 * rnormax) && fabs(rnorm_loc) > 1e-10)
                    outv = filtdat_loc / rnorm_loc;
            } else {
                if (fabs(rnorm_loc) >= 0.0001 * anormax)
                    outv = filtdat_loc;
            }
            
            // 安全写入输出
            if (jidx < output_n_columns) {
                Out->data[left_node + jidx] = outv;
            }
        }

        // 更新缓冲区（根据需要读取新行）
        // 注意：这里需要根据原来的逻辑更新缓冲区
        // 由于OpenMP并行化，我们需要确保缓冲区更新在并行区域之外
        if (ic + idec < ydim && ic + ibuff - yarr2 > ylen) {
            // 移动缓冲区内容
            int move_lines = ibuff - yarr2;
            if (move_lines > 0) {
                for (i = 0; i < move_lines; i++) {
                    for (j = 0; j < xdim; j++) {
                        buffer[j + xdim * i] = buffer[j + xdim * (i + yarr2)];
                    }
                }
            }
            
            // 读取新数据
            int new_lines = min(ibuff, ydim - ic - ic0);
            if (new_lines > 0) {
                if (format_flag == 1) {
                    read_float(indat, xdim, f_input, move_lines, buffer, new_lines);
                }
                // 其他格式的读取...
            }
            ylen = move_lines + new_lines;
        }
    }

    // 写入输出文件
    if (GMT_Write_Data(API, GMT_IS_GRID, GMT_IS_FILE, GMT_IS_SURFACE,
                       GMT_GRID_ALL, NULL, output_name, Out) != GMT_NOERROR) {
        fprintf(stderr, "Error writing output file\n");
    }

    // 清理内存
    if (filter != NULL) free(filter);
    if (buffer != NULL) free(buffer);
    if (indat != NULL) free(indat);
    if (cindat != NULL) free(cindat);
    if (cfdat != NULL) free(cfdat);
    
    GMT_Destroy_Data(API, &Out);
    if (In != NULL) GMT_Destroy_Data(API, &In);
    GMT_Destroy_Session(API);
    
    return (EXIT_SUCCESS);
}