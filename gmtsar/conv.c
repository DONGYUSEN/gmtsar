/*  conv_omp.c : OpenMP-accelerated version of conv.c (GMTSAR) */
/*  Based on original conv.c by D. Sandwell, with safe OpenMP parallelization */

#include "gmtsar.h"
#include <omp.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define min(x, y) (((x) < (y)) ? (x) : (y))
#define max(x, y) (((x) > (y)) ? (x) : (y))

char *USAGE = "conv [GMTSAR] - 2-D image convolution (OpenMP)\n\n"
              "Usage: conv idec jdec filter_file input output \n"
              "   idec           - row decimation factor \n"
              "   jdec           - column decimation factor \n"
              "   filter_file    - eg. filters/gauss17x5 \n"
              "   input          - name of file to be filtered (I*2 or R*4) \n"
              "   output         - name of filtered output file (R*4 only) \n\n";

int input_file_type, format_flag;

/*-------------------------------------------------------------*/
int determine_file_type(char *name, int *input_file_type) {
    int n, m;
    char tail[8];

    *input_file_type = 1;

    n = (int)strlen(name);
    m = n - 3;
    strncpy(&tail[0], &name[m], 4);
    if (verbose)
        fprintf(stderr, " name %s tail %s \n", name, tail);

    if ((strncmp(tail, "PRM", 3) == 0) || (strncmp(tail, "prm", 3) == 0)) {
        if (verbose)
            fprintf(stderr, " input: PRM file\n");
        *input_file_type = 2;
    }

    if (*input_file_type == 1)
        if (verbose)
            fprintf(stderr, " input: GMT binary\n");

    return (EXIT_SUCCESS);
}

/*-------------------------------------------------------------*/
FILE *read_PRM_file(char *prmfilename, char *input_file_name, struct PRM p, int *xdim, int *ydim) {
    FILE *f_input_prm, *f_input;

    if (verbose)
        fprintf(stderr, " reading PRM file %s\n", prmfilename);
    if ((f_input_prm = fopen(prmfilename, "r")) == NULL)
        die("Can't open input header", prmfilename);
    null_sio_struct(&p);
    get_sio_struct(f_input_prm, &p);
    strcpy(input_file_name, p.SLC_file);
    format_flag = 2;
    if (strncmp(p.dtype, "c", 1) == 0)
        format_flag = 3;
    if (verbose)
        fprintf(stderr, " reading PRM file %s\n", input_file_name);
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

/*-------------------------------------------------------------*/
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

/*-------------------------------------------------------------*/
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
/* Optimized single-point convolution kernel (no size check inside) */
static inline void conv2d2(float * restrict rdat, int ni, int nj,
                           float * restrict filt, int nif, int njf,
                           float *fdat, int ic, int jc, float *rnorm)
{
    int nif2 = nif / 2;
    int njf2 = njf / 2;

    int i0 = max(0, ic - nif2);
    int i1 = min(ni - 1, ic + nif2);
    int j0 = max(0, jc - njf2);
    int j1 = min(nj - 1, jc + njf2);

    float sum = 0.0f;
    float norm = 0.0f;

    for (int i = i0; i <= i1; i++) {
        int iflt = i - ic + nif2;
        float *pf = filt + iflt * njf;
        float *pr = rdat + i * nj;
        for (int j = j0; j <= j1; j++) {
            int jflt = j - jc + njf2;
            float w = pf[jflt];
            sum  += w * pr[j];
            norm += w;
        }
    }

    *fdat  = sum;
    *rnorm = norm;
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
    float filtin, filtdat, rnorm, rnormax, anormax;
    FILE *f_filter = NULL, *f_input = NULL;
    struct PRM p;
    void *API = NULL;
    struct GMT_GRID *Out = NULL;
    struct GMT_GRID *In = NULL;

    if (argc < 6)
        die("\n", USAGE);

    if ((API = GMT_Create_Session(argv[0], 0U, 0U, NULL)) == NULL)
        return EXIT_FAILURE;
		
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

    ibuff = 512;
    verbose = 0;

    null_sio_struct(&p);
    input_file_type = 1;
    format_flag = 1;

    idec = atoi(argv[1]);
    jdec = atoi(argv[2]);
    if (idec <= 0 || jdec <= 0)
        die("idec and jdec should be positive integers.", "");

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
        if ((c = strstr(input_name, "=bf"))) c[0] = '\0';
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

    if (fscanf(f_filter, "%d%d", &xarr, &yarr) != 2 ||
        xarr < 1 || yarr < 1 || (xarr & 1) == 0 || (yarr & 1) == 0)
        die("filter incomplete or not odd-sized", "");

    if (ibuff < yarr)
        die("increase dimension of ibuff", "");

    /* output size */
    iout = jout = 0;
    for (ic = 0; ic < ydim; ic += idec) iout++;
    for (jc = 0; jc < xdim; jc += jdec) jout++;

    inc[GMT_X] = round(xmax / (double)jout);
    inc[GMT_Y] = round(ymax / (double)iout);
    jout = floor(xmax / inc[GMT_X]);
    iout = floor(ymax / inc[GMT_Y]);

    wesn[GMT_XLO] = 0.0;
    wesn[GMT_XHI] = inc[GMT_X] * jout;
    wesn[GMT_YLO] = 0.0;
    wesn[GMT_YHI] = inc[GMT_Y] * iout;

    if ((Out = GMT_Create_Data(API, GMT_IS_GRID, GMT_IS_SURFACE, GMT_GRID_ALL,
                               NULL, wesn, inc, GMT_GRID_PIXEL_REG, 0, NULL)) == NULL)
        die("could not allocate output grid", "");

    narr = xarr * yarr;
    yarr2 = yarr / 2;
    if (ydim < ibuff) ibuff = ydim;
    imove = ibuff - yarr;
    nbuff = xdim * ibuff;

    filter = (float *)malloc(sizeof(float) * narr);
    buffer = (float *)malloc(sizeof(float) * nbuff);

    if (!filter || !buffer)
        die("memory allocation", "");

    if (format_flag == 1) indat = (float *)malloc(4 * xdim);
    if (format_flag == 2) cindat = (short *)malloc(4 * xdim);
    if (format_flag == 3) cfdat = (float *)malloc(8 * xdim);

    anormax = rnormax = 0.0f;
    for (i = 0; i < narr; i++) {
        if (fscanf(f_filter, "%f", &filtin) == EOF)
            die("filter incomplete", "");
        filter[i] = filtin;
        anormax += fabsf(filter[i]);
        rnormax += filter[i];
    }

    norm = 0;
    if (fabs(rnormax) > 0.05 * anormax)
        norm = 1;

    ic0 = 0;
    iend = ylen = ibuff;

    if (format_flag == 1) read_float(indat, xdim, f_input, 0, buffer, ibuff);
    if (format_flag == 2) read_SLC_int(cindat, xdim, f_input, 0, buffer, DFACT, ibuff);
    if (format_flag == 3) read_SLC_float(cfdat, xdim, f_input, 0, buffer, DFACT, ibuff);

    for (ic = 0, row = 0; ic < iout * idec; ic += idec) {

        if ((ic + yarr2) >= iend && (ic + yarr2) < (ydim - 1)) {

            for (i = 0; i < yarr; i++)
                for (j = 0; j < xdim; j++)
                    buffer[j + xdim * i] = buffer[j + xdim * (i + (ylen - yarr))];

            iread = MIN(imove, (ydim - iend));
            ic0 = iend - yarr;
            iend = iend + iread;
            ylen = iread + yarr;

            if (format_flag == 1) read_float(indat, xdim, f_input, yarr, buffer, iread);
            if (format_flag == 2) read_SLC_int(cindat, xdim, f_input, yarr, buffer, DFACT, iread);
            if (format_flag == 3) read_SLC_float(cfdat, xdim, f_input, yarr, buffer, DFACT, iread);
        }

        left_node = GMT_Get_Index(API, Out->header, row, 0);
        ic1 = ic - ic0;

        int njc_max = (int)floor(xmax / inc[GMT_X]) * jdec;
        int nj_out = 0;
        for (jc = 0; jc < njc_max; jc += jdec) nj_out++;

        /* =================== OpenMP PARALLEL REGION =================== */
        #pragma omp parallel for schedule(static)
        for (int jj = 0; jj < nj_out; jj++) {

            int jc_loc = jj * jdec;
            float filtdat_loc = 0.0f;
            float rnorm_loc = 0.0f;

            conv2d2(buffer, ylen, xdim, filter, yarr, xarr,
                    &filtdat_loc, ic1, jc_loc, &rnorm_loc);

            float outv = 0.0f;
            if (norm > 0) {
                if (fabsf(rnorm_loc) > (0.01f * rnormax))
                    outv = filtdat_loc / rnorm_loc;
            }
            else {
                if (fabsf(rnorm_loc) < 0.0001f * anormax)
                    outv = filtdat_loc;
            }

            Out->data[left_node + jj] = outv;
        }
        /* =============================================================== */

        row++;
    }

    fclose(f_input);

    if (GMT_Write_Data(API, GMT_IS_GRID, GMT_IS_FILE, GMT_IS_SURFACE,
                       GMT_GRID_ALL, NULL, output_name, Out)) {
        die("Failed to write output grid", "");
    }

    if (GMT_Destroy_Session(API))
        return EXIT_FAILURE;

    return (EXIT_SUCCESS);
}