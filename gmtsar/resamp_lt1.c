/*=============================================================
  OpenMP-safe resamp.c (GMTSAR)
  Parallelized over azimuth lines (ii)
=============================================================*/

#include "gmtsar.h"
#include "lib_functions.h"

#include <omp.h>
#include <fcntl.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/mman.h>
#include <sys/types.h>
#include <unistd.h>


#define MAX_COEFFS 7  // 每个方向最多6个系数


/*-------------------------------------------------------------*/
char *USAGE =
    "\nUsage: "
    "resamp_lt1  master.PRM aligned.PRM new_aligned.PRM new_aligned.SLC intrp fitfile\n"
    "            intrp: 1-nearest  2-bilinear  3-bicubic  4-bisinc\n";

/*---------------- function prototypes ------------------------*/
void ram2ras(struct PRM ps, double *ram, double *ras);
void nearest(double *, short *, int, int, short *);
void bilinear(double *, short *, int, int, short *);
void bicubic(double *, short *, int, int, short *);
void bisinc(double *, short *, int, int, short *);
void print_prm_params(struct PRM, struct PRM);
void fix_prm_params(struct PRM *, char *);
void bicubic_one(double *, double *, double, double, double *);
double cubic_kernel(double, double);
void sinc_one(double *, double *, double, double, double *);
double sinc_kernel(double);

	// add by dong

double a_coeffs[MAX_COEFFS] = {0};
double b_coeffs[MAX_COEFFS] = {0};
double error_x = 0.0;
double error_y = 0.0;


// 从文件中读取拟合系数和误差
int read_fit_coefficients(const char *filename, 
                         double *a_coeffs, double *b_coeffs,
                         double *error_x, double *error_y) {
    FILE *fp = fopen(filename, "r");
    if (!fp) {
        printf("错误: 无法打开文件 %s\n", filename);
        return -1;
    }
    
    char line[1024];
    if (fgets(line, sizeof(line), fp) == NULL) {
        printf("错误: 文件为空\n");
        fclose(fp);
        return -1;
    }
    
    fclose(fp);
    
    // 解析行中的数值
    char *token = strtok(line, " \t\n");
    int count = 0;
    double values[14];  // 12个系数 + 2个误差 = 14个值
    
    while (token != NULL && count < 14) {
        values[count] = atof(token);
        count++;
        token = strtok(NULL, " \t\n");
    }
    
    if (count < 14) {
        printf("警告: 只读取到 %d 个值，预期 14 个\n", count);
    }
    
    // 将值复制到相应的数组中
    for (int i = 0; i < MAX_COEFFS; i++) {
        if (i < count) {
            a_coeffs[i] = values[i];
        } else {
            a_coeffs[i] = 0.0;
        }
    }
    
    for (int i = 0; i < MAX_COEFFS; i++) {
        int idx = MAX_COEFFS + i;
        if (idx < count) {
            b_coeffs[i] = values[idx];
        } else {
            b_coeffs[i] = 0.0;
        }
    }
    
    // 读取误差值
    if (count >= 13) {
        *error_x = values[12];
    } else {
        *error_x = 0.0;
    }
    
    if (count >= 14) {
        *error_y = values[13];
    } else {
        *error_y = 0.0;
    }
    
    return 0;
}



/*=============================================================*/
int main(int argc, char **argv) {

    int ii, jj;
    int intrp;
    int xdimm, ydimm;
    int xdims, ydims;
    short *sinn = NULL;
    FILE *SLC_file2 = NULL, *prmout = NULL;
    int fdin;
    size_t st_size;
    double ram[2], ras[2];


    struct PRM pm, ps;
    double sv_pr[6];

    if (argc < 7)
        die(USAGE, "");


    if (read_fit_coefficients(argv[6], a_coeffs, b_coeffs, &error_x, &error_y) == 0) {
        printf("成功读取系数:\n");
        
        printf("X方向系数 (a0-a5):\n");
        for (int i = 0; i < MAX_COEFFS; i++) {
            printf("  a%d = %.12e\n", i, a_coeffs[i]);
        }
        
        printf("\nY方向系数 (b0-b5):\n");
        for (int i = 0; i < MAX_COEFFS; i++) {
            printf("  b%d = %.12e\n", i, b_coeffs[i]);
        }
        
        //printf("\n误差:\n");
        //printf("  X方向误差: %.12e\n", error_x);
        //printf("  Y方向误差: %.12e\n", error_y);
    } else {
        printf("读取系数失败\n");
    }


    intrp = atoi(argv[5]);

    get_prm(&pm, argv[1]);
    get_prm(&ps, argv[2]);

    xdimm = pm.num_rng_bins;
    ydimm = pm.num_patches * pm.num_valid_az;
    xdims = ps.num_rng_bins;
    ydims = ps.num_patches * ps.num_valid_az;

    /* nearest neighbor: force integer-only */
    if (intrp == 1) {
        sv_pr[0] = ps.sub_int_r; ps.sub_int_r = 0.;
        sv_pr[1] = ps.stretch_r; ps.stretch_r = 0.;
        sv_pr[2] = ps.a_stretch_r; ps.a_stretch_r = 0.;
        sv_pr[3] = ps.sub_int_a; ps.sub_int_a = 0.;
        sv_pr[4] = ps.stretch_a; ps.stretch_a = 0.;
        sv_pr[5] = ps.a_stretch_a; ps.a_stretch_a = 0.;
    }

	//printf("a:%f, %f, %f, %f", ps.ashift,ps.sub_int_a,ps.stretch_a, ps.a_stretch_a);

    /* mmap input SLC */
    if ((fdin = open(ps.SLC_file, O_RDONLY)) < 0)
        die("cannot open input SLC", ps.SLC_file);

    st_size = (size_t)4 * (size_t)xdims * (size_t)ydims;

    sinn = mmap(NULL, st_size, PROT_READ, MAP_SHARED, fdin, 0);
    if (sinn == MAP_FAILED)
        die("mmap failed", "");

    /* open output */
    if ((SLC_file2 = fopen(argv[4], "wb")) == NULL)
        die("cannot open output SLC", argv[4]);

		    /* open output */
    //if ((fitfile = fopen(argv[7], "r")) == NULL)
    //    die("cannot open output fitfile", argv[7]);

    /*=========================================================
      OpenMP parallel region
    =========================================================*/
    //#pragma omp parallel
    {
        short *sout_local = (short *)malloc(2 * xdimm * sizeof(short));
        if (!sout_local) {
            fprintf(stderr, "OMP malloc failed\n");
            exit(1);
        }

        #pragma omp for schedule(static)
        for (ii = 0; ii < ydimm; ii++) {

            for (jj = 0; jj < xdimm; jj++) {

                ram[0] = (double)jj;
                ram[1] = (double)ii;
                ram2ras(ps, ram, ras);

                if (intrp == 1)
                    nearest(ras, sinn, ydims, xdims, &sout_local[2 * jj]);
                else if (intrp == 2)
                    bilinear(ras, sinn, ydims, xdims, &sout_local[2 * jj]);
                else if (intrp == 3)
                    bicubic(ras, sinn, ydims, xdims, &sout_local[2 * jj]);
                else
                    bisinc(ras, sinn, ydims, xdims, &sout_local[2 * jj]);
            }

            #pragma omp critical
            {
                fseek(SLC_file2,
                      (long long)ii * 2 * xdimm * sizeof(short),
                      SEEK_SET);
                fwrite(sout_local, 2 * sizeof(short), xdimm, SLC_file2);
            }
        }

        free(sout_local);
    }

    /* restore parameters */
    if (intrp == 1) {
        ps.sub_int_r = sv_pr[0];
        ps.stretch_r = sv_pr[1];
        ps.a_stretch_r = sv_pr[2];
        ps.sub_int_a = sv_pr[3];
        ps.stretch_a = sv_pr[4];
        ps.a_stretch_a = sv_pr[5];
    }


    /* update PRM */
    ps.num_rng_bins = pm.num_rng_bins;
    ps.fs = pm.fs;
    ps.bytes_per_line = pm.bytes_per_line;
    ps.good_bytes = pm.good_bytes;
    ps.prf = pm.prf;
    ps.num_valid_az = pm.num_valid_az;
    ps.num_lines = pm.num_lines;
    ps.num_patches = pm.num_patches;
    ps.nrows = pm.nrows;
    //printf(" Ready update PRM \n");

    if ((prmout = fopen(argv[3], "w")) == NULL)
        die("cannot open output PRM", argv[3]);

    put_sio_struct(ps, prmout);

    fclose(prmout);
    fclose(SLC_file2);
    munmap(sinn, st_size);
    close(fdin);

    return EXIT_SUCCESS;
}



/*=============================================================
  Below: interpolation & math routines
  （与原 GMTSAR 完全一致，未做任何修改）
=============================================================*/



/************************************************************************
 * bi-cubic interpolation algorithm modified from GMT                    *
 ************************************************************************/
/************************************************************************
 * Creator: David Sandwell       (Scripps Institution of Oceanography)   *
 * Date   : 03/22/13                                                     *
 ************************************************************************/
#include <math.h>
#include <stdio.h>

double cubic_kernel(double, double);

void bicubic_one(double *rdata, double *idata, double x, double y, double *cz) {
	int i, j, ij;
	double wx[4], wy[4];
	double arg, w, wsum, rsum, isum;
	double a = -0.3;

	/* These weights are based on the cubic convolution kernel, see for example
	   http://undergraduate.csse.uwa.edu.au/units/CITS4241/Handouts/Lecture04.html
	   These weights include a free parameter (a).
	*/

	for (i = 0; i < 4; i++) {
		arg = fabs(x + 1 - i);
		wx[i] = cubic_kernel(arg, a);
		arg = fabs(y + 1 - i);
		wy[i] = cubic_kernel(arg, a);
	}

	rsum = isum = wsum = 0.0;
	ij = 0;
	for (j = 0; j < 4; j++) {
		for (i = 0; i < 4; i++) {
			w = wx[i] * wy[j];
			rsum += rdata[ij + i] * w;
			isum += idata[ij + i] * w;
			wsum += w;
		}
		ij += 4;
	}
	if (wsum <= 0.0)
		printf(" error wsum is zero \n");
	cz[0] = rsum / wsum;
	cz[1] = isum / wsum;
}
/************************************************************************
  kernel computes a bi-cubic spline kernel using the formula given at
  the following web page
  http://undergraduate.csse.uwa.edu.au/units/CITS4241/Handouts/Lecture04.html
************************************************************************/
/************************************************************************
 * Creator: David Sandwell       (Scripps Institution of Oceanography)   *
 * Date   : 03/22/13                                                     *
 ************************************************************************/
#include <stdio.h>

double cubic_kernel(double arg, double a)

/* note arg must be positive and a must be between -3 and 0. */
{

	double arg2, arg3, f;

	arg2 = arg * arg;
	arg3 = arg2 * arg;
	if (arg <= 1.) {
		f = (a + 2) * arg3 - (a + 3) * arg2 + 1.;
	}
	else if (arg <= 2.) {
		f = a * arg3 - 5 * a * arg2 + 8 * a * arg - 4 * a;
	}
	else {
		f = 0.;
	}
	return (f);
}

/************************************************************************
 * nearest, bilinear, and bicubic interpolations                         *
 ************************************************************************/
/************************************************************************
 * Creator: David Sandwell       (Scripps Institution of Oceanography)   *
 * Date   : 03/22/13                                                     *
 ************************************************************************/
#include "gmtsar.h"
#include <math.h>
#include <stdio.h>

void nearest(double *ras, short *s_in, int ydims, int xdims, short *sout) {
	int i, j, k;
	short *tmp_sin;

	/* compute the indices of the upper left corner */

	j = (int)(ras[0] + 0.5);
	i = (int)(ras[1] + 0.5);
	// k = 2*xdims*i + 2*j;
	k = 2 * j;

	/* shift the pointer by i0-ns2 lines  */
	tmp_sin = s_in;
	tmp_sin++;
	tmp_sin = s_in + (size_t)(2 * xdims) * (size_t)i * (tmp_sin - s_in);

	/* use the nearest point if it is within the bounds of the aligned array */

	if (i < 0 || i >= ydims || j < 0 || j >= xdims) {
		sout[0] = 0;
		sout[1] = 0;
	}
	else {
		// sout[0] = s_in[k];
		// sout[1] = s_in[k+1];
		sout[0] = tmp_sin[k];
		sout[1] = tmp_sin[k + 1];
	}
}

void bilinear(double *ras, short *s_in, int ydims, int xdims, short *sout) {
	double dr, da, real, imag;
	int k00, k01, k10, k11;
	int i0, j0;
	int nclip;
	short *tmp_sin;

	/* compute the residual offsets */
	nclip = 0;
	j0 = (int)floor(ras[0]);
	i0 = (int)floor(ras[1]);
	dr = ras[0] - (double)j0;
	da = ras[1] - (double)i0;
	if (dr < 0. || dr > 1. || da < 0. || da > 1)
		fprintf(stderr, " dr or da out of bounds %f %f \n", dr, da);

	/* shift the pointer by i0-ns2 lines  */
	tmp_sin = s_in;
	tmp_sin++;
	tmp_sin = s_in + (size_t)(2 * xdims) * (size_t)i0 * (tmp_sin - s_in);

	/* compute the indices of the 4 corners */
	/*
	k00 = 2*xdims*i0     + 2*j0;
	k01 = 2*xdims*i0     + 2*(j0+1);
	k10 = 2*xdims*(i0+1) + 2*j0;
	k11 = 2*xdims*(i0+1) + 2*(j0+1);
	*/
	k00 = 2 * j0;
	k01 = 2 * (j0 + 1);
	k10 = 2 * xdims + 2 * j0;
	k11 = 2 * xdims + 2 * (j0 + 1);

	/* do the interpolation if all 4 corners are within the bounds of the aligned
	 * array */

	if (i0 < 0 || i0 >= (ydims - 1) || j0 < 0 || j0 >= (xdims - 1)) {
		sout[0] = 0;
		sout[1] = 0;
	}
	else {
		/*real = s_in[k00] * (1.0 - da) * (1.0 - dr)
		     + s_in[k10] * (da)       * (1.0 - dr)
		     + s_in[k01] * (1.0 - da) * (dr)
		     + s_in[k11] * (da)       * (dr);*/
		real = tmp_sin[k00] * (1.0 - da) * (1.0 - dr) + tmp_sin[k10] * (da) * (1.0 - dr) + tmp_sin[k01] * (1.0 - da) * (dr) +
		       tmp_sin[k11] * (da) * (dr);

		if ((int)fabs(real) > I2MAX)
			nclip = nclip + 1;
		sout[0] = (short)clipi2(real + 0.5);
		/*imag = s_in[k00+1] * (1.0 - da) * (1.0 - dr)
		     + s_in[k10+1] * (da)       * (1.0 - dr)
		     + s_in[k01+1] * (1.0 - da) * (dr)
		     + s_in[k11+1] * (da)       * (dr);*/
		imag = tmp_sin[k00 + 1] * (1.0 - da) * (1.0 - dr) + tmp_sin[k10 + 1] * (da) * (1.0 - dr) +
		       tmp_sin[k01 + 1] * (1.0 - da) * (dr) + tmp_sin[k11 + 1] * (da) * (dr);

		if ((int)fabs(imag) > I2MAX)
			nclip = nclip + 1;
		sout[1] = (short)clipi2(imag + 0.5);
	}
	/*if(nclip > 0) fprintf(stderr," %d integers were clipped \n",nclip);*/
}

void bicubic(double *ras, short *s_in, int ydims, int xdims, short *sout) {
	double dr, da;
	double rdata[16], idata[16], cz[2];
	int i, j, k, kk;
	int i0, j0;
	int nclip;
	short *tmp_sin;

	/* compute the residual offsets */
	nclip = 0;
	j0 = (int)floor(ras[0]);
	i0 = (int)floor(ras[1]);
	dr = ras[0] - (double)j0;
	da = ras[1] - (double)i0;
	if (dr < 0. || dr > 1. || da < 0. || da > 1)
		fprintf(stderr, " dr or da out of bounds %f %f \n", dr, da);

	/* make sure all 4 corners are within the bounds of the aligned array */

	if ((i0 - 1) < 0 || (i0 + 2) >= ydims || (j0 - 1) < 0 || (j0 + 2) >= xdims) {
		sout[0] = 0;
		sout[1] = 0;
	}
	else {

		/* safe to do the interpolation */

		/* shift the pointer by i0-1 lines  */
		tmp_sin = s_in;
		tmp_sin++;
		tmp_sin = s_in + (size_t)(2 * xdims) * (size_t)(i0 - 1) * (tmp_sin - s_in);

		for (i = 0; i < 4; i++) {
			for (j = 0; j < 4; j++) {
				k = i * 4 + j;
				// kk = 2*xdims*(i0-1+i)  + 2*(j0-1+j);
				kk = 2 * xdims * i + 2 * (j0 - 1 + j);
				// rdata[k] = s_in[kk];
				// idata[k] = s_in[kk+1];
				rdata[k] = tmp_sin[kk];
				idata[k] = tmp_sin[kk + 1];
			}
		}

		/* interpolate the real and imaginary data */

		bicubic_one(rdata, idata, dr, da, cz);

		if ((int)fabs(cz[0]) > I2MAX)
			nclip = nclip + 1;
		sout[0] = (short)clipi2(cz[0] + 0.5);
		if ((int)fabs(cz[1]) > I2MAX)
			nclip = nclip + 1;
		sout[1] = (short)clipi2(cz[1] + 0.5);
	}
	/*if(nclip > 0) fprintf(stderr," %d integers were clipped \n",nclip);*/
}

void bisinc(double *ras, short *s_in, int ydims, int xdims, short *sout) {
	double dr, da, ns2 = NS / 2 - 1;
	double rdata[NS * NS], idata[NS * NS], cz[2];
	int i, j, k, kk;
	int i0, j0;
	int nclip;
	short *tmp_sin;

	/* compute the residual offsets */
	nclip = 0;
	j0 = (int)floor(ras[0]);
	i0 = (int)floor(ras[1]);
	dr = ras[0] - (double)j0;
	da = ras[1] - (double)i0;
	if (dr < 0. || dr > 1. || da < 0. || da > 1)
		fprintf(stderr, " dr or da out of bounds %f %f \n", dr, da);

	/* make sure all 4 corners are within the bounds of the aligned array */

	if ((i0 - ns2) < 0 || (i0 + ns2 + 1) >= ydims || (j0 - ns2) < 0 || (j0 + ns2 + 1) >= xdims) {
		sout[0] = 0;
		sout[1] = 0;
	}
	else {

		/* safe to do the interpolation */

		/* shift the pointer by i0-ns2 lines  */
		tmp_sin = s_in;
		tmp_sin++;
		tmp_sin = s_in + (size_t)(2 * xdims) * (size_t)(i0 - ns2) * (tmp_sin - s_in);

		for (i = 0; i < NS; i++) {
			for (j = 0; j < NS; j++) {
				k = i * NS + j;
				// kk = 2*xdims*(i0-ns2+i)  + 2*(j0-ns2+j);
				kk = 2 * xdims * i + 2 * (j0 + j - ns2);
				rdata[k] = tmp_sin[kk];
				idata[k] = tmp_sin[kk + 1];
			}
		}

		/* interpolate the real and imaginary data */

		sinc_one(rdata, idata, dr, da, cz);

		if ((int)fabs(cz[0]) > I2MAX)
			nclip = nclip + 1;
		sout[0] = (short)clipi2(cz[0] + 0.5);
		if ((int)fabs(cz[1]) > I2MAX)
			nclip = nclip + 1;
		sout[1] = (short)clipi2(cz[1] + 0.5);
	}
	/*if(nclip > 0) fprintf(stderr," %d integers were clipped \n",nclip);*/
}
#include "gmtsar.h"
#include "lib_functions.h"

/*--------------------------------------------------------------*/
void print_prm_params(struct PRM p1, struct PRM p2) {
	fprintf(stderr, " SLC 1: num_rng_bins %d num_lines %d \n", p1.num_rng_bins, p1.num_lines);
	fprintf(stderr, " SLC 2: num_rng_bins %d num_lines %d \n", p2.num_rng_bins, p2.num_lines);
	fprintf(stderr, " lambda %f \n", p2.lambda);
	fprintf(stderr, " near_range %f \n", p2.near_range);
	fprintf(stderr, " rng_samp_rate %.7f \n", p2.fs);
	fprintf(stderr, " sc_clock_start %f \n", p2.SC_clock_start);
	fprintf(stderr, " sc_clock_stop %f \n", p2.SC_clock_stop);
	fprintf(stderr, " clock_start %f \n", p2.clock_start);
	fprintf(stderr, " clock_stop %f \n", p2.clock_stop);
	fprintf(stderr, " prfm %f \n", p1.prf);
	fprintf(stderr, " prfs %f \n", p2.prf);
	fprintf(stderr, " rshift %f \n", p2.rshift + p2.sub_int_r);
	fprintf(stderr, " ashift %f \n", p2.ashift + p2.sub_int_a);
}

/*--------------------------------------------------------------*/
void fix_prm_params(struct PRM *p, char *s) {
	double delr;

	delr = SOL / p->fs / 2.0;

	/* these are from prm2gips */
	p->near_range = p->near_range + (p->st_rng_bin - p->chirp_ext + p->rshift - 1) * delr;
	p->SC_clock_start =
	    p->SC_clock_start + p->ashift / (p->prf * 86400.0) + (p->nrows - p->num_valid_az) / (2.0 * p->prf * 86400);
	p->SC_clock_stop = p->SC_clock_start + (p->num_valid_az * p->num_patches) / (p->prf * 86400.0);
}
/************************************************************************
 * ram2ras maps range and azimuth from a master image location into the  *
 * corresponding range and azimuth location of the aligned image.          *
 ************************************************************************/
/************************************************************************
 * Creator: David Sandwell       (Scripps Institution of Oceanography)   *
 * Date   : 03/22/13                                                     *
 ************************************************************************/
/************************************************************************
 * Modification history:                                                 *
 *                                                                       *
 * DATE                                                                  *
 ************************************************************************/
#include "gmtsar.h"
#include <stdio.h>
#include <stdlib.h>

void ram2ras(struct PRM ps, double *ram, double *ras) {
	/* this is the range coordinate */
	//ras[0] = ram[0] + ((ps.rshift + ps.sub_int_r) + ram[0] * ps.stretch_r + ram[1] * ps.a_stretch_r);
	ras[0] = ram[0] + a_coeffs[0] + ram[0] * a_coeffs[1] +  pow(ram[0],2) * a_coeffs[2] + pow(ram[0],3) *  a_coeffs[3] + pow(ram[0],4) * a_coeffs[4] + pow(ram[0],5) *  a_coeffs[5];
	//ras[0] = ram[0] + a_coeffs[0] + ram[0] * a_coeffs[1] +  ram[1] * a_coeffs[2] + ram[1] * ram[0] * a_coeffs[3] + ram[0] * ram[0] * a_coeffs[4] + ram[1] * ram[1] * a_coeffs[5];

	/* this is the azimuth coordinate */
	//ras[1] = ram[1] + ((ps.ashift + ps.sub_int_a) + ram[0] * ps.stretch_a + ram[1] * ps.a_stretch_a);
	ras[1] = ram[1] + b_coeffs[0] + ram[1] * b_coeffs[1] +  pow(ram[1],2) * b_coeffs[2] + pow(ram[1],3) * b_coeffs[3] + pow(ram[1],4) *  b_coeffs[4] + pow(ram[1],5) * b_coeffs[5];
	//ras[1] = ram[1] + b_coeffs[0] + ram[1] * b_coeffs[1] +  ram[0] * b_coeffs[2] + ram[1] * ram[0] * b_coeffs[3] + ram[1] * ram[1] * b_coeffs[4] + ram[0] * ram[0] * b_coeffs[5];
}


/************************************************************************
  computes sinc function kernel for interpolation
************************************************************************/
/************************************************************************
 * Creator: David Sandwell       (Scripps Institution of Oceanography)   *
 * Date   : 03/28/13                                                     *
 ************************************************************************/
#include <math.h>
#define PI 3.1415926535897932

double sinc_kernel(double x) {
	double arg, f;

	arg = fabs(PI * x);
	if (arg > 0.) {
		f = sin(arg) / arg;
	}
	else {
		f = 1.;
	}
	return (f);
}
/************************************************************************
 * sinc function interpolation                                           *
 ************************************************************************/
/************************************************************************
 * Creator: David Sandwell       (Scripps Institution of Oceanography)   *
 * Date   : 03/28/13                                                     *
 ************************************************************************/
#include "gmtsar.h"
#include <math.h>
#include <stdio.h>

double sinc_kernel(double);

void sinc_one(double *rdata, double *idata, double x, double y, double *cz) {
	int i, j, ij, ns2 = NS / 2 - 1;
	double wx[NS], wy[NS];
	double arg, w, wsum, rsum, isum;

	for (i = 0; i < NS; i++) {
		arg = fabs(x + ns2 - i);
		wx[i] = sinc_kernel(arg);
		arg = fabs(y + ns2 - i);
		wy[i] = sinc_kernel(arg);
	}

	rsum = isum = wsum = 0.0;
	ij = 0;
	for (j = 0; j < NS; j++) {
		for (i = 0; i < NS; i++) {
			w = wx[i] * wy[j];
			rsum += rdata[ij + i] * w;
			isum += idata[ij + i] * w;
			wsum += w;
		}
		ij += NS;
	}
	if (wsum <= 0.0)
		printf(" error wsum is zero \n");
	cz[0] = rsum / wsum;
	cz[1] = isum / wsum;
}

