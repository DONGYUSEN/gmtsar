/* SAT_llt2rat2_omp.c */
/* OpenMP parallel version - optimized for performance and memory safety */

#include "gmtsar.h"
#include "llt2xyz.h"
#include "orbit.h"

#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>

#define R 0.61803399
#define C (1.0 - R)
#define TOL 2
#define SOL 299792458.0
#define INIT_CAPACITY 1000000

/* 优化轨道缓存结构 */
typedef struct {
    double time;
    double x, y, z;
} OrbPos;

/* 预计算的结构 */
typedef struct {
    double inv_dr;
    double prf;
    double t1;
    double near_range;
    double dr;
    double ra;
    double fll;
    double RE;
} CalcParams;

/* 外部函数声明 */
void read_orb(FILE *, struct SAT_ORB *);
void set_prm_defaults(struct PRM *);
void hermite_c(double *, double *, double *, int, int, double, double *, int *);
void null_sio_struct(struct PRM *);
void get_sio_struct(FILE *, struct PRM *);
void plh2xyz(double *, double *, double, double);

/* 快速距离计算 */
static inline double dist3_fast(double x, double y, double z, const OrbPos *o) {
    double dx = x - o->x;
    double dy = y - o->y;
    double dz = z - o->z;
    return sqrt(dx*dx + dy*dy + dz*dz);
}

/* 优化的黄金搜索算法 - 修复边界检查 */
static int golden_search(const OrbPos *orb, int a, int b, 
                         double x, double y, double z, 
                         double *rng, double *tm) {
    // 确保索引有效
    if (a >= b) {
        *rng = dist3_fast(x, y, z, &orb[a]);
        *tm = orb[a].time;
        return a;
    }
    
    int x0 = a, x3 = b;
    // 确保x1, x2在有效范围内
    int x1 = a + (int)((b - a) * C);
    int x2 = b - (int)((b - a) * C);
    
    // 确保索引不越界
    if (x1 < a) x1 = a;
    if (x1 > b) x1 = b;
    if (x2 < a) x2 = a;
    if (x2 > b) x2 = b;
    
    double f1 = dist3_fast(x, y, z, &orb[x1]);
    double f2 = dist3_fast(x, y, z, &orb[x2]);
    
    int iterations = 0;
    const int max_iter = 100;
    
    while ((x3 - x0) > TOL && iterations++ < max_iter) {
        if (f2 < f1) {
            x0 = x1; 
            x1 = x2;
            x2 = (int)(R * x3 + C * x1);
            if (x2 < x1) x2 = x1 + 1;
            if (x2 > x3) x2 = x3;
            f1 = f2;
            f2 = dist3_fast(x, y, z, &orb[x2]);
        } else {
            x3 = x2; 
            x2 = x1;
            x1 = (int)(R * x0 + C * x2);
            if (x1 < x0) x1 = x0;
            if (x1 > x2) x1 = x2;
            f2 = f1;
            f1 = dist3_fast(x, y, z, &orb[x1]);
        }
    }
    
    // 最终确定最小值点
    int xmin;
    if (f1 < f2) {
        xmin = x1; 
        *rng = f1;
    } else {
        xmin = x2; 
        *rng = f2;
    }
    
    // 确保索引有效
    if (xmin < a) xmin = a;
    if (xmin > b) xmin = b;
    
    *tm = orb[xmin].time;
    return xmin;
}

/* 优化的轨道计算函数 - 修复内存管理 */
static int calculate_orbit_positions(struct SAT_ORB *orb, OrbPos *op, 
                                     double ts, double t1, int nrec, int npad) {
    int i, k;
    int nval = 6;
    int N = nrec + 2 * npad;
    
    // 检查轨道数据是否足够
    if (orb->nd < nval) {
        fprintf(stderr, "Error: Not enough orbit data points (%d < %d)\n", orb->nd, nval);
        return -1;
    }
    
    // 一次性分配所有临时数组
    double *pt = malloc(sizeof(double) * orb->nd);
    double *px = malloc(sizeof(double) * orb->nd);
    double *py = malloc(sizeof(double) * orb->nd);
    double *pz = malloc(sizeof(double) * orb->nd);
    double *pvx = malloc(sizeof(double) * orb->nd);
    double *pvy = malloc(sizeof(double) * orb->nd);
    double *pvz = malloc(sizeof(double) * orb->nd);
    
    if (!pt || !px || !py || !pz || !pvx || !pvy || !pvz) {
        if (pt) free(pt);
        if (px) free(px);
        if (py) free(py);
        if (pz) free(pz);
        if (pvx) free(pvx);
        if (pvy) free(pvy);
        if (pvz) free(pvz);
        return -1;
    }
    
    // 预计算时间基准
    double pt0 = 86400.0 * orb->id + orb->sec;
    
    // 注意：hermite_c可能不是线程安全的，改为串行
    for (k = 0; k < orb->nd; k++) {
        pt[k] = pt0 + k * orb->dsec;
        px[k] = orb->points[k].px;
        py[k] = orb->points[k].py;
        pz[k] = orb->points[k].pz;
        pvx[k] = orb->points[k].vx;
        pvy[k] = orb->points[k].vy;
        pvz[k] = orb->points[k].vz;
    }
    
    // 串行计算轨道位置（hermite_c可能不是线程安全的）
    for (i = 0; i < N; i++) {
        double time = t1 - npad * ts + i * ts;
        double xs, ys, zs;
        int ir = 0;
        
        // 插值计算位置
        hermite_c(pt, px, pvx, orb->nd, nval, time, &xs, &ir);
        hermite_c(pt, py, pvy, orb->nd, nval, time, &ys, &ir);
        hermite_c(pt, pz, pvz, orb->nd, nval, time, &zs, &ir);
        
        // 存储结果
        op[i].time = time;
        op[i].x = xs;
        op[i].y = ys;
        op[i].z = zs;
    }
    
    // 清理临时内存
    free(pt);
    free(px);
    free(py);
    free(pz);
    free(pvx);
    free(pvy);
    free(pvz);
    
    return N;
}

/* 优化的主处理函数 */
static void process_points(double *in, double *out, int npt,
                          const OrbPos *orb_pos, int orb_count,
                          const CalcParams *params) {
    int i;
    
    #pragma omp parallel for schedule(dynamic, 1024) private(i)
    for (i = 0; i < npt; i++) {
        double rp[3], xp[3];
        
        // 输入数据顺序：经度、纬度、高度
        rp[1] = in[3 * i + 0];     // 经度
        rp[0] = in[3 * i + 1];     // 纬度
        rp[2] = in[3 * i + 2];     // 高度
        
        // 将地理坐标转换为地心直角坐标
        plh2xyz(rp, xp, params->ra, params->fll);

        /* compute the topography due to the difference between the local radius and
         * center radius  修正高度 */
        rp[2] = sqrt(xp[0] * xp[0] + xp[1] * xp[1] + xp[2] * xp[2]) - params->RE;
        in[3 * i + 2] = rp[2];

        // 确保orb_count至少为1
        if (orb_count <= 0) {
            out[2 * i + 0] = 0.0;
            out[2 * i + 1] = 0.0;
            continue;
        }
        
        // 使用黄金搜索找到最近的点
        double rng0, tm;
        golden_search(orb_pos, 0, orb_count - 1, 
                     xp[0], xp[1], xp[2], &rng0, &tm);
        
        // 计算距离门和方位时间
        out[2 * i + 0] = (rng0 - params->near_range) * params->inv_dr;
        out[2 * i + 1] = params->prf * (tm - params->t1);
    }
}

/* 安全的内存重新分配 */
static void* safe_realloc(void *ptr, size_t new_size) {
    if (new_size == 0) {
        free(ptr);
        return NULL;
    }
    
    void *new_ptr = realloc(ptr, new_size);
    if (!new_ptr) {
        fprintf(stderr, "Error: Memory reallocation failed (requested %zu bytes)\n", new_size);
        free(ptr);
        exit(EXIT_FAILURE);
    }
    return new_ptr;
}

/* 主函数 */
int main(int argc, char **argv) {
    if (argc < 3) {
        fprintf(stderr, "Usage: %s prm_file precise [-bos|-bod]\n", argv[0]);
        fprintf(stderr, "  -bos: binary output (short integers)\n");
        fprintf(stderr, "  -bod: binary output (double precision)\n");
        return EXIT_FAILURE;
    }
    
    int binary_short = 0;
    int binary_double = 0;
    
    if (argc >= 4) {
        if (strcmp(argv[3], "-bos") == 0) binary_short = 1;
        else if (strcmp(argv[3], "-bod") == 0) binary_double = 1;
    }
    
    // 读取PRM参数
    struct PRM prm;
    null_sio_struct(&prm);
    set_prm_defaults(&prm);
    
    FILE *f = fopen(argv[1], "r");
    if (!f) {
        perror("Error opening PRM file");
        return EXIT_FAILURE;
    }
    get_sio_struct(f, &prm);
    fclose(f);
    
    int precise = atoi(argv[2]);
    if (precise <= 0) {
        fprintf(stderr, "Error: precise value must be positive\n");
        return EXIT_FAILURE;
    }
    
    // 读取轨道数据
    FILE *fo = fopen(prm.led_file, "r");
    if (!fo) {
        perror("Error opening LED file");
        return EXIT_FAILURE;
    }
    
    struct SAT_ORB *orb = malloc(sizeof(struct SAT_ORB));
    if (!orb) {
        perror("Memory allocation failed for orbit");
        fclose(fo);
        return EXIT_FAILURE;
    }
    
    read_orb(fo, orb);
    fclose(fo);
    
    // 检查轨道数据是否足够
    if (orb->nd < 6) {
        fprintf(stderr, "Error: insufficient orbit data points (%d)\n", orb->nd);
        free(orb->points);
        free(orb);
        return EXIT_FAILURE;
    }
    
    // 计算时间参数
    double t1 = 86400.0 * prm.clock_start;
    double t2 = t1 + prm.num_patches * prm.num_valid_az / prm.prf;
    double ts = (prm.prf < 600.0) ? (1.0 / prm.prf) : (2.0 / prm.prf);
    int npad = (prm.prf < 600.0) ? 20000 : 8000;
    
    int nrec = (int)((t2 - t1) / ts + 0.5);
    int total_orb_pos = nrec + 2 * npad;
    
    // 分配轨道位置缓存
    OrbPos *orb_pos = malloc(sizeof(OrbPos) * total_orb_pos);
    if (!orb_pos) {
        perror("Memory allocation failed for orbit positions");
        free(orb->points);
        free(orb);
        return EXIT_FAILURE;
    }
    
    // 计算轨道位置
    int orb_count = calculate_orbit_positions(orb, orb_pos, ts, t1, nrec, npad);
    if (orb_count < 0) {
        fprintf(stderr, "Error calculating orbit positions\n");
        free(orb_pos);
        free(orb->points);
        free(orb);
        return EXIT_FAILURE;
    }
    
    // 读取输入数据（来自stdin）
    size_t capacity = INIT_CAPACITY;
    int npt = 0;
    double *in = malloc(sizeof(double) * 3 * capacity);
    if (!in) {
        perror("Memory allocation failed for input buffer");
        free(orb_pos);
        free(orb->points);
        free(orb);
        return EXIT_FAILURE;
    }
    
    double lon, lat, h;
    while (scanf("%lf %lf %lf", &lon, &lat, &h) == 3) {
        if (npt >= capacity) {
            capacity *= 2;
            in = safe_realloc(in, sizeof(double) * 3 * capacity);
        }
        in[3 * npt + 0] = lon;
        in[3 * npt + 1] = lat;
        in[3 * npt + 2] = h;
        npt++;
    }
    
    // 检查是否有有效数据
    if (npt == 0) {
        fprintf(stderr, "Error: No valid input data read from stdin\n");
        free(in);
        free(orb_pos);
        free(orb->points);
        free(orb);
        return EXIT_FAILURE;
    }
    
    // 输出缓冲区
    double *out = malloc(sizeof(double) * 2 * npt);
    if (!out) {
        perror("Memory allocation failed for output buffer");
        free(in);
        free(orb_pos);
        free(orb->points);
        free(orb);
        return EXIT_FAILURE;
    }
    
    // 预计算参数
    CalcParams params;
    params.dr = 0.5 * SOL / prm.fs;
    params.inv_dr = 1.0 / params.dr;
    params.prf = prm.prf;
    params.t1 = t1;
    params.near_range = prm.near_range;
    params.ra = prm.ra;
    params.fll = (prm.ra - prm.rc) / prm.ra;
    params.RE = prm.RE;
    
    // 处理所有点
    process_points(in, out, npt, orb_pos, orb_count, &params);
    
    // 输出结果
    for (int i = 0; i < npt; i++) {
        double rg = out[2 * i + 0];
        double az = out[2 * i + 1];
        double h  = in[3 * i + 2];
        double lat = in[3 * i + 1];
        double lon = in[3 * i + 0];
        
        if (binary_short) {
            int16_t buf[5];
            buf[0] = (int16_t)round(rg);
            buf[1] = (int16_t)round(az);
            buf[2] = (int16_t)round(h);
            buf[3] = (int16_t)round(lon * 1000);
            buf[4] = (int16_t)round(lat * 1000);
            fwrite(buf, sizeof(int16_t), 5, stdout);
        } else if (binary_double) {
            double buf[5] = {rg, az, h, lon, lat};
            fwrite(buf, sizeof(double), 5, stdout);
        } else {
            printf("%.6f %.6f %.6f %.6f %.6f\n", rg, az, h, lon, lat);
        }
    }
    
    // 清理内存
    free(in);
    free(out);
    free(orb_pos);
    if (orb->points) free(orb->points);
    free(orb);
    
    return EXIT_SUCCESS;
}