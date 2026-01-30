#  special code for LT1
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <float.h>

#define MAX_POINTS 5000
#define MAX_DEGREE 6  // 最高支持6阶

// 数据点结构
typedef struct {
    double x;      // 原始X坐标
    double dx;     // X方向位移
    double y;      // 原始Y坐标
    double dy;     // Y方向位移
    double r;      // 相关系数
    double weight; // 权重
    int valid;     // 是否有效点
} DataPoint;

// 矩阵结构
typedef struct {
    int rows;
    int cols;
    double **data;
} Matrix;

// 智能四舍五入函数，考虑数值精度
long smart_round(double value) {
    if (fabs(value) < 1e-9) {
        return 0;
    }
    
    // 计算这个值的不确定度（基于double的机器精度）
    double uncertainty = fabs(value) * DBL_EPSILON;
    
    // 如果不确定度小于0.1，说明可以安全四舍五入
    if (uncertainty < 0.1) {
        // 标准四舍五入
        return (long)round(value);
    }
    
    // 如果不确定度较大，直接四舍五入
    return (long)round(value);
}

// 创建矩阵
Matrix* create_matrix(int rows, int cols) {
    Matrix *mat = (Matrix*)malloc(sizeof(Matrix));
    mat->rows = rows;
    mat->cols = cols;
    mat->data = (double**)malloc(rows * sizeof(double*));
    for (int i = 0; i < rows; i++) {
        mat->data[i] = (double*)malloc(cols * sizeof(double));
        memset(mat->data[i], 0, cols * sizeof(double));
    }
    return mat;
}

// 释放矩阵
void free_matrix(Matrix *mat) {
    for (int i = 0; i < mat->rows; i++) {
        free(mat->data[i]);
    }
    free(mat->data);
    free(mat);
}

// 高斯消元法
int gauss_elimination(Matrix *A, double *b, double *x, int n) {
    double **aug = (double**)malloc(n * sizeof(double*));
    for (int i = 0; i < n; i++) {
        aug[i] = (double*)malloc((n + 1) * sizeof(double));
        for (int j = 0; j < n; j++) {
            aug[i][j] = A->data[i][j];
        }
        aug[i][n] = b[i];
    }
    
    for (int i = 0; i < n; i++) {
        // 寻找主元
        int max_row = i;
        for (int k = i + 1; k < n; k++) {
            if (fabs(aug[k][i]) > fabs(aug[max_row][i])) {
                max_row = k;
            }
        }
        
        if (max_row != i) {
            double *temp = aug[i];
            aug[i] = aug[max_row];
            aug[max_row] = temp;
        }
        
        if (fabs(aug[i][i]) < 1e-15) {
            for (int j = 0; j < n; j++) free(aug[j]);
            free(aug);
            return 0;
        }
        
        double pivot = aug[i][i];
        for (int j = i; j <= n; j++) {
            aug[i][j] /= pivot;
        }
        
        for (int k = i + 1; k < n; k++) {
            double factor = aug[k][i];
            for (int j = i; j <= n; j++) {
                aug[k][j] -= factor * aug[i][j];
            }
        }
    }
    
    for (int i = n - 1; i >= 0; i--) {
        x[i] = aug[i][n];
        for (int j = i + 1; j < n; j++) {
            x[i] -= aug[i][j] * x[j];
        }
    }
    
    for (int i = 0; i < n; i++) free(aug[i]);
    free(aug);
    return 1;
}

// 加权多项式拟合
int weighted_polyfit(double *x, double *y, double *weights, int n, int degree, double *coeffs) {
    if (n <= degree) return 0;
    
    int m = degree + 1;
    Matrix *A = create_matrix(m, m);
    double *b = (double*)calloc(m, sizeof(double));
    
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < m; j++) {
            double sum = 0.0;
            for (int k = 0; k < n; k++) {
                sum += weights[k] * pow(x[k], i + j);
            }
            A->data[i][j] = sum;
        }
    }
    
    for (int i = 0; i < m; i++) {
        double sum = 0.0;
        for (int k = 0; k < n; k++) {
            sum += weights[k] * y[k] * pow(x[k], i);
        }
        b[i] = sum;
    }
    
    int success = gauss_elimination(A, b, coeffs, m);
    
    free_matrix(A);
    free(b);
    return success;
}

// 读取数据
DataPoint* read_data(const char *filename, int *n_points) {
    FILE *fp = fopen(filename, "r");
    if (!fp) {
        return NULL;
    }
    
    DataPoint *points = (DataPoint*)malloc(MAX_POINTS * sizeof(DataPoint));
    *n_points = 0;
    
    while (*n_points < MAX_POINTS) {
        DataPoint p;
        if (fscanf(fp, "%lf %lf %lf %lf %lf", 
                   &p.x, &p.dx, &p.y, &p.dy, &p.r) == 5) {
            p.valid = 1;
            p.weight = fabs(p.r);
            if (p.weight < 0.01) p.weight = 0.01;
            
            points[*n_points] = p;
            (*n_points)++;
        } else {
            break;
        }
    }
    
    fclose(fp);
    return points;
}

// 综合评估数据点质量并剔除异常值
void assess_and_remove_outliers(DataPoint *points, int n_points, 
                               double *coeffs_x, double *coeffs_y, 
                               int degree_x, int degree_y) {
    
    double *residuals_x = (double*)malloc(n_points * sizeof(double));
    double *weighted_residuals_x = (double*)malloc(n_points * sizeof(double));
    double *scores_x = (double*)malloc(n_points * sizeof(double));
    
    double *residuals_y = (double*)malloc(n_points * sizeof(double));
    double *weighted_residuals_y = (double*)malloc(n_points * sizeof(double));
    double *scores_y = (double*)malloc(n_points * sizeof(double));
    
    double *combined_scores = (double*)malloc(n_points * sizeof(double));
    
    // 计算X方向残差
    for (int i = 0; i < n_points; i++) {
        if (!points[i].valid) continue;
        
        double dx_fit = 0.0;
        for (int j = 0; j <= degree_x; j++) {
            dx_fit += coeffs_x[j] * pow(points[i].x, j);
        }
        residuals_x[i] = points[i].dx - dx_fit;
        weighted_residuals_x[i] = residuals_x[i] / points[i].weight;
    }
    
    // 计算Y方向残差
    for (int i = 0; i < n_points; i++) {
        if (!points[i].valid) continue;
        
        double dy_fit = 0.0;
        for (int j = 0; j <= degree_y; j++) {
            dy_fit += coeffs_y[j] * pow(points[i].y, j);
        }
        residuals_y[i] = points[i].dy - dy_fit;
        weighted_residuals_y[i] = residuals_y[i] / points[i].weight;
    }
    
    // 计算统计量
    double mean_x = 0.0, std_x = 0.0, count_x = 0.0;
    double mean_y = 0.0, std_y = 0.0, count_y = 0.0;
    
    for (int i = 0; i < n_points; i++) {
        if (!points[i].valid) continue;
        mean_x += weighted_residuals_x[i];
        mean_y += weighted_residuals_y[i];
        count_x += 1.0;
        count_y += 1.0;
    }
    mean_x /= count_x;
    mean_y /= count_y;
    
    for (int i = 0; i < n_points; i++) {
        if (!points[i].valid) continue;
        std_x += (weighted_residuals_x[i] - mean_x) * (weighted_residuals_x[i] - mean_x);
        std_y += (weighted_residuals_y[i] - mean_y) * (weighted_residuals_y[i] - mean_y);
    }
    std_x = sqrt(std_x / (count_x - 1));
    std_y = sqrt(std_y / (count_y - 1));
    
    // 计算综合得分
    for (int i = 0; i < n_points; i++) {
        if (!points[i].valid) {
            combined_scores[i] = 0.0;
            continue;
        }
        
        double z_x = (weighted_residuals_x[i] - mean_x) / std_x;
        double z_y = (weighted_residuals_y[i] - mean_y) / std_y;
        double r_score = 1.0 / (fabs(points[i].r) + 0.1);
        
        scores_x[i] = fabs(z_x) * r_score;
        scores_y[i] = fabs(z_y) * r_score;
        combined_scores[i] = (scores_x[i] > scores_y[i]) ? scores_x[i] : scores_y[i];
    }
    
    // 找出得分的分布
    double *sorted_scores = (double*)malloc(count_x * sizeof(double));
    int idx = 0;
    for (int i = 0; i < n_points; i++) {
        if (points[i].valid) {
            sorted_scores[idx++] = combined_scores[i];
        }
    }
    
    // 排序（简单冒泡）
    for (int i = 0; i < idx - 1; i++) {
        for (int j = 0; j < idx - i - 1; j++) {
            if (sorted_scores[j] > sorted_scores[j + 1]) {
                double temp = sorted_scores[j];
                sorted_scores[j] = sorted_scores[j + 1];
                sorted_scores[j + 1] = temp;
            }
        }
    }
    
    // 计算阈值
    double q1 = sorted_scores[(int)(idx * 0.25)];
    double q3 = sorted_scores[(int)(idx * 0.75)];
    double iqr = q3 - q1;
    double threshold = q3 + 1.5 * iqr;
    
    // 标记异常值
    for (int i = 0; i < n_points; i++) {
        if (!points[i].valid) continue;
        if (combined_scores[i] > threshold) {
            points[i].valid = 0;
        }
    }
    
    free(residuals_x);
    free(weighted_residuals_x);
    free(scores_x);
    free(residuals_y);
    free(weighted_residuals_y);
    free(scores_y);
    free(combined_scores);
    free(sorted_scores);
}

// 获取有效数据
void get_valid_data(DataPoint *points, int n_points, 
                   double **x, double **dx, double **weights_x, int *n_x,
                   double **y, double **dy, double **weights_y, int *n_y) {
    
    *n_x = 0;
    *n_y = 0;
    for (int i = 0; i < n_points; i++) {
        if (points[i].valid) {
            (*n_x)++;
            (*n_y)++;
        }
    }
    
    *x = (double*)malloc(*n_x * sizeof(double));
    *dx = (double*)malloc(*n_x * sizeof(double));
    *weights_x = (double*)malloc(*n_x * sizeof(double));
    *y = (double*)malloc(*n_y * sizeof(double));
    *dy = (double*)malloc(*n_y * sizeof(double));
    *weights_y = (double*)malloc(*n_y * sizeof(double));
    
    int idx = 0;
    for (int i = 0; i < n_points; i++) {
        if (points[i].valid) {
            (*x)[idx] = points[i].x;
            (*dx)[idx] = points[i].dx;
            (*weights_x)[idx] = points[i].weight;
            (*y)[idx] = points[i].y;
            (*dy)[idx] = points[i].dy;
            (*weights_y)[idx] = points[i].weight;
            idx++;
        }
    }
}

// 计算误差统计
void calculate_errors(double *x_data, double *dx_data, double *weights_x,
                     double *y_data, double *dy_data, double *weights_y,
                     double *coeffs_x, double *coeffs_y, int n, int degree,
                     double *rms_error_x, double *rms_error_y,
                     double *weighted_rms_error_x, double *weighted_rms_error_y) {
    
    double total_error_x = 0.0;
    double total_weighted_error_x = 0.0;
    double total_weight_x = 0.0;
    
    double total_error_y = 0.0;
    double total_weighted_error_y = 0.0;
    double total_weight_y = 0.0;
    
    for (int i = 0; i < n; i++) {
        // X方向拟合值
        double dx_fit = 0.0;
        for (int j = 0; j <= degree; j++) {
            dx_fit += coeffs_x[j] * pow(x_data[i], j);
        }
        
        // Y方向拟合值
        double dy_fit = 0.0;
        for (int j = 0; j <= degree; j++) {
            dy_fit += coeffs_y[j] * pow(y_data[i], j);
        }
        
        // X方向误差
        double error_x = dx_data[i] - dx_fit;
        total_error_x += error_x * error_x;
        total_weighted_error_x += weights_x[i] * error_x * error_x;
        total_weight_x += weights_x[i];
        
        // Y方向误差
        double error_y = dy_data[i] - dy_fit;
        total_error_y += error_y * error_y;
        total_weighted_error_y += weights_y[i] * error_y * error_y;
        total_weight_y += weights_y[i];
    }
    
    *rms_error_x = sqrt(total_error_x / n);
    *rms_error_y = sqrt(total_error_y / n);
    *weighted_rms_error_x = sqrt(total_weighted_error_x / total_weight_x);
    *weighted_rms_error_y = sqrt(total_weighted_error_y / total_weight_y);
}

// 按照指定格式输出系数和误差
void output_coefficients_and_errors(double *coeffs_x, double *coeffs_y, int degree,
                                   double rms_error_x, double rms_error_y,
                                   double weighted_rms_error_x, double weighted_rms_error_y) {
    // 输出X方向系数 a0-a5
    for (int i = 0; i <= MAX_DEGREE; i++) {
        if (i <= degree) {
            if (i == 0) {
                // a0输出整数
                //printf("%ld ", smart_round(coeffs_x[i]));
                printf("%.12e ", coeffs_x[i]);
            } else {
                // a1-a5输出浮点数
                printf("%.12e ", coeffs_x[i]);
            }
        } else {
            // 高阶项为0
            printf("0 ");
        }
    }
    
    // 输出Y方向系数 b0-b5
    for (int i = 0; i <= MAX_DEGREE; i++) {
        if (i <= degree) {
            if (i == 0) {
                // b0输出整数
                //printf("%ld ", smart_round(coeffs_y[i]));
                printf("%.12e ", coeffs_y[i]);
            } else {
                // b1-b5输出浮点数
                printf("%.12e ", coeffs_y[i]);
            }
        } else {
            // 高阶项为0
            printf("0 ");
        }
    }
    
    // 输出误差
    printf("%.12e %.12e\n", weighted_rms_error_x, weighted_rms_error_y);
}

int main(int argc, char *argv[]) {
    if (argc != 3) {
        fprintf(stderr, "For Lt1 Dataset, can be used for other data. \n");
        fprintf(stderr, "Usage: %s <degree> <data_file>\n", argv[0]);
        fprintf(stderr, "Degree must be between 0 and %d\n", MAX_DEGREE);
        return 1;
    }
    
    int degree = atoi(argv[1]);
    const char *filename = argv[2];
    
    if (degree < 0 || degree > MAX_DEGREE) {
        fprintf(stderr, "Error: Degree must be between 0 and %d\n", MAX_DEGREE);
        return 1;
    }
    
    // 读取原始数据
    int n_total;
    DataPoint *points = read_data(filename, &n_total);
    if (!points || n_total == 0) {
        fprintf(stderr, "Error: No data read from file\n");
        return 1;
    }
    
    double *coeffs_x = (double*)calloc(MAX_DEGREE + 1, sizeof(double));
    double *coeffs_y = (double*)calloc(MAX_DEGREE + 1, sizeof(double));
    
    // 迭代拟合（5次）
    for (int iteration = 1; iteration <= 5; iteration++) {
        // 获取当前有效数据
        double *x_data, *dx_data, *weights_x;
        double *y_data, *dy_data, *weights_y;
        int n_x, n_y;
        get_valid_data(points, n_total, &x_data, &dx_data, &weights_x, &n_x,
                       &y_data, &dy_data, &weights_y, &n_y);
        
        if (n_x <= degree) {
            free(x_data); free(dx_data); free(weights_x);
            free(y_data); free(dy_data); free(weights_y);
            break;
        }
        
        // X方向加权拟合
        weighted_polyfit(x_data, dx_data, weights_x, n_x, degree, coeffs_x);
        
        // Y方向加权拟合
        weighted_polyfit(y_data, dy_data, weights_y, n_y, degree, coeffs_y);
        
        // 如果不是最后一次迭代，评估并剔除异常值
        if (iteration < 5) {
            assess_and_remove_outliers(points, n_total, coeffs_x, coeffs_y, degree, degree);
        }
        
        // 计算误差
        double rms_error_x, rms_error_y, weighted_rms_error_x, weighted_rms_error_y;
        calculate_errors(x_data, dx_data, weights_x, 
                        y_data, dy_data, weights_y,
                        coeffs_x, coeffs_y, n_x, degree,
                        &rms_error_x, &rms_error_y,
                        &weighted_rms_error_x, &weighted_rms_error_y);
        
        // 释放临时数组
        free(x_data); free(dx_data); free(weights_x);
        free(y_data); free(dy_data); free(weights_y);
        
        // 最后一次迭代输出结果
        if (iteration == 5) {
            // 重新获取最终的有效数据进行误差计算
            get_valid_data(points, n_total, &x_data, &dx_data, &weights_x, &n_x,
                           &y_data, &dy_data, &weights_y, &n_y);
            
            // 最终误差计算
            calculate_errors(x_data, dx_data, weights_x, 
                            y_data, dy_data, weights_y,
                            coeffs_x, coeffs_y, n_x, degree,
                            &rms_error_x, &rms_error_y,
                            &weighted_rms_error_x, &weighted_rms_error_y);
            
            // 输出最终系数和误差
            output_coefficients_and_errors(coeffs_x, coeffs_y, degree,
                                          rms_error_x, rms_error_y,
                                          weighted_rms_error_x, weighted_rms_error_y);
            
            free(x_data); free(dx_data); free(weights_x);
            free(y_data); free(dy_data); free(weights_y);
        }
    }
    
    // 清理内存
    free(coeffs_x);
    free(coeffs_y);
    free(points);
    
    return 0;
}