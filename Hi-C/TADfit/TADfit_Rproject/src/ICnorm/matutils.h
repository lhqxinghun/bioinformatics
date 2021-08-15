/*
 * matutils.h
 *
 *  Created on: Jun 19, 2014
 *      Author: Wenyuan Li
 */

#ifndef MATUTILS_H_
#define MATUTILS_H_

/*---------------------------
  Data Types
---------------------------*/
typedef int BOOL;
#define TRUE 1
#define FALSE 0

/*---------------------------
  Constants
---------------------------*/
#define MAXCHAR 1024
#define MAX_LINE_LENGTH 1000000
#define MAXCHAR_ID 64
#define MEM_NOT_CONSIDERED 0  // in MB

/*---------------------------
  Data structure
---------------------------*/
typedef struct _vector_i {
	int n;
	int* v;
} VEC_INT;

typedef struct _vector_i2 {
	int n;
	int* v1;
	int* v2;
} VEC_INT2;

typedef struct _vector_f {
	int n;
	float* v;
} VEC_FLOAT;

typedef struct _vector_d {
	int n;
	double* v;
} VEC_DOUBLE;

typedef struct _matrix_f {
	int start_row;
	int end_row;
	int cnt_row;
	int total_column;
	int real_cnt_row;
	float** v;
} MATRIX_FLOAT;

typedef struct _matrix_d {
	int start_row;
	int end_row;
	int cnt_row;
	int total_column;
	double** v;
} MATRIX_DOUBLE;

/*---------------------------------
    Routines for VEC_DOUBLE
-----------------------------------*/
void create_VEC_DOUBLE(int n, VEC_DOUBLE* vec);
void create_ones_VEC_DOUBLE(int n, VEC_DOUBLE* vec);
void init_VEC_DOUBLE(VEC_DOUBLE* vec);
void free_VEC_DOUBLE(VEC_DOUBLE* vec);
void elementwise_multi_vector_VEC_DOUBLE(VEC_DOUBLE* dst, int start_idx, int end_idx, VEC_DOUBLE* src);

/*---------------------------------
    Routines for MATRIX_DOUBLE
-----------------------------------*/
void init_MATRIX_DOUBLE(MATRIX_DOUBLE* mat);
void create_MATRIX_DOUBLE(int cnt_row, int total_column, MATRIX_DOUBLE* mat);
void copy_MATRIX_DOUBLE(double* psourdata, int width, int height, MATRIX_DOUBLE* mat);
void copy_DEST_DOUBLE(MATRIX_DOUBLE* mat, int width, int height, double* pdestdata);
void free_MATRIX_DOUBLE(MATRIX_DOUBLE* mat);
double normp_double(double* v, int len, double p);
double norm1_double(double* v, int len);
double max_double(double* v, int len);
void norm_rows_MATRIX_DOUBLE(MATRIX_DOUBLE* mat, float p, VEC_DOUBLE* norms);
void elementwise_div_value_MATRIX_DOUBLE(MATRIX_DOUBLE* mat, double v);
void elementwise_div_vector_MATRIX_DOUBLE(MATRIX_DOUBLE* mat, VEC_DOUBLE* norms);

#endif
