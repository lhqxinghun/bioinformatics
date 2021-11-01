/*
 * matutils.c
 *
 *  Created on: Jun 19, 2014
 *      Author: Wenyuan Li
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "errlog.h"
#include "matutils.h"

/*---------------------------------
    Routines for VEC_DOUBLE
-----------------------------------*/
void create_VEC_DOUBLE(int n, VEC_DOUBLE* vec)
{
    vec->n = n;
    vec->v = (double*)malloc( n*sizeof(double) );
    if ( vec->v == NULL )
    {
        Errlog timerr;
        timerr.errlog("Application for memory failure;");
        exit(EXIT_FAILURE);
    }
}

void create_ones_VEC_DOUBLE(int n, VEC_DOUBLE* vec)
{
    int i;
    vec->n = n;
    vec->v = (double*)malloc( n*sizeof(double) );
    if ( vec->v == NULL )
    {
        Errlog timerr;
        timerr.errlog("Application for memory failure;");
        exit(EXIT_FAILURE);
    }
    for ( i=0; i<n; i++ )
        vec->v[i] = 1;
}

void init_VEC_DOUBLE(VEC_DOUBLE* vec)
{
    vec->n = 0;
    vec->v = NULL;
}

void free_VEC_DOUBLE(VEC_DOUBLE* vec)
{
    if (vec->v!=NULL)
        free( vec->v );
}

void elementwise_multi_vector_VEC_DOUBLE(VEC_DOUBLE* dst, int start_idx, int end_idx, VEC_DOUBLE* src)
{
    int i, cnt_elements=end_idx-start_idx+1;
    if (cnt_elements!=src->n)
    {        
        Errlog timerr;
        char msg[MAXCHAR];
        sprintf(msg,"Err: Elementwise_multi_vector_VEC_FLOAT:\n\tlength of src (%d) doesn't match with (start_idx-end_idx+1)=%d;", src->n, cnt_elements);
        timerr.errlog(msg);
    }
    for (i=0; i<cnt_elements; i++)
        dst->v[i+start_idx-1] *= src->v[i];
}
/*---------------------------------
    Routines for MATRIX_DOUBLE
-----------------------------------*/
void init_MATRIX_DOUBLE(MATRIX_DOUBLE* mat)
{
    mat->start_row = 0;
    mat->end_row = 0;
    mat->cnt_row = 0;
    mat->total_column = 0;
    mat->v = NULL;
}

void create_MATRIX_DOUBLE(int cnt_row, int total_column, MATRIX_DOUBLE* mat)
{
    int i;
    mat->start_row = 1;
    mat->end_row = cnt_row;
    mat->cnt_row = cnt_row;
    mat->total_column = total_column;
    mat->v = (double**)malloc( cnt_row*sizeof(double*) );
    if ( mat->v == NULL )
    {
        Errlog timerr;
        timerr.errlog("Application for memory failure;");
        exit(EXIT_FAILURE);
    }
    for ( i=0; i<cnt_row; i++ )
    {
        mat->v[i] = (double*)malloc( total_column*sizeof(double) );
        if ( mat->v[i] == NULL )
        {
            Errlog timerr;
            timerr.errlog("Application for memory failure;");
            exit(EXIT_FAILURE);
        }
    }
}

void copy_MATRIX_DOUBLE(double* psourdata, int width, int height, MATRIX_DOUBLE* mat)
{
    int i, j, total_column, cnt_row;

    total_column = mat->total_column;
    cnt_row = mat->cnt_row;
    if(width == total_column && height == cnt_row)
    {
        for ( i=0; i<cnt_row; i++ )
        {
            for ( j=0; j<total_column; j++ )
                mat->v[i][j] = *(psourdata + i * total_column + j);
        }
    }
}

void copy_DEST_DOUBLE(MATRIX_DOUBLE* mat, int width, int height, double* pdestdata)
{
    int i, j;

    if(width == mat->total_column && height == mat->cnt_row)
    {
        for ( i=0; i<height; i++ )
        {
            for ( j=0; j<width; j++ )
                *(pdestdata + i * width + j) = mat->v[i][j];
        }
    }
}

void free_MATRIX_DOUBLE(MATRIX_DOUBLE* mat)
{
    int i;
    if (mat->v!=NULL)
    {
        for ( i=0; i<mat->cnt_row; i++ )
            if (mat->v[i]!=NULL) free( mat->v[i] );
        free( mat->v );
    }
}

double normp_double(double* v, int len, double p)
{
    int i;
    double ret=0;
    for (i=0; i<len; i++)
        ret += pow(v[i], p);
    return pow(ret,1/p);
}

double norm1_double(double* v, int len)
{
    int i;
    double ret=0;
    for (i=0; i<len; i++)
    {
        ret+=v[i];
    }
    return ret;
}

double max_double(double* v, int len)
{
    int i;
    double max=v[0];
    for (i=0; i<len; i++)
        if (v[i]>max) max=v[i];
    return max;
}

void norm_rows_MATRIX_DOUBLE(MATRIX_DOUBLE* mat, float p, VEC_DOUBLE* norms)
{
    int i, cnt_zeros = 0;
    if (norms->n != mat->cnt_row)
    {
        free_VEC_DOUBLE(norms);
        init_VEC_DOUBLE(norms);
        create_VEC_DOUBLE(mat->cnt_row, norms);
    }
    for (i = 0; i < mat->cnt_row; i++)
    {
        if (p == -1)
            norms->v[i] = max_double(mat->v[i], mat->total_column);
        else if (p==1)
            norms->v[i] = norm1_double(mat->v[i], mat->total_column);
        else if (p>0)
            norms->v[i] = normp_double(mat->v[i], mat->total_column, p);
        else
        {
            Errlog timerr;
            timerr.errlog("Error: p<0 and p!=-1 is wrong;");
            exit(EXIT_FAILURE);
        }
        if (norms->v[i]==0)
            cnt_zeros++;
    }
    if (cnt_zeros!=0)
    {
        Errlog timerr;
        char msg[MAXCHAR];
        sprintf(msg,"Warning: %d rows are all zeros;",cnt_zeros);
        timerr.errlog(msg);
    }
}

void elementwise_div_value_MATRIX_DOUBLE(MATRIX_DOUBLE* mat, double v)
{
    int i,j;
    if (v==0) return;
    for (i=0; i<mat->cnt_row; i++)
    {
        for (j=0; j<mat->total_column; j++)
            mat->v[i][j] /= v;
    }
}

void elementwise_div_vector_MATRIX_DOUBLE(MATRIX_DOUBLE* mat, VEC_DOUBLE* norms)
{
    int i,j;
    if(mat->total_column != norms->n)
    {
        Errlog timerr;
        char msg[MAXCHAR];
        sprintf(msg,"Error 'elementwise_div_vector_MATRIX_FLOAT'\n\tLengths of matrix (%d) and vector (%d) doesn't match.\n", mat->total_column, norms->n);
        timerr.errlog(msg);
        exit(EXIT_FAILURE);
    }
    for (i=0; i<mat->cnt_row; i++)
    {
        if (norms->v[i+mat->start_row-1]!=0)
        {
            for (j=0; j<mat->total_column; j++)
                if (norms->v[j]!=0)
                    mat->v[i][j] /= norms->v[i+mat->start_row-1]*norms->v[j];
        }
    }
}
