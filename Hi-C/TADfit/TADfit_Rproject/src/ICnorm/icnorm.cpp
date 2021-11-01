/*
 * Normalization using iterative correction (abbreviated as IC) algorithm
 *    
 *
 * 
 */
#include "icnorm.h"
#include "matutils.h"

void normalizationbyic(double *psourdata, int width, int height, int maxiter, double experowsum, double *pdestdata)
{
    int i;
    double p = 1.0;
    double rowsum;
    MATRIX_DOUBLE mat;
    VEC_DOUBLE bias_vec, norms_vec;

    create_ones_VEC_DOUBLE(width, &bias_vec);
    create_ones_VEC_DOUBLE(width, &norms_vec);
    init_MATRIX_DOUBLE(&mat);
    create_MATRIX_DOUBLE(height, width, &mat);
    copy_MATRIX_DOUBLE(psourdata, width, height, &mat);
    for(i = 0; i < maxiter; i++)
    {
        // get new norms from updated matrix: t = normp(O,p);
        norm_rows_MATRIX_DOUBLE(&mat, p, &norms_vec);
        // update matrix: O = diag(1./t)*O*diag(1./t);
        elementwise_div_vector_MATRIX_DOUBLE(&mat, &norms_vec);
        // update the bias vector: b = b.*t;
        elementwise_multi_vector_VEC_DOUBLE(&bias_vec, 1, width, &norms_vec);
    }
    copy_MATRIX_DOUBLE(psourdata, width, height, &mat);
    // update matrix: B = diag(1./d_save)*A*diag(1./d_save);
    elementwise_div_vector_MATRIX_DOUBLE(&mat, &bias_vec);
    rowsum = norm1_double(mat.v[0], mat.total_column) / experowsum;
    elementwise_div_value_MATRIX_DOUBLE(&mat, rowsum);
    copy_DEST_DOUBLE(&mat, width, height, pdestdata);
    // free space
    free_MATRIX_DOUBLE(&mat);
    free_VEC_DOUBLE(&bias_vec);
    free_VEC_DOUBLE(&norms_vec);
}
