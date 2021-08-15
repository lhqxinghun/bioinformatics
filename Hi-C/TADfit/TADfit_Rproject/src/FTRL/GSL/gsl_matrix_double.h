#ifndef __GSL_MATRIX_DOUBLE_H__
#define __GSL_MATRIX_DOUBLE_H__

#include <stdlib.h>
#include "gsl_types.h"
#include "gsl_block_double.h"
#include "gsl_errno.h"
#include "gsl_inline.h"
#include "gsl_check_range.h"

#undef __BEGIN_DECLS
#undef __END_DECLS
#ifdef __cplusplus
# define __BEGIN_DECLS extern "C" {
# define __END_DECLS }
#else
# define __BEGIN_DECLS /* empty */
# define __END_DECLS /* empty */
#endif

__BEGIN_DECLS

typedef struct 
{
  size_t size1;
  size_t size2;
  size_t tda;
  double * data;
  gsl_block * block;
  int owner;
} gsl_matrix;

/* Allocation */

gsl_matrix * 
gsl_matrix_alloc (const size_t n1, const size_t n2);

gsl_matrix * 
gsl_matrix_calloc (const size_t n1, const size_t n2);

void gsl_matrix_free (gsl_matrix * m);

/* Operations */

void gsl_matrix_set_zero (gsl_matrix * m);
void gsl_matrix_set_identity (gsl_matrix * m);
void gsl_matrix_set_all (gsl_matrix * m, double x);

/* inline functions if you are using GCC */

INLINE_DECL double   gsl_matrix_get(const gsl_matrix * m, const size_t i, const size_t j);
INLINE_DECL void    gsl_matrix_set(gsl_matrix * m, const size_t i, const size_t j, const double x);

#ifdef HAVE_INLINE
INLINE_DECL
double
gsl_matrix_get(const gsl_matrix * m, const size_t i, const size_t j)
{
#if GSL_RANGE_CHECK
  if (GSL_RANGE_COND(1)) 
    {
      if (i >= m->size1)
        {
          GSL_ERROR_VAL("first index out of range", GSL_EINVAL, 0) ;
        }
      else if (j >= m->size2)
        {
          GSL_ERROR_VAL("second index out of range", GSL_EINVAL, 0) ;
        }
    }
#endif
  return m->data[i * m->tda + j] ;
} 

INLINE_DECL
void
gsl_matrix_set(gsl_matrix * m, const size_t i, const size_t j, const double x)
{
#if GSL_RANGE_CHECK
  if (GSL_RANGE_COND(1)) 
    {
      if (i >= m->size1)
        {
          GSL_ERROR_VOID("first index out of range", GSL_EINVAL) ;
        }
      else if (j >= m->size2)
        {
          GSL_ERROR_VOID("second index out of range", GSL_EINVAL) ;
        }
    }
#endif
  m->data[i * m->tda + j] = x ;
}

#endif

__END_DECLS

#endif /* __GSL_MATRIX_DOUBLE_H__ */
