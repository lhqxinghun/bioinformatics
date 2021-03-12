#include "gsl_errno.h"
#include "errlog.h"
#include <stdlib.h>

void GSL_ERROR_VAL(const char * reason, int gsl_errno, int value)
{
    Errlog timerr;
    char msg[MAXCHAR];
    sprintf(msg,"Err: GSL_ERRNO: %d, GSL_REASON: %s, GSL_VALUE: %d;", gsl_errno, reason, value);
    timerr.errlog(msg);
    exit(EXIT_FAILURE);
}

void GSL_ERROR_VOID(const char * reason, int gsl_errno)
{
    Errlog timerr;
    char msg[MAXCHAR];
    sprintf(msg,"Err: GSL_ERRNO: %d, GSL_REASON: %s;", gsl_errno, reason);
    timerr.errlog(msg);
    exit(EXIT_FAILURE);
}
