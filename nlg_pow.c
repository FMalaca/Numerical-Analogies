/* Copyright (c) 2024, Yves Lepage */
#define MODULE "nlg_pow"
#define TRACE 0
/* extern int tracet ; */

/*
 * Compute the analogical power.
 *
 * For instance:
 *     5 : 7 :: 10 : 12  =>  p =  1.0   arithmetic analogy
 *     2 : 4 ::  8 : 16  =>  p =  0.0   geometric analogy
 *   8/5 : 2 ::  4 :  8  =>  p = -1.0   harmonic analogy
 *     2 : 2 ::  3 :  4  =>  p = -100   minimum analogy
 *
 * Method:
 *  dichotomic search
 *
 * start from minus infinity and plus infinity,
 *  i.e., large negative and positive powers
 * compute the generalized means of the extremes and the means
 * compare them
 * and recursively progress with the middle of the powers,
 * stop when the difference is smaller than a small value EPSILON.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
/*
#include "usrlib.h"
#include "usrmath.h"
*/
#include "nlg_pow.h"

/*
 * REPLACES #include "usrlib.h"
 * macro for trace:
 * use:
 * trace((string))
 * in C code not followed by semi-colon
 */

#define ytrace(x) { printf x ; fflush(stderr) ; fflush(stdout) ; }
#define ntrace(x)

#if TRACE
/* #define trace(x) { ytrace(x) ; if ( iscgiuse() ) printf("<br>\n") ; fflush(stdout) ; } */
#define trace(x) { ytrace(x) ; }
#else
#define trace(x)
#endif

/*
 * Exit on error
 */

#define ERR_MODULE_LABEL       " in module:   "
#define ERR_FUNCTION_LABEL     " in function: "
#define ERR_DIAGNOSTIC_LABEL   " diagnostic:  "

/*
 * Error, do not exit
 */

void warning(char *module_name, char *function_name, char *msg)
{
   fprintf(stderr,
          "\n#"
          "\n#" ERR_MODULE_LABEL     "%s"
          "\n#" ERR_FUNCTION_LABEL   "%s"
          "\n#" ERR_DIAGNOSTIC_LABEL "%s"
          "\n#"
          "\n",
          module_name,function_name,msg) ;
   fflush(stderr) ;
}

/*
 * Exit on error
 */

void error(char *module_name, char *function_name, char *msg)
{
    warning(module_name,function_name,msg) ;
    abort() ;
}


/*
 * Constants
 */

#define EPSILON          0.0001
#define MINUS_INFINITY  -100.0
#define PLUS_INFINITY   +100.0
#define ALL_R           -200.0

/*
 * Trace macros
 */

#define F_NLGFMT    "%f : %f :: %f : %f"

/*
 * Local functions
 */

static float generalized_mean(float, float,
                              float) ;
static float difference(float, float, float, float,
                        float) ;
static float inner_dichotomic_search(float, float, float, float,
                                      float, float,
                                      float, float) ;
static float dichotomic_search(float a, float b, float c, float d) ;

/*
 * Compute the generalized mean of two non-negative real numbers.
 */

static float generalized_mean(float a, float d, float p)
{
    float result = 0.0 ;

trace(("in  generalized_mean(%f, %f, p=%f)\n", a, d, p))

    if (a == d)
        result = a ;
    else if ( p <= MINUS_INFINITY )
        result = fmin(a, d) ;
    else if ( 0.0 == p )
        result = pow((a * d), 0.5) ;
    else if ( p >= PLUS_INFINITY )
        result = fmax(a, d) ;
    else
    {
        if ( 0.0 == a )
            result = d * pow(0.5, 1/p) ;
        else if ( 0.0 == d )
            result = a * pow(0.5, 1/p) ;
        else
            result = d * pow(0.5 * (pow(a/d, p) + 1), 1/p) ;
    } ;

trace(("out generalized_mean(%f, %f, p=%f) = %f\n", a, d, p, result))

    return result ;
}

/*
 * Compute the difference of the means of the extremes and the means
 * for an analogy between non-negative real numbers.
 */

static float difference(float a, float b, float c, float d, float p)
{
    float result = 0.0 ;

trace(("in  difference(" F_NLGFMT ", p=%f)\n", a, b, c, d, p))

    /* Error when the terms are not sorted. */
    if ( ! ( a <= b && b <= c && c <= d ) )
        error("nlg_pow", "difference", "terms not sorted.") ;
    result = generalized_mean(a, d, p) - generalized_mean(b, c, p) ;

trace(("out difference(" F_NLGFMT ", p=%f) = %f\n", a, b, c, d, p, result))

    return result ;
}

/*
 * Compute the analogical power between four ordered non-negative real numbers.
 * Recursive function.
 */

static float inner_dichotomic_search(float a, float b, float c, float d,
                              float i, float j,
                              float di, float dj)
{
    float result = 0 ;
    float m = 0.0 ;
    float dm = 0.0 ;

trace(("in  inner_dichotomic_search(" F_NLGFMT ",\n\t%f, %f,\n\t%f, %f)\n", a, b, c, d, i, j, di, dj))

    if ( ! ( i < j ) )
        error("nlg_pow", "inner_dichotomic_search", "i >= j.") ;
    if ( ! ( di <= 0 && 0 <= dj ) )
        error("nlg_pow", "inner_dichotomic_search", "di >= 0 or dj <= 0.") ;
    if ( fabs(i-j) < EPSILON ) {
        if ( fabs(di) < fabs(dj) )
            result = i ;
        else
            result = j ;
    }
    else
    {
        m = (i+j) / 2 ;
        dm = difference(a, b, c, d, m) ;
        trace(("mid inner_dichotomic_search(...) %f=%f < %f=%f < %f=%f\n",i,di,m,dm,j,dj))
        if ( 0.0 < dm )
            result = inner_dichotomic_search(a, b, c, d, i, m, di, dm) ;
        else
            result = inner_dichotomic_search(a, b, c, d, m, j, dm, dj) ;
    } ;

trace(("out inner_dichotomic_search(" F_NLGFMT ",\n\t%f, %f,\n\t%f, %f) = %f\n", a, b, c, d, i, j, di, dj, result))

    return result ;
}

/*
 * Compute the analogical power between four ordered non-negative real numbers.
 * This is only for the general case where p is not 0, or the infinities.
 */

static float dichotomic_search(float a, float b, float c, float d)
{
    float result = 0.0 ;
    float i = MINUS_INFINITY, j = PLUS_INFINITY ;
    float di = 0.0, dj = 0.0 ;

trace(("in  dichotomic_search(" F_NLGFMT ")\n", a, b, c, d))

    di = difference(a, b, c, d, i) ;
    dj = difference(a, b, c, d, j) ;
    result = inner_dichotomic_search(a, b, c, d, i, j, di, dj) ;

trace(("out dichotomic_search(" F_NLGFMT ") = %f\n", a, b, c, d, result))

    return result ;
}

/*
 * Compute the analogical power between four ordered non-negative real numbers.
 * Particular cases are dealt with first.
 * If not a particular case, call the dichotomic search.
 */

float power(float a, float b, float c, float d)
{
    float result = 0.0 ;

trace(("in  power(" F_NLGFMT ")\n", a, b, c, d))

    if ( a == b && c == d )
        /* all values are possible, R and even infinities */
        result = ALL_R ;
    else if ( a == b ) /* but c != d */
        /*minimum analogy */
        result = MINUS_INFINITY ;
    else if ( c == d ) /* but a != b */
        /* maximum analogy */
        result = PLUS_INFINITY ;
    else if ( a * d == b * c )
        /* geometric analogy */
        result = 0.0 ;
    else if ( a - b == c - d )
        /* arithmetic analogy */
        result = 1.0 ;
    else if ( fabs( ((a + d) * (b * c)) - ((b + c) * (a * d)) ) < 0.000001 )
        /* harmonic analogy */
        result = -1.0 ;
    else
        result = dichotomic_search(a, b, c, d) ;

trace(("out power(" F_NLGFMT ") = %f\n", a, b, c, d, result))

    return result ;
}



