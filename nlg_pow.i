/* File : nlg_pow.i */
%module nlg_pow
%{
#include "nlg_pow.h"
/* Dichotomic search for the computation of analogical power */
extern float power(float, float, float, float) ;
%}
extern float power(float, float, float, float) ;


