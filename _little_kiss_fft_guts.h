// LITTLE KISS FFT, μcontroller-friendly port by Bram Stolk of KISS FFT.

/*
 *  Copyright (c) 2003-2010, Mark Borgerding. All rights reserved.
 *  This file is part of KISS FFT - https://github.com/mborgerding/kissfft
 *
 *  SPDX-License-Identifier: BSD-3-Clause
 *  See COPYING file for more information.
 */

#ifndef _kiss_fft_guts_h
#define _kiss_fft_guts_h

#include "little_kiss_fft.h"
#include <limits.h>
#include <stdint.h>


#define MAXFACTORS 32
/* e.g. an fft of length 128 has 4 factors
 as far as kissfft is concerned
 4*4*4*2
 */

struct kiss_fft_state{
    int nfft;
    int inverse;
    int factors[2*MAXFACTORS];
    kiss_fft_cpx twiddles[1];
};

/*
  Explanation of macros dealing with complex math:

   C_MUL(m,a,b)         : m = a*b
   C_FIXDIV( c , div )  : if a fixed point impl., c /= div. noop otherwise
   C_SUB( res, a,b)     : res = a - b
   C_SUBFROM( res , a)  : res -= a
   C_ADDTO( res , a)    : res += a
 * */

#if (FIXED_POINT==32)
#	define FRACBITS 30 // was: 31 (Somehow this does not work!)
#	define SAMPPROD int64_t
#	define SAMP_MAX  1073741823 // was: INT32_MAX
#	define SAMP_MIN -1073741824 // was: INT32_MIN
#else
#	define FRACBITS 15
#	define SAMPPROD int32_t
#	define SAMP_MAX INT16_MAX
#	define SAMP_MIN INT16_MIN
#endif

#define smul(a,b) ( (SAMPPROD)(a)*(b) )
#define sround( x )  (kiss_fft_scalar)( ( (x) + (1<<(FRACBITS-1)) ) >> FRACBITS )

#define S_MUL(a,b) sround( smul(a,b) )

#define C_MUL(m,a,b) \
      do{ (m).r = sround( smul((a).r,(b).r) - smul((a).i,(b).i) ); \
          (m).i = sround( smul((a).r,(b).i) + smul((a).i,(b).r) ); }while(0)

#define DIVSCALAR(x,k) \
    (x) = sround( smul(  x, SAMP_MAX/k ) )

#define C_FIXDIV(c,div) \
    do {    DIVSCALAR( (c).r , div);  \
        DIVSCALAR( (c).i  , div); }while (0)

#define C_MULBYSCALAR( c, s ) \
    do{ (c).r =  sround( smul( (c).r , s ) ) ;\
        (c).i =  sround( smul( (c).i , s ) ) ; }while(0)


#define  C_ADD( res, a,b)\
    do { \
        (res).r=(a).r+(b).r;  (res).i=(a).i+(b).i; \
    }while(0)
#define  C_SUB( res, a,b)\
    do { \
        (res).r=(a).r-(b).r;  (res).i=(a).i-(b).i; \
    }while(0)
#define C_ADDTO( res , a)\
    do { \
        (res).r += (a).r;  (res).i += (a).i;\
    }while(0)

#define C_SUBFROM( res , a)\
    do {\
        (res).r -= (a).r;  (res).i -= (a).i; \
    }while(0)


#define KISS_FFT_COS(phase)  floorf(.5+SAMP_MAX * cosf (phase))
#define KISS_FFT_SIN(phase)  floorf(.5+SAMP_MAX * sinf (phase))
#define HALF_OF(x) ((x)>>1)

#define  kf_cexp(x,phase) \
    do{ \
        (x)->r = KISS_FFT_COS(phase);\
        (x)->i = KISS_FFT_SIN(phase);\
    }while(0)


#endif /* _little_kiss_fft_guts_h */
