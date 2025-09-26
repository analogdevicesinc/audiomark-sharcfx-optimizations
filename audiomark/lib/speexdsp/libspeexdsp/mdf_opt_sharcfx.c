/* Copyright(c) 2025 Analog Devices, Inc. All Rights Reserved. 
   Copyright (C) 2003 Epic Games (written by Jean-Marc Valin)
   Copyright (C) 2004-2006 Epic Games

   File: preprocess.c
   Preprocessor with denoising based on the algorithm by Ephraim and Malah

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions are
   met:

   1. Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.

   2. Redistributions in binary form must reproduce the above copyright
   notice, this list of conditions and the following disclaimer in the
   documentation and/or other materials provided with the distribution.

   3. The name of the author may not be used to endorse or promote products
   derived from this software without specific prior written permission.

   THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR
   IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
   OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
   DISCLAIMED. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT,
   INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
   (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
   SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
   HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
   STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
   ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
   POSSIBILITY OF SUCH DAMAGE.
*/


#include <sys/platform.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include <complex.h>
#include <matrix.h>
#include <filter.h>
#include <vector.h>

#include "math_fixedpoint_vec.h"
/* Cross-platform data type definitions. */
#include "libdsp_types.h"

//#pragma once
#include <inttypes.h>

#ifdef __XTENSA__
#include <xtensa/sim.h>
#include <xtensa/tie/xt_pdxn.h>
#endif

#include "arch.h"

#if defined(FLOATING_POINT)

#define MIN(X, Y) (((X) < (Y)) ? (X) : (Y))
#define MAX(X, Y) (((X) > (Y)) ? (X) : (Y))

#include <libs/DSP/Include/dsp/basic_math_functions.h>



#ifdef OVERRIDE_MDF_INNER_PROD
static spx_word32_t mdf_inner_prod(const spx_word16_t * x, const spx_word16_t * y, int len)
{
    spx_word32_t    sum;
	// float vector dot product
    sum = vecdotf(x, y, len);
    return sum;
}
#endif

#ifdef OVERRIDE_MDF_POWER_SPECTRUM
static void power_spectrum(const spx_word16_t * X, spx_word32_t * ps, int N)
{
	//Compute power spectrum of a complex vector
    ps[0] = MULT16_16(X[0], X[0]);
	// X[i]*X[i] + X[i+1]*X[i+1]
    adi_sharcfx_cmplx_mag_squared_f32(&X[1], ps + 1, N - 1);
}
#endif



#ifdef OVERRIDE_MDF_VEC_SUB
static void vect_sub(const spx_word16_t * pSrcA, const spx_word16_t * pSrcB, spx_word16_t * pDst, uint32_t blockSize)
{
	// float vector subtraction
	adi_sharcfx_sub_f32(pSrcA, pSrcB, pDst, blockSize);
}

#endif

#ifdef OVERRIDE_MDF_VEC_ADD
static void vect_add(const spx_word16_t * pSrcA, const spx_word16_t * pSrcB, spx_word16_t * pDst, uint32_t blockSize)
{
	// float vector addition
    adi_sharcfx_add_f32(pSrcA, pSrcB, pDst, blockSize);
}

#endif

#ifdef OVERRIDE_MDF_VEC_MULT
static void vect_mult(const spx_word16_t * pSrcA, const spx_word16_t * pSrcB, spx_word16_t * pDst, uint32_t blockSize)
{
	// float vector multiplication
    adi_sharcfx_mult_f32(pSrcA, pSrcB, pDst, blockSize);
}
#endif


#ifdef OVERRIDE_MDF_VEC_SCALE
static void vect_scale(const spx_word16_t * pSrc, spx_word16_t scale, spx_word16_t * pDst, uint32_t blockSize)
{
	// float vector scaling
    adi_sharcfx_scale_f32(pSrc, scale, pDst, blockSize);
}
#endif


#ifdef OVERRIDE_MDF_SMOOTHED_ADD
static void smoothed_add(const spx_word16_t * pSrc1, const spx_word16_t * pWin1,
                  const spx_word16_t * pSrc2, const spx_word16_t * pWin2, spx_word16_t * pDst, uint16_t frame_size, uint16_t nbChan, uint16_t N)

{
	//Blend error and echo residual to apply a smooth transition to avoid introducing blocking artifacts.
	for (int chan = 0; chan < nbChan; chan++) 
	{
		xb_vecMxf32 Src1, Src2, Win1, Win2;
		xb_vecMxf32 out, out1, out2;

		xb_vecMxf32 *restrict ppSrc1 = (xb_vecMxf32 *)(pSrc1 + chan*N);
		xb_vecMxf32 *restrict ppSrc2 = (xb_vecMxf32 *)(pSrc2 + chan*N);
		xb_vecMxf32 *restrict ppWin1 = (xb_vecMxf32 *)(pWin1);
		xb_vecMxf32 *restrict ppWin2 = (xb_vecMxf32 *)(pWin2);
		valign ppSrc1a = PDX_LA_MXF32_PP(ppSrc1);
		valign ppSrc2a = PDX_LA_MXF32_PP(ppSrc2);
		valign ppWin1a = PDX_LA_MXF32_PP(ppWin1);
		valign ppWin2a = PDX_LA_MXF32_PP(ppWin2);

		xb_vecMxf32* __restrict pvOut  = (xb_vecMxf32 *)(pDst + chan*N);
		valign pvOuta = PDX_LA_MXF32_PP(pvOut);
		//better performance with loop unrolling : 2 loops unrolled
		for (int count = 0; count < frame_size; count+=16)
		{
			PDX_LA_MXF32_IP(Src1, ppSrc1a, ppSrc1);
			PDX_LA_MXF32_IP(Src2, ppSrc2a, ppSrc2);
			PDX_LA_MXF32_IP(Win1, ppWin1a, ppWin1);
			PDX_LA_MXF32_IP(Win2, ppWin2a, ppWin2);
			out = PDX_ADD_MXF32(PDX_MUL_MXF32(Src1, Win1), PDX_MUL_MXF32(Src2, Win2)); //PDX_MULA_MXF32(out, Src2,Win2);
			PDX_SAV_MXF32_XP (out, pvOuta, pvOut, 32);

			PDX_LA_MXF32_IP(Src1, ppSrc1a, ppSrc1);
			PDX_LA_MXF32_IP(Src2, ppSrc2a, ppSrc2);
			PDX_LA_MXF32_IP(Win1, ppWin1a, ppWin1);
			PDX_LA_MXF32_IP(Win2, ppWin2a, ppWin2);
			out = PDX_ADD_MXF32(PDX_MUL_MXF32(Src1, Win1), PDX_MUL_MXF32(Src2, Win2));
			PDX_SAV_MXF32_XP (out, pvOuta, pvOut, 32);
		}
		PDX_SAPOS_MXF32_FP(pvOuta,pvOut);
	}
}
#endif


#ifdef OVERRIDE_MDF_POWER_SPECTRUM_ACCUM
static inline void power_spectrum_accum(const spx_word16_t *X, spx_word32_t *ps, int M)
{
	//Compute power spectrum of a complex vector and accumulate the result
	ps[0] += X[0]* X[0];
	int n;
	int N = M/2;
	xb_vecMxf32 Xre, Xim, Zout, temp, Zin;
	xb_vecMx32 temp0, exp_re, exp_im, exp_abs;
	vboolM b_cond;
	valign X_va, Z_va;

	xb_vecMxf32 *restrict x = (xb_vecMxf32*)(X+1);
	xb_vecMxf32 *z = (      xb_vecMxf32 *)(ps+1);
	X_va = PDX_LA_MXF32_PP(x);
	Z_va = PDX_LA_MXF32_PP(z);
	xb_vecMxf32 * Zp = (  xb_vecMxf32 *)(ps+1);
	valign Za =PDX_LA_MXF32_PP(Zp);
	//ps[i] += Xre[i]^2 + Xim[i]^2
	for ( n = 0; n < M-1; n+=16 )
	{
		int byte = MIN(32, 4*(M/2 - n/2 ));
		int byte1 = MIN(32, 4*(M-1-n));
		int byte2 = MIN(32, 4*(M-1-n - 8));
		PDX_LAV_MXF32_XP( Xre, X_va, x, byte1 );
		PDX_LAV_MXF32_XP( Xim, X_va, x, byte2  );
		PDX_LAV_MXF32_XP(Zout, Za, Zp, byte );
		PDX_DSELI_MXF32(Xim, Xre, Xim, Xre, PDX_DSELI_32B_DEINTERLEAVE_1); // separating out real and imaginary values and saving it in Xim and Xre vectors. 
		PDX_MULA_MXF32(Zout, Xre, Xre); // Xre*Xre
		PDX_MULA_MXF32(Zout, Xim, Xim); //Xim*Xim
		PDX_SAV_MXF32_XP( Zout, Z_va, z, byte);
	}
	PDX_SAPOS_MXF32_FP( Z_va, z );
}
#endif


#ifdef OVERRIDE_MDF_PREEMPH_FLT
static inline int32_t mdf_preemph(spx_word16_t * input, spx_word16_t * out, spx_word16_t preemph, int frame_size, spx_word16_t * mem)
{
	//Copy spx_word16_(int16_t) input data to buffer and apply pre-emphasis filter.
	spx_word32_t tmp32;
	out[0] = input[0] -  preemph*mem[0];
	int32_t saturated = 0;
	xb_vecMxf32 *__restrict Xip = ( xb_vecMxf32 *)(input);
	xb_vecMxf32 *__restrict Xp = ( xb_vecMxf32 *)(input+1);
	xb_vecMxf32 *Zp = (xb_vecMxf32 *)(out+1);
	valign Xia = PDX_LA_MXF32_PP(Xip);
	valign Xa = PDX_LA_MXF32_PP(Xp);
	valign Za = PDX_LA_MXF32_PP(Zp);
	xb_vecMxf32 Vpreemph = -preemph;
	xb_vecMxf32 Xi, X, Z;
	mem[0] = input[frame_size-1];

	int loopend = (frame_size/8)*8;
	for (int i=0; i<frame_size;i+=8)
	{
		PDX_LA_MXF32_IP( Xi, Xia, Xip);// mem = input
		PDX_LA_MXF32_IP( X, Xa, Xp); //input+1
		Z = PDX_ADD_MXF32(X,Vpreemph*Xi );//Z = X + Vpreemph*Xi;
		PDX_SA_MXF32_IP( Z, Za, Zp);
	}
	int bytes = MIN(32, 4*(frame_size-loopend));
	PDX_LAV_MXF32_XP( Xi, Xia, Xip, bytes);
	PDX_LAV_MXF32_XP( X, Xa, Xp,  (bytes-4) );
	Z = PDX_ADD_MXF32(X,Vpreemph*Xi );//Z = X + Vpreemph*Xi;
	PDX_SAV_MXF32_XP( Z, Za, Zp, bytes );

	PDX_SAPOS_MXF32_FP( Za, Zp);

	return saturated;
}
#endif


#ifdef OVERRIDE_MDF_SMOOTH_FE_NRG
static inline void smooth_fe_nrg(spx_word32_t * in1, spx_word16_t c1, spx_word32_t * in2, spx_word16_t c2, spx_word32_t * pDst, uint16_t frame_size)
{
	//Smooth far end energy estimate over time.
	xb_vecMxf32 *restrict X1p = ( xb_vecMxf32 *)(in1);
	xb_vecMxf32 *restrict X2p = ( xb_vecMxf32 *)(in2);
	xb_vecMxf32 *restrict Zp = (xb_vecMxf32 *)(pDst);
	valign X1a = PDX_LA_MXF32_PP(X1p);
	valign X2a = PDX_LA_MXF32_PP(X2p);
	valign Za = PDX_LA_MXF32_PP(Zp);
	xb_vecMxf32 vc1 = c1;
	xb_vecMxf32 vc2 = c2;
	xb_vecMxf32 one = (xb_vecMxf32)1;
	xb_vecMxf32 X1, X2, Z, out, out2, X3, X4;
	
	//pDst[i] = c1*in1[i] + 1 +c2*in2[i];
	for (int i=0; i<frame_size;i+=16)
	{
		// unrolled the loop to process two 8way-32bit vectors
		PDX_LA_MXF32_IP( X1, X1a, X1p );
		PDX_LA_MXF32_IP( X2, X2a, X2p );
		PDX_LA_MXF32_IP( X3, X1a, X1p );
		PDX_LA_MXF32_IP( X4, X2a, X2p );
		out = one; out2 = one;
		PDX_MULA_MXF32(out, X1, vc1);
		PDX_MULA_MXF32(out, X2, vc2);
		PDX_MULA_MXF32(out2, X3, vc1);
		PDX_MULA_MXF32(out2, X4, vc2);
		int bytes =  MIN(32, 4*(frame_size-i));
		int bytes2 =  MIN(32, 4*(frame_size-i - 8)); //the left out elements
		PDX_SAV_MXF32_XP( out, Za, Zp,bytes );
		PDX_SAV_MXF32_XP( out2, Za, Zp,bytes2 );
	}
	PDX_SAPOS_MXF32_FP( Za, Zp);
}
#endif

#ifdef OVERRIDE_MDF_SPECTRAL_MUL_ACCUM
#include <vector.h>
static inline void spectral_mul_accum(const spx_word16_t *X, const spx_word32_t *Y, spx_word16_t *acc, int N, int M)
{
	//Compute cross-power spectrum of a complex vectors and accumulate.
	for (int i=0;i<N;i++)
	 	acc[i] = 0; //setting to zero first

	for (int j=0;j<M;j++)
	{
		int n, byte;
		xb_vecMxf32 X0, Y0, Z0, Z1, Zout;
		acc[0] += X[0 + j*N]*Y[0 + j*N];
		const xb_vecMxf32 * __restrict px = (const xb_vecMxf32 *)(X+1 + j*N);
		const xb_vecMxf32 * __restrict py = (const xb_vecMxf32 *)(Y+1 + j*N);
		valign ax = PDX_LA_MXF32_PP(px);
		valign ay = PDX_LA_MXF32_PP(py);

		xb_vecMxf32* pz = (xb_vecMxf32 *)(acc+1);
		valign az = PDX_LA_MXF32_PP(pz);
		xb_vecMxf32* psz = (xb_vecMxf32 *)(acc+1);
		valign asz = PDX_LA_MXF32_PP(pz);

		for (n = 0; n<N-2; n +=8)
		{
			int byte = MIN(8*4, (N-2-n)*sizeof(*X) );
			PDX_LAV_MXF32_XP(X0, ax, px, byte);
			PDX_LAV_MXF32_XP(Y0, ay, py, byte);
			PDX_LAV_MXF32_XP(Zout, az, pz, byte);
			Z0 = PDX_MULM_MXF32( X0, Y0, 4);   //00 01 00       x.re*y.re     x.re*y.im
			Z1 = PDX_MULM_MXF32( X0, Y0, 27);  //01 10 11      -x.im*y.im     x.im*y.re
			Zout += PDX_ADD_MXF32(Z0, Z1); //Zout += (Z0 + Z1);
			PDX_SAV_MXF32_XP(Zout, asz, psz, byte);
		}
		PDX_SAPOS_MXF32_FP(asz, psz);
		acc[N-1] += X[N-1 + j*N]*Y[N-1 + j*N];
	}
}
#endif


#ifdef OVERRIDE_MDF_DEEMPH
static int mdf_deemph(const spx_int16_t * restrict micin, spx_word16_t * restrict input,
        spx_word16_t * restrict e, spx_int16_t * restrict out,
        spx_word16_t preemph, int frame_size, spx_word16_t * restrict mem, int stride)
{
	//Compute error signal, check input saturation and convert / saturate strided output to spx_int16_t (int16_t).
    int32_t saturated = 0; spx_word32_t tmp_out;
    int i;
    spx_word16_t preemph_mem = preemph * mem[0];
    for (i = 0; i < frame_size; i++) 
	{
        tmp_out = input[i] - e[i];
        tmp_out += preemph_mem;
        out[i * stride] = WORD2INT(tmp_out);
        preemph_mem = preemph * tmp_out;
        if (micin[i * stride] <= -32000 || micin[i * stride] >= 32000) 
		{
			if (saturated == 0)
				saturated = 1;
        }
    }
    *mem = tmp_out;
    return saturated;
}
#endif


#ifdef OVERRIDE_MDF_VEC_CLEAR // working fine - can call fill function too, but need to include header
static inline void vect_clear(spx_word16_t * pDst, uint32_t blockSize)
{
	xb_vecMxf32 In = 0;
	xb_vecMxf32 *restrict pvOut = (xb_vecMxf32 *)(pDst);
	valign pvOuta = PDX_LA_MXF32_PP(pvOut);

    for (int i = 0; i < blockSize; i += 8)
    {
    	int byte = (blockSize-i)*sizeof(*pDst);
    	PDX_SAV_MXF32_XP(In, pvOuta, pvOut, byte); // clearing vector by filling it with zeroes
    }
   PDX_SAPOS_MXF32_FP(pvOuta,pvOut);
}
#endif


#ifdef OVERRIDE_MDF_VEC_SUB_INT16 // Optimized - no change
static void vect_sub16(const spx_int16_t * pSrcA, const spx_int16_t * pSrcB, spx_word16_t * pDst, uint32_t blockSize)
{
    xb_vec2Mx16 *restrict pA = (xb_vec2Mx16 *)(pSrcA);
    xb_vec2Mx16 *restrict pB = (xb_vec2Mx16 *)(pSrcB);
    valign pAa = PDX_LA_2MX16_PP(pA);
    valign pBa = PDX_LA_2MX16_PP(pB);

    xb_vec2Mx16 *restrict pvOut = (xb_vec2Mx16 *)(pDst);
    valign pvOuta = PDX_LA_2MX16_PP(pvOut);
    xb_vec2Mx16 A, B, out;

    for (int i = 0; i < blockSize; i+=16)
    {
    	int byte = (blockSize-i)*sizeof(*pDst);
		PDX_LAV_2MX16_XP(A, pAa, pA, byte); // loading 16way 16bit
		PDX_LAV_2MX16_XP(B, pBa, pB, byte); // loading 16way 16bit
    	out = PDX_SUBS_2MX16(A, B);
    	PDX_SAV_2MX16_XP(out, pvOuta, pvOut, byte);
    }
    PDX_SAPOS_2MX16_FP(pvOuta,pvOut);
}
#endif


#ifdef OVERRIDE_MDF_WEIGHT_SPECT_MUL_CONJ
#include "pseudofloat.h"
static void weighted_spectral_mul_conj(const spx_float_t * w, const spx_float_t p, const spx_word16_t * X, const spx_word16_t * Y, spx_word32_t * prod, int N)
{
	//Compute weighted cross-power spectrum of a complex vector with conjugate. 
	spx_word16_t weigh[512];
	int j =0;
	for (int i = 0;i<N-1; i++)
	{
		if(i%2 ==0) j++;
		weigh[i] = p*w[j]; //replicating w[0-N/2] to weigh[0-N] for multiplying it with real and imag parts
	}

	int n, byte;
	xb_vecMxf32 X0, Y0, Z0, W0, Z1, Zout, Zo;
	prod[0] =  p*w[0]*X[0]*Y[0];
	const xb_vecMxf32 *restrict px = (const xb_vecMxf32 *)(X+1);
	const xb_vecMxf32 *restrict py = (const xb_vecMxf32 *)(Y+1);
	const xb_vecMxf32 *restrict pw = (const xb_vecMxf32 *)(weigh);
	valign ax = PDX_LA_MXF32_PP(px);
	valign ay = PDX_LA_MXF32_PP(py);
	valign aw = PDX_LA_MXF32_PP(pw);

	xb_vecMxf32 *pz = (xb_vecMxf32 *)(prod+1);
	valign az = PDX_LA_MXF32_PP(pz);

	for (n = 0; n<N-2; n +=8)
	{
		int byte = MIN(8*4, (N-2-n)*sizeof(*X) );
		PDX_LAV_MXF32_XP(X0, ax, px, byte);
		PDX_LAV_MXF32_XP(Y0, ay, py, byte);
		PDX_LAV_MXF32_XP(W0, aw, pw, byte);
		Z0 = PDX_MULM_MXF32( X0, Y0, 4);   //00 01 00       x.re*y.re     x.re*y.im
		Z1 = PDX_MULM_MXF32( X0, Y0, 27);  //01 10 11      -x.im*y.im     x.im*y.re
		Zo = PDX_SUB_MXF32(Z0, Z1);
		Zout = PDX_MUL_MXF32(W0, Zo); //multiplying weights with the result vector
		PDX_SAV_MXF32_XP(Zout, az, pz, byte);
	}
	PDX_SAPOS_MXF32_FP(az, pz);
	prod[N-1] =  p*w[N/2]*X[N-1]*Y[N-1];

}
#endif


#ifdef OVERRIDE_MDF_VEC_COPY //giving error- therefore code not included
static void vect_copy(spx_word16_t *dst, spx_word16_t *src, uint32_t blockSize)
{
	#if defined (SHARC_FXx)
		//memcpy(dst, src, blockSize);
		xb_vecMxf32 *psrc = (xb_vecMxf32 *)(src);
		valign psrca = PDX_LA_MXF32_PP(psrc);
		xb_vecMxf32 *pvOut = (xb_vecMxf32 *)(dst);
		valign pvOuta = PDX_LA_MXF32_PP(pvOut);
		xb_vecMxf32 In;
		for (int i = 0; i < blockSize; i += 8)
		{
			int byte = (blockSize-i)*sizeof(*dst);
			PDX_LAV_MXF32_XP(In, psrca, psrc, byte); //copy
			PDX_SAV_MXF32_XP(In, pvOuta, pvOut, byte); //paste
		}
		PDX_SAPOS_MXF32_FP(pvOuta,pvOut);
	#else
		for (int i = 0; i < blockSize; i++)
		{
			dst[i] = src[i];
		}
	#endif
}
#endif


#ifdef OVERRIDE_MDF_FILTERED_SPEC_AD_XCORR
static void filtered_spectra_cross_corr(spx_word32_t * pRf, spx_word32_t * pEh, spx_word32_t * pYf, spx_word32_t * pYh,
                                        spx_float_t * Pey, spx_float_t * Pyy, spx_word16_t spec_average, uint16_t frame_size)
{
	// Compute filtered spectra and (cross-)correlations.
	xb_vecMxf32 E, Y, Rf, Yf, Eh, Yh;
	xb_vecMxf32 spec1 = spec_average;  xb_vecMxf32 spec2 = 1 -spec_average;
	xb_vecMxf32 *ppRf = (xb_vecMxf32 *)(pRf);
	xb_vecMxf32 *ppYf = (xb_vecMxf32 *)(pYf);
	valign pRfa = PDX_LA_MXF32_PP(ppRf);
	valign pYfa = PDX_LA_MXF32_PP(ppYf);

	xb_vecMxf32 *pE = (xb_vecMxf32 *)(pEh);
	xb_vecMxf32 *pY = (xb_vecMxf32 *)(pYh);
	valign pEa = PDX_LA_MXF32_PP(pE);
	valign pYa = PDX_LA_MXF32_PP(pY);

	xb_vecMxf32 *pE2 = (xb_vecMxf32 *)(pEh);
	xb_vecMxf32 *pY2 = (xb_vecMxf32 *)(pYh);
	valign pE2a = PDX_LA_MXF32_PP(pE);
	valign pY2a = PDX_LA_MXF32_PP(pY);

	// pEh[j] = spec2 * pEh[j] + spec1 * pRf[j];
	// pYh[j] = spec2 * pYh[j] + spec1 * pYf[j];
	xb_vecMxf32 accumE = 0; xb_vecMxf32 accumY = 0;;
	for (int n = 0; n<frame_size+1 ; n +=8)
		{
			int byte = MIN(8*4, (frame_size+1-n)*sizeof(*pEh) );
			PDX_LAV_MXF32_XP(E, pEa, pE, byte);
			PDX_LAV_MXF32_XP(Y, pYa, pY, byte);
			PDX_LAV_MXF32_XP(Rf, pRfa, ppRf, byte);
			PDX_LAV_MXF32_XP(Yf, pYfa, ppYf, byte);
			Eh = (Rf - E); Yh = (Yf - Y);
			PDX_MULA_MXF32(accumE, Eh, Yh);
			PDX_MULA_MXF32(accumY, Yh, Yh);
			xb_vecMxf32 Ze = 0, Zy = 0;
			PDX_MULA_MXF32(Ze, spec1, Rf );
			PDX_MULA_MXF32(Ze, spec2, E );
			PDX_MULA_MXF32(Zy, spec1, Yf );
			PDX_MULA_MXF32(Zy, spec2, Y );
			PDX_SAV_MXF32_XP(Ze, pE2a, pE2, byte);
			PDX_SAV_MXF32_XP(Zy, pY2a, pY2, byte);
		}
	PDX_SAPOS_MXF32_FP(pE2a, pE2);
	PDX_SAPOS_MXF32_FP(pY2a, pY2);
	*Pey+= PDX_RADD_MXF32(accumE);
	*Pyy+= PDX_RADD_MXF32(accumY);
}
#endif

#ifdef OVERRIDE_MDF_ADJUST_PROP
// Adjust properties function
static inline void mdf_adjust_prop(const spx_word32_t *restrict W, int N, int M, int P, spx_word16_t *prop) {
	//Computes filter adaptation rate, proportional to inverse of weight filter energy.
    xb_vecMxf32 W0, W1, W2, W3, acc1, acc2, acc3, acc0, vprop;
    xb_vecMxf32 W4, W5, W6, W7, acc4, acc5, acc6, acc7, vprop2;
	int i;
	spx_word16_t max_sum = 1;
	xb_vecMxf32 vmax_sum = max_sum;

	xb_vecMxf32 *pz = (xb_vecMxf32 *)(prop);
	valign az = PDX_LA_MXF32_PP(pz);
	//tmp = 1;
	// tmp +=  W[p * N * M + i * N + j]^2;
	// prop[i] = sqrt(tmp);
	for (i = 0; i < M; i+=8) {
		acc0 = 0; acc1 = 0; acc2 = 0; acc3 = 0;
		acc4 = 0; acc5 = 0; acc6 = 0; acc7 = 0;
		spx_word32_t tmp = 1;
		for (int p = 0; p < P; p++)
		{
			const xb_vecMxf32 *restrict pw = (const xb_vecMxf32 *)(W + p * N * M + i * N);
			valign aw = PDX_LA_MXF32_PP(pw);
			const xb_vecMxf32 *restrict pw1 = (const xb_vecMxf32 *)(W + p * N * M + (i+1) * N);
			valign aw1 = PDX_LA_MXF32_PP(pw1);
			const xb_vecMxf32 *restrict pw2 = (const xb_vecMxf32 *)(W + p * N * M + (i+2) * N);
			valign aw2 = PDX_LA_MXF32_PP(pw2);
			const xb_vecMxf32 *restrict pw3 = (const xb_vecMxf32 *)(W + p * N * M + (i+3) * N);
			valign aw3 = PDX_LA_MXF32_PP(pw3);

			const xb_vecMxf32 *restrict pw4 = (const xb_vecMxf32 *)(W + p * N * M + (i+4) * N);
			valign aw4 = PDX_LA_MXF32_PP(pw4);
			const xb_vecMxf32 *restrict pw5 = (const xb_vecMxf32 *)(W + p * N * M + (i+5) * N);
			valign aw5 = PDX_LA_MXF32_PP(pw5);
			const xb_vecMxf32 *restrict pw6 = (const xb_vecMxf32 *)(W + p * N * M + (i+6) * N);
			valign aw6 = PDX_LA_MXF32_PP(pw6);
			const xb_vecMxf32 *restrict pw7 = (const xb_vecMxf32 *)(W + p * N * M + (i+7) * N);
			valign aw7 = PDX_LA_MXF32_PP(pw7);
			for (int j = 0; j<N; j +=8)
			{
				int byte = MIN(8*4, (N-j)*4 );
				PDX_LAV_MXF32_XP(W0, aw, pw, byte);
				PDX_LAV_MXF32_XP(W1, aw1, pw1, byte);
				PDX_LAV_MXF32_XP(W2, aw2, pw2, byte);
				PDX_LAV_MXF32_XP(W3, aw3, pw3, byte);
				PDX_MULAN_MXF32(acc0, W0, W0);
				PDX_MULAN_MXF32(acc1, W1, W1);
				PDX_MULAN_MXF32(acc2, W2, W2);
				PDX_MULAN_MXF32(acc3, W3, W3);

				PDX_LAV_MXF32_XP(W4, aw4, pw4, byte);
				PDX_LAV_MXF32_XP(W5, aw5, pw5, byte);
				PDX_LAV_MXF32_XP(W6, aw6, pw6, byte);
				PDX_LAV_MXF32_XP(W7, aw7, pw7, byte);
				PDX_MULAN_MXF32(acc4, W4, W4);
				PDX_MULAN_MXF32(acc5, W5, W5);
				PDX_MULAN_MXF32(acc6, W6, W6);
				PDX_MULAN_MXF32(acc7, W7, W7);

			}
		}
		prop[i]   = sqrtf( 1 + PDX_RADD_MXF32(acc0) ); // 1 added, because tmp is initialized with 1
		prop[i+1]   = sqrtf( 1 + PDX_RADD_MXF32(acc1) );//RADD adds all the elements in the vector register
		prop[i+2]   = sqrtf( 1 + PDX_RADD_MXF32(acc2) );
		prop[i+3]   = sqrtf( 1 + PDX_RADD_MXF32(acc3) );

		prop[i+4]   = sqrtf( 1 + PDX_RADD_MXF32(acc4) );
		prop[i+5]   = sqrtf( 1 + PDX_RADD_MXF32(acc5) );
		prop[i+6]   = sqrtf( 1 + PDX_RADD_MXF32(acc6) );
		prop[i+7]   = sqrtf( 1 + PDX_RADD_MXF32(acc7) );

		xb_vecMxf32 vprop;
		xb_vecMxf32 *pprop = (xb_vecMxf32 *)(prop);
		valign propa = PDX_LA_MXF32_PP(pprop);
		PDX_LAV_MXF32_XP(vprop, propa, pprop, 8); //loading the prop value calculated above
		vmax_sum = PDX_MAX_MXF32(vprop, vmax_sum); 
	}
	max_sum = PDX_RMAX_MXF32(vmax_sum);

	float prop_sum = 1;
	for (i = 0; i < M; i++) 
	{
		prop[i] += 0.1f * max_sum;
		prop_sum += prop[i];
	}
	for (i = 0; i < M; i++) 
	{
		prop[i] = 0.99f * prop[i] / prop_sum;
	}
}
#endif


#ifdef OVERRIDE_MDF_NORM_LEARN_RATE_CALC
static void mdf_nominal_learning_rate_calc(spx_word32_t *restrict pRf, spx_word32_t *restrict power,
                                           spx_word32_t *restrict pYf, spx_float_t *restrict power_1, spx_word16_t leak_estimate, spx_word16_t RER, uint16_t frame_size)
{
//Normal learning rate calculation once we're past the minimal adaptation phase.
//giving error currently, will modify later
#if defined (SHARC_FXx)
    xb_vecMxf32 Yf,Rf, R, E, limit, Vpower, out1, out;
    xb_vecMxf32 one = 1; xb_vecMxf32 half = 0.5;
    xb_vecMxf32 seven = 0.7; xb_vecMxf32 three = 0.3*RER;  xb_vecMxf32 ten = 10;
    xb_vecMxf32 Vleak_estimate = leak_estimate;
	xb_vecMxf32 *ppYf = (xb_vecMxf32 *)(pYf);
	valign pYfa = PDX_LA_MXF32_PP(ppYf);
	xb_vecMxf32 *ppRf = (xb_vecMxf32 *)(pRf);
	valign pRfa = PDX_LA_MXF32_PP(ppRf);

	xb_vecMxf32 *pPower = (xb_vecMxf32 *)(power);
	valign pPowera = PDX_LA_MXF32_PP(pPower);
	xb_vecMxf32 *pPower_1 = (xb_vecMxf32 *)(power_1);
	valign pPower_1a = PDX_LA_MXF32_PP(pPower_1);

    for (int i = 0; i < frame_size; i+=8) {
        float    r, e;
        int bytes = MIN(8*4, 4*(frame_size-i));
        PDX_LAV_MXF32_XP(Yf, pYfa,ppYf, bytes );
        PDX_LAV_MXF32_XP(Rf, pRfa,ppRf, bytes );
        PDX_LAV_MXF32_XP(Vpower, pPowera, pPower, bytes );
        R  = Vleak_estimate*Yf;
        E = Rf + one;
        limit = half*E;
        PDX_MAX_MXF32(R, limit);
        R = seven*R +  three*E;
        out1 = E*( Vpower + ten );
        out = PDX_DIV_MXF32(R, out1);
        PDX_SAV_MXF32_XP( out, pPower_1a, pPower_1, bytes);
    }
    PDX_SAPOS_MXF32_FP(pPower_1a, pPower_1);

#else

    for (int i = 0; i < frame_size; i++) {
        spx_word32_t    r, e;
        r = leak_estimate*pYf[i];
        e = pRf[i] + 1;
        if (r > .5 * e)
            r = .5 * e;
        r = 0.7*r + 0.3*RER*e;
        power_1[i] = r/(e*(power[i] + 10));
    }

#endif
}

#endif

#ifdef OVERRIDE_MDF_CONVERG_LEARN_RATE_CALC
static void mdf_non_adapt_learning_rate_calc(spx_word32_t * power, spx_float_t * power_1, spx_word16_t adapt_rate, uint16_t frame_size)
{
	//Part of the process of the computing the adaption rate when filter is not yet adapted enough. This routine divides the adaptation rate by Far End power over the whole subframe.
	const xb_vecMxf32 *restrict ppwr = (const xb_vecMxf32 *)(power);
	valign apwr = PDX_LA_MXF32_PP(ppwr);
	xb_vecMxf32 *restrict ppwr_1 = (xb_vecMxf32 *)(power_1);
	valign appwr_1 = PDX_LA_MXF32_PP(ppwr_1);
	xb_vecMxf32 vpwr, zout;
	xb_vecMxf32 vrate = adapt_rate;
	xb_vecMxf32 ten = 10;
	//power_1[i] = (adapt_rate)/(power[i]+10));
	for (int j = 0; j<frame_size; j +=8)
	{
		int byte = MIN(8*4, (frame_size-j)*4 );
		PDX_LAV_MXF32_XP(vpwr, apwr, ppwr, byte);
		zout = PDX_ADD_MXF32(vpwr, ten);
		zout = PDX_DIV_MXF32(vrate,zout);
		PDX_SAV_MXF32_XP(zout, appwr_1, ppwr_1, byte);
	}
	PDX_SAPOS_MXF32_FP(appwr_1, ppwr_1);
}

#endif

#ifdef OVERRIDE_MDF_DC_NOTCH
static void filter_dc_notch16(const spx_int16_t * in, spx_word16_t radius, spx_word16_t * out, int len, spx_mem_t * mem, int stride)
{
//Notch filter with strided spx_int16_t (int16_t) type input and spx_word16_t (floating point) output.
    int             i;
    spx_word16_t    den2;
    float32_t       mem1 = mem[1];
    float32_t       mem0 = mem[0];
    xb_vecMxf32 vradius = radius;
    den2 = radius * radius + .7f * (1.0f - radius) * (1.0f - radius);


    xb_vecMxf32 * pout = (xb_vecMxf32 *)(out);
    valign pouta = PDX_LA_MXF32_PP(pout);
	spx_word32_t vecout[300];
    xb_vecMxf32 * pvecout = (xb_vecMxf32 *)(vecout);
    valign pvecouta = PDX_LA_MXF32_PP(pvecout);
    spx_word32_t    vout; xb_vecMxf32 vec;
    for (i = 0; i < len ; i+=8) {
        xb_vecMxf32     voutF;
		//processing 8 values in the for loop to get it vectorized
        for (int j = 0; j < 8; j++) {
        	vout = mem0 + in[i + j];
            vecout[i + j] = vout;
            mem0 = mem1 +  2 * (-in[i + j] + radius * vout);
            mem1 = in[i + j] - den2 * vout;
        }
        int byte = MIN(32, 4*(len-i));
		PDX_LAV_MXF32_XP(vec, pvecouta, pvecout, byte);
		vec = vec*vradius;
		PDX_SAV_MXF32_XP(vec, pouta, pout, byte);
    }

    PDX_SAPOS_MXF32_FP(pouta, pout);

    mem[1] = mem1;
    mem[0] = mem0;
}

#endif

#else

/* FIXED_POINT not needed for EEMBC AudioMark */
#error "Fixed Point Optimization is not available"

#endif //defined(FLOATING_POINT) end


