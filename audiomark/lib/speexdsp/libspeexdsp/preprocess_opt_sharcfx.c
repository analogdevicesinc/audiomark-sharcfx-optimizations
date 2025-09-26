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



#include <math.h>
#include <float.h>
#include <complex.h>
#include <vector.h>

#ifdef __XTENSA__
#include <xtensa/sim.h>
#include <xtensa/tie/xt_pdxn.h>
#endif

#include <stdint.h>
#include "arch.h"
#include "preprocess_opt_sharcfx.h"

#define MIN(X, Y) (((X) < (Y)) ? (X) : (Y))
#define MAX(X, Y) (((X) > (Y)) ? (X) : (Y))

#ifdef OVERRIDE_ANR_VEC_MUL
static void vect_mult(const spx_word16_t * pSrcA, const spx_word16_t * pSrcB, spx_word16_t * pDst, uint32_t blockSize)
{
	//multiplication of two float vectors pSrcA, pSrcB of size blockSize and returning the result in pDst vector
	vecvmltf(pSrcA, pSrcB, pDst, blockSize);
}
#endif


#ifdef OVERRIDE_ANR_VEC_CONV_FROM_INT16
static void vect_conv_from_int16(const spx_int16_t * pSrc, spx_word16_t * pDst, uint32_t blockSize)
{
	//Convert spx_int16_t (16-bit signed integer) vector to spx_word16_t (floating-point).
	const xb_vecMx16 * __restrict psrc = (const xb_vecMx16 *)(pSrc); // source 8x16bit vector
	xb_vecMxf32 *pdst = ( xb_vecMxf32 *)(pDst); // destination 8x32float vector
	valign adst = PDX_LA_MXF32_PP(pdst);
	valign asrc = PDX_LA_MX16_PP(psrc);
	xb_vecMx32 vA;
	for (int i = 0; i < blockSize; i+=8)
	{
		int byte = MIN(16, 2*(blockSize-i));
		PDX_LA32_MX16_XP(vA, asrc, psrc, byte);
		xb_vecMxf32 vout = (xb_vecMxf32)vA; // converting 16bit int to 32bit float
		PDX_SAV_MXF32_XP(vout, adst, pdst, 2*byte);
	}
	PDX_SAPOS_MXF32_FP( adst, pdst);
}
#endif


#ifdef OVERRIDE_ANR_OLA
static void vect_ola(const spx_word16_t * pSrcA, const spx_word16_t * pSrcB, spx_int16_t * pDst, uint32_t blockSize)
{
	//spx_word16_t vector overlap and add
	const xb_vecMxf32 * __restrict psrca = (const xb_vecMxf32 *)(pSrcA);
	const xb_vecMxf32 * __restrict psrcb = (const xb_vecMxf32 *)(pSrcB);
	xb_vecMx16 *pdst = ( xb_vecMx16 *)(pDst);
	valign asrcb = PDX_LA_MXF32_PP(psrcb);
	valign asrca = PDX_LA_MXF32_PP(psrca);
	valign adst = PDX_LA_MX16_PP(pdst);
	xb_vecMxf32 vA, vB, vout;

    for (int i = 0; i < blockSize; i+=8)
    {
		int byte = MIN(32, 4*(blockSize-i));
		PDX_LAV_MXF32_XP(vA, asrca, psrca, byte); // loading 8-32bit float values into vA vector
		PDX_LAV_MXF32_XP(vB, asrcb, psrcb, byte); // loading 8-32bit float values into vA vector
		vout = vA + vB; // adding
		vout = PDX_FIRINT_MXF32(vout); //  // floating point round with mode to integral
		xb_vecMx32 intout = (xb_vecMx32)vout; // converting float to int
		PDX_SAV32_MX16_XP(intout, adst, pdst, byte); // saving 32bit int as 16bit int
    }
    PDX_SAPOS_MX16_FP( adst, pdst);
}
#endif


#ifdef OVERRIDE_ANR_VEC_COPY //less opt
static void vect_copy(void *dst, const void *src, uint32_t blockSize)
{
	//Generic vector copy.
    memcpy(dst, src, blockSize);
}
#endif


#ifdef OVERRIDE_ANR_COMPUTE_GAIN_FLOOR //no update
static void compute_gain_floor(int noise_suppress, int effective_echo_suppress, spx_word32_t * noise, spx_word32_t * echo, spx_word16_t * gain_floor, int len)
{
	// Compute the gain floor based on different floors for the background noise and residual echo
	//gain_floor[i] = 1.f*sqrt(noise_floor*noise[i]+echo_floor*echo[i]) / sqrt(1+noise[i]+echo[i]);
	xb_vecMxf32 vnoise_floor, vecho_floor;
	vnoise_floor = expf(.2302585f * noise_suppress);
	vecho_floor = expf(.2302585f * effective_echo_suppress);

	const xb_vecMxf32 * __restrict pn = (const xb_vecMxf32 *)(noise);
	const xb_vecMxf32 * __restrict pe = (const xb_vecMxf32 *)(echo);
	valign an = PDX_LA_MXF32_PP(pn);
	valign ae = PDX_LA_MXF32_PP(pe);

	xb_vecMxf32* pgf = (xb_vecMxf32 *)(gain_floor);
	valign agf = PDX_LA_MXF32_PP(pgf);

	xb_vecMxf32 Xn, Xe,nom,denom, Zout;
	xb_vecMxf32 onef = 1.f; xb_vecMxf32 one = 1;
	int loopend = (len/8)*8;
	//computation loop till multiple of 8, vector computation
	for (int n = 0; n<loopend; n+=8)
	{
		nom = 0;
		PDX_LA_MXF32_IP(Xn, an, pn); 
		PDX_LA_MXF32_IP(Xe, ae, pe); 
		PDX_MULAN_MXF32 (nom, vnoise_floor, Xn);
		PDX_MULAN_MXF32 (nom, vecho_floor, Xe);
		denom = PDX_ADD_MXF32(Xe , Xn);
		denom = PDX_ADD_MXF32(denom , one);
		nom = PDX_SQRT_MXF32(nom); //numerator for gain_floor computation
		denom = PDX_SQRT_MXF32(denom); // denominator for gain_floor computation
		Zout = PDX_DIV_MXF32(nom, denom);
		PDX_SA_MXF32_IP(Zout, agf, pgf);
	}
	//computation for last left out values
	int byte = 4*(len-loopend);
	PDX_LAV_MXF32_XP(Xn, an, pn, byte);
	PDX_LAV_MXF32_XP(Xe, ae, pe, byte);
	nom = vnoise_floor * Xn + vecho_floor * Xe;
	nom = onef*PDX_SQRT_MXF32(nom); 
	denom = one + Xn + Xe;
	denom = PDX_SQRT_MXF32(denom);
	Zout = PDX_DIV_MXF32(nom, denom);
	PDX_SAV_MXF32_XP(Zout, agf, pgf, byte);
	PDX_SAPOS_MXF32_FP( agf, pgf);
}
#endif


#ifdef OVERRIDE_ANR_POWER_SPECTRUM //no update
static void power_spectrum(spx_word16_t * ft, spx_word32_t * ps, int N)
{
	//Complex magnitude squared of a spx_word16_t (floating-point) vector.
	ps[0] = ft[0]* ft[0];
	int M = N*2;
	xb_vecMxf32 Xre, Xim, Zout;
	valign X_va, Z_va;

	xb_vecMxf32 *restrict x = (xb_vecMxf32*)(ft+1);
	xb_vecMxf32 *z = (xb_vecMxf32 *)(ps+1);
	X_va = PDX_LA_MXF32_PP(x);
	Z_va = PDX_LA_MXF32_PP(z);

	//ps[i] = ft[2*i-1]*ft[2*i-1] +ft[2*i]*ft[2*i];
	for ( int n = 0; n < M-1; n+=16 )
	{
		Zout = 0;
		int byte = MIN(32, 4*(M/2 - n/2 ));
		int byte1 = MIN(32, 4*(M-1-n));
		int byte2 = MIN(32, 4*(M-1-n - 8));
		PDX_LAV_MXF32_XP( Xre, X_va, x, byte1 );
		PDX_LAV_MXF32_XP( Xim, X_va, x, byte2  );
		PDX_DSELI_MXF32(Xim, Xre, Xim, Xre, PDX_DSELI_32B_DEINTERLEAVE_1); //separating ft[2*i],  ft[2*i-1] in Xre, Xim vector registors
		PDX_MULA_MXF32(Zout, Xre, Xre); //ft[2*i-1]*ft[2*i-1]
		PDX_MULA_MXF32(Zout, Xim, Xim); //ft[2*i]*ft[2*i]
		PDX_SAV_MXF32_XP( Zout, Z_va, z, byte);
	}
	PDX_SAPOS_MXF32_FP( Z_va, z );
}
#endif


#ifdef OVERRIDE_ANR_UPDATE_NOISE_ESTIMATE
static void update_noise_estimate(SpeexPreprocessState2 * st, spx_word16_t beta, spx_word16_t beta_1)
{
	//Update noise estimates
	//  if (!st->update_prob[i] || st->ps[i] < st->noise[i]) then update st->noise[i] = MAX32(0, beta_1*st->noise[i] + beta*st->ps[i] );
    int             N = st->ps_size;
    xb_vecMxf32 *restrict ps = (xb_vecMxf32 *)(st->ps);
    xb_vecMxf32 *noise = (xb_vecMxf32 *)(st->noise);
    xb_vecMx32 *restrict prob = (xb_vecMx32 *)(st->update_prob);
    xb_vecMxf32 vbeta = beta; xb_vecMxf32 vbeta1 = beta_1;
    xb_vecMxf32 out, update, vps, vnoise; xb_vecMxf32 zero =0;
    xb_vecMx32 one =1; xb_vecMx32 vprob;
    valign psa = PDX_LA_MXF32_PP(ps);
    valign noisea = PDX_LA_MXF32_PP(noise);
    valign proba = PDX_LA_MX32_PP(prob);
    xb_vecMxf32 *pout = (xb_vecMxf32 *)(st->noise);
    valign aout = PDX_LA_MXF32_PP(pout);
    vboolM lthan, flse, or;

    int loopend =(N/8)*8; int byte;
    for (int i = 0; i < N-1; i+=8) 
	{
		byte = MIN(32, 4*(N-i));
		PDX_LAV_MXF32_XP( vps, psa, ps, byte);
		PDX_LAV_MXF32_XP( vnoise, noisea, noise, byte);
		PDX_LAV_MX32_XP( vprob, proba, prob, byte);
		update = vnoise;
		flse =  PDX_LT_MX32(vprob, one); //!st->update_prob[i]
		lthan = PDX_OLT_MXF32(vps, vnoise); //bits =1 when st->ps[i] < st->noise[i]
		or = PDX_OR_BM(lthan, flse);
		PDX_MUL_MXF32_T(update, vbeta1, vnoise, or); //beta_1*st->noise[i]
		PDX_MULA_MXF32_T(update, vbeta, vps, or); //beta*st->ps[i]
		out = PDX_MAX_MXF32( update, zero); //max(0,update)
		PDX_SAV_MXF32_XP( out, aout, pout, byte);
    }
    PDX_SAPOS_MXF32_FP( aout, pout );
}
#endif

#ifdef OVERRIDE_ANR_APOSTERIORI_SNR
static void aposteriori_snr(SpeexPreprocessState2 * st)
{
	//Compute A-posteriori / A-priori SNRs.

    int             N = st->ps_size;
    int             M = st->nbands;
    spx_word32_t   *ps = st->ps;
    int             i;

    xb_vecMxf32 *restrict pnoise = (xb_vecMxf32 *)(st->noise);
    xb_vecMxf32 *restrict pecho = (xb_vecMxf32 *)(st->echo_noise);
    xb_vecMxf32 *restrict preverb = (xb_vecMxf32 *)(st->reverb_estimate);
    xb_vecMxf32 *restrict poldps = (xb_vecMxf32 *)(st->old_ps);
    xb_vecMxf32 *restrict pps = (xb_vecMxf32 *)(st->ps);
    xb_vecMxf32 *restrict pprior = (xb_vecMxf32 *)(st->prior);
    xb_vecMxf32 *restrict ppost = (xb_vecMxf32 *)(st->post);
    valign anoise = PDX_LA_MXF32_PP(pnoise);
    valign aecho = PDX_LA_MXF32_PP(pecho);
    valign areverb = PDX_LA_MXF32_PP(preverb);
    valign aoldps = PDX_LA_MXF32_PP(poldps);
    valign aps = PDX_LA_MXF32_PP(pps);
    valign aprior = PDX_LA_MXF32_PP(pprior);
    valign apost = PDX_LA_MXF32_PP(ppost);
 
    xb_vecMxf32 hundred = 100.f; xb_vecMxf32 one = 1.f; xb_vecMxf32 pointone = .1f; xb_vecMxf32 eightnine = 0.89f; xb_vecMxf32 zero = 0.0f;
    xb_vecMxf32 gamma, tot_noise, vnoise, vecho, vreverb, vps, voldps, val2,val3, val1, out1, out2;

	//total noise = 1 + noise + reverb + echo
	//posteriori SNR = MIN(100, ps/totnoise - 1 )
	// update gamma = .1 + .89*(old/(old+totnoise))^2
	//priori SNR  = gamma*max(0,post) + (1-gamma)*old/totnoise
    for (i = 0; i < N + M; i+=8) 
	{
    	int byte = MIN(32, 4*(N+M-i));
    	PDX_LAV_MXF32_XP( vnoise, anoise, pnoise, byte);
    	PDX_LAV_MXF32_XP( vecho, aecho, pecho, byte);
    	PDX_LAV_MXF32_XP( vreverb, areverb, preverb, byte);
    	PDX_LAV_MXF32_XP( vps, aps, pps, byte);
    	PDX_LAV_MXF32_XP( voldps, aoldps, poldps, byte);
    	tot_noise = PDX_ADD_MXF32(one, vnoise);
    	tot_noise += PDX_ADD_MXF32(vecho, vreverb);
    	out1 = PDX_DIV_MXF32(vps, tot_noise) - one;
    	out1 =  PDX_MIN_MXF32( out1, hundred); //post SNR
    	val1 = PDX_ADD_MXF32(voldps, tot_noise);
    	val2 =  PDX_DIV_MXF32(voldps, val1);
    	val3 = eightnine*val2*val2;
        gamma = PDX_ADD_MXF32(pointone, val3);
        val1 = gamma*PDX_MAX_MXF32(zero, out1);
        val2 = (one - gamma)*PDX_DIV_MXF32(voldps, tot_noise);
        out2 = PDX_ADD_MXF32(val1 , val2);
        out2 =  PDX_MIN_MXF32( out2, hundred); //prior SNR
        PDX_SAV_MXF32_XP( out1,apost, ppost, byte);
        PDX_SAV_MXF32_XP( out2, aprior, pprior, byte);
    }
    PDX_SAPOS_MXF32_FP( aprior, pprior );
    PDX_SAPOS_MXF32_FP( apost, ppost );
}
#endif


#ifdef OVERRIDE_ANR_UPDATE_ZETA
static void preprocess_update_zeta(SpeexPreprocessState2 * st)
{
//Update Smoothed a priori SNR.
    int             N = st->ps_size;
    int             M = st->nbands;
    int             i;

    xb_vecMxf32 * pzeta = (xb_vecMxf32 *)(st->zeta+1);
    valign azeta = PDX_LA_MXF32_PP(pzeta);
    xb_vecMxf32 *pzetaout = (xb_vecMxf32 *)(st->zeta+1);
    valign azetaout = PDX_LA_MXF32_PP(pzetaout);

    xb_vecMxf32 * pprior0 = (xb_vecMxf32 *)(st->prior-1);
    xb_vecMxf32 * pprior1 = (xb_vecMxf32 *)(st->prior);
    xb_vecMxf32 * pprior2 = (xb_vecMxf32 *)(st->prior+1);
    valign aprior0 = PDX_LA_MXF32_PP(pprior0);
    valign aprior1 = PDX_LA_MXF32_PP(pprior1);
    valign aprior2 = PDX_LA_MXF32_PP(pprior2);

    xb_vecMxf32 pointseven = 0.7; xb_vecMxf32 pointthree = 0.3; xb_vecMxf32 pointsevenfive = 0.075f; xb_vecMxf32 pointfifteen = 0.15f;
    xb_vecMxf32 vzeta, vprior0, vprior1, vprior2, vzetaout;

    st->zeta[0] = 0.7f*st->zeta[0] + 0.3f* st->prior[0];
   // st->zeta[i] = 0.7f*st->zeta[i] + 0.15f*st->prior[i] + 0.075f*st->prior[i-1] + 0.075f*st->prior[i+1];
    for (i = 1; i < N - 1; i+=8)
    {
    	int byte = MIN(32, 4*(N-1-i));
		PDX_LAV_MXF32_XP( vzeta, azeta, pzeta, byte);
		PDX_LAV_MXF32_XP( vprior0, aprior0, pprior0, byte);
		PDX_LAV_MXF32_XP( vprior1, aprior1, pprior1, byte);
		PDX_LAV_MXF32_XP( vprior2, aprior2, pprior2, byte);
		vzetaout = pointseven*vzeta + pointfifteen*vprior1 + pointsevenfive*vprior0 + pointsevenfive*vprior2;
		PDX_SAV_MXF32_XP( vzetaout, azetaout, pzetaout, byte);
    }
    PDX_SAPOS_MXF32_FP( azetaout, pzetaout );

    xb_vecMxf32 * pzeta2 = (xb_vecMxf32 *)(st->zeta+N-1);
    valign azeta2 = PDX_LA_MXF32_PP(pzeta2);
    xb_vecMxf32 * pzetaout2 = (xb_vecMxf32 *)(st->zeta+N-1);
    valign azetaout2 = PDX_LA_MXF32_PP(pzetaout2);
    xb_vecMxf32 * pprior11 = (xb_vecMxf32 *)(st->prior+N-1);
    valign aprior11 = PDX_LA_MXF32_PP(pprior11);

	// st->zeta[i] = 0.7f*st->zeta[i] + 0.3f*st->prior[i];
    for (i = N-1; i < N+M; i+=8)
    {
    	int byte = MIN(32, 4*(N+M-i));
		PDX_LAV_MXF32_XP( vzeta, azeta2, pzeta2, byte);
		PDX_LAV_MXF32_XP( vprior1, aprior11, pprior11, byte);
		vzetaout = pointseven*vzeta + pointthree*vprior1;
		PDX_SAV_MXF32_XP( vzetaout, azetaout2, pzetaout2, byte);
    }
    PDX_SAPOS_MXF32_FP( azetaout2, pzetaout2 );
}

#endif


#ifdef OVERRIDE_ANR_UPDATE_GAINS_CRITICAL_BANDS
// Update gains in critical bands (MEL scale).
#include "math.h"
#include "__fenv.h"
#include <math_floating_vec.h>
#include "expf_tbl.h"
#include "inff_tbl.h"
#include "nanf_tbl.h"

//finding exponent of all the values of vector register xin
static inline xb_vecMxf32 register_exp(xb_vecMxf32 xin)
{
	xb_vecMxf32 scl0, scl1, gf, zero_f, zout;
	xb_vecMx32 xin_i, fr, g, t, exp, exp0, exp1, zero_i;
	xb_vecMx80 w0, w1;
	xb_vecMxu32 flags, flags_n, maskx;
	vboolM b_nan, b_snan, b_ovfl, b_edom, b_fe_ovfl, b_fe_inv;
	__fenv_t fenv;
	b_edom = b_fe_ovfl = b_fe_inv = PDX_MOVCI_BM( PDX_MOVBI_N_ZERO );
	__feholdexcept( &fenv );
	zero_i = PDX_ZERO_MX32();
	zero_f = PDX_MOV_MXF32_FROM_4MX8( PDX_MOV_4MX8_FROM_MX32( zero_i ) );
	flags = PDX_ZERO_MXU32();

	xin_i = PDX_TRUNC32_MXF32( xin, 24 );
	w0 = PDX_MULW_MX32( xin_i, invln2_Q30 );
	w1 = PDX_SRAI_MX80( w0, 22 );
	fr = PDX_PACKV_MX80( w1 );
	fr = PDX_SRLI_MX32( fr, 1 );
	w1 = PDX_SRAI_MX80( w0, 54 );
	exp = PDX_PACKV_MX80( w1 );
	{
		xb_vecMx32 y,y1,y2,c1,c2,f2;
		f2=PDX_PACKSIV_MX80(PDX_MULW_MX32( fr,fr ),31);
		y1=PDX_LSR_32_I((const int*)expftbl_Q30,0*sizeof(int32_t));
		y2=PDX_LSR_32_I((const int*)expftbl_Q30,1*sizeof(int32_t));
		c1=PDX_LSR_32_I((const int*)expftbl_Q30,2*sizeof(int32_t)); t = PDX_PACKSIV_MX80(PDX_MULW_MX32(f2, y1), 31);  y1 = PDX_ADD_MX32(c1, t);
		c2=PDX_LSR_32_I((const int*)expftbl_Q30,3*sizeof(int32_t)); t = PDX_PACKSIV_MX80(PDX_MULW_MX32(f2, y2), 31);  y2 = PDX_ADD_MX32(c2, t);
		c1=PDX_LSR_32_I((const int*)expftbl_Q30,4*sizeof(int32_t)); t = PDX_PACKSIV_MX80(PDX_MULW_MX32(f2, y1), 31);  y1 = PDX_ADD_MX32(c1, t);
		c2=PDX_LSR_32_I((const int*)expftbl_Q30,5*sizeof(int32_t)); t = PDX_PACKSIV_MX80(PDX_MULW_MX32(f2, y2), 31);  y2 = PDX_ADD_MX32(c2, t);
		c1=PDX_LSR_32_I((const int*)expftbl_Q30,6*sizeof(int32_t)); t = PDX_PACKSIV_MX80(PDX_MULW_MX32(f2, y1), 31);  y1 = PDX_ADD_MX32(c1, t);
																	t = PDX_PACKSIV_MX80(PDX_MULW_MX32(fr, y2), 31);  y  = PDX_ADD_MX32(y1, t);
		g=y;
	}
	gf = PDX_FLOATF32_MX32( g, 30 );
	exp1 = PDX_SRAI_MX32( exp, 1 );
	exp0 = PDX_SUB_MX32( exp, exp1 );
	exp0 = PDX_ADD_MX32( 127, exp0 );
	exp1 = PDX_ADD_MX32( 127, exp1 );
	exp0 = PDX_SLLI_MX32( exp0, 23 );
	exp1 = PDX_SLLI_MX32( exp1, 23 );
	scl0 = PDX_MOV_MXF32_FROM_4MX8( PDX_MOV_4MX8_FROM_MX32( exp0 ) );
	scl1 = PDX_MOV_MXF32_FROM_4MX8( PDX_MOV_4MX8_FROM_MX32( exp1 ) );
	zout = PDX_MUL_MXF32( gf  , scl0 );
	zout = PDX_MUL_MXF32( zout, scl1 );
	xin_i = PDX_MOV_MX32_FROM_4MX8(PDX_MOV_4MX8_FROM_MXF32( xin ) );
	flags_n = PDX_CLSFY_MXF32(xin);
	flags = PDX_OR_MXU32(flags, flags_n);
	b_nan  = PDX_UN_MXF32( xin, xin );
	b_ovfl = PDX_OLE_MXF32( expfminmax[1].f, xin );
	b_ovfl = PDX_OLT_MXF32_T( xin, plusInff.f, b_ovfl );
	b_fe_ovfl = PDX_OR_BM( b_fe_ovfl, b_ovfl );
	zout = PDX_MOV_MXF32_T( qNaNf.f, zout, b_nan );
	return zout;
}


static inline xb_vecMxf32 hypergeom_gain_vector(xb_vecMxf32 xx, int byte)
{
	//compute hypergeometric function value for all the elements in the vector
	
	static const float32_t table[21] = {
		0.82157f, 1.02017f, 1.20461f, 1.37534f, 1.53363f, 1.68092f, 1.81865f,1.94811f,
		2.07038f, 2.18638f, 2.29688f, 2.40255f, 2.50391f, 2.60144f, 2.69551f, 2.78647f,
		2.87458f, 2.96015f, 3.04333f, 3.12431f, 3.20326f };

	xb_vecMxf32 *ptable = (xb_vecMxf32 *)(table);
	valign atable = PDX_LA_MXF32_PP(ptable);
	xb_vecMxf32 vtable1, vtable2, vtable3,  outtable3;
	int foo = (1 << byte) - 1;
	vboolM vbyte = *((vboolM *) &foo); //boolean vector to get how many bytes in the input vector needs computation (ignore rest)

	PDX_LAV_MXF32_XP( vtable1, atable, ptable, 32); //table's 1st row
	PDX_LAV_MXF32_XP( vtable2, atable, ptable, 32); //table's 2nd row
	PDX_LAV_MXF32_XP( vtable3, atable, ptable, 20); //table's 3rd row

	//if ind<0, hypergeom gain = 1
	//if ind>19, hypergeom gain = (1 + .1296 / x)
	//if 0<ind<19, hypergeom gain =  FRAC_SCALING * ((1 - frac) * table[ind] + frac * table[ind + 1]) / sqrt(x + .0001f);

	xb_vecMxf32 outsmall, outbig, out1, out2, out3, out4;
	xb_vecMxf32 two = 2.0f; xb_vecMxf32 one = 1.f; xb_vecMxf32 pointone = .0001f;
	xb_vecMx32 seven = 7; xb_vecMx32 fifteen = 15; xb_vecMx32 zero = 0; xb_vecMx32 intone = 1; xb_vecMx32 twenty = 20;
	xb_vecMxf32 point1296 = .1296f;

	xb_vecMxf32 x = PDX_MUL_MXF32(one,xx);
	xb_vecMxf32     intg = x*two;
	vboolM ind10, ind11, ind1, ind2, ind20, ind21, ind3;
	intg = PDX_FIFLOOR_MXF32(intg);
	xb_vecMx32 ind = (xb_vecMx32)intg;
	xb_vecMxf32 frac = two*x - intg;
	xb_vecMxf32  out_corner, out_main;
	vboolM vbool_big, vbool_small, vbool_corner, vbool_main1, vbool_main2, vbool_main;
	outsmall = one;
	outbig = one*(PDX_ADD_MXF32(1,PDX_DIV_MXF32(point1296,x)) );
	vbool_big = PDX_GE_MX32(ind,twenty); //ind>=20
	vbool_small = PDX_GT_MX32(zero, ind); //ind<0
	out_corner = PDX_MOV_MXF32_T(outsmall, zero, vbool_small); //if vbool=1, outsmall=1, else 0
	out_corner = PDX_MOV_MXF32_T(outbig, out_corner, vbool_big); //if vbool=1, outcorner=outbig, else outsmall
	ind10 = PDX_GE_MX32(ind, zero);
	ind11 = PDX_GE_MX32(seven, ind);
	ind1 = ind10 & ind11  & vbyte; //0=<ind<=7 &vbyte to avoid 0 index values being recognized
	ind20 = PDX_GE_MX32(fifteen, ind);
	ind21 = PDX_GT_MX32(ind, seven);
	ind2 = ind20 & ind21; // 7<ind<=15
	ind3 = PDX_GT_MX32(ind, fifteen); // ind>15
	xb_vecMxf32 outtable = 0, outtable2 = 0;
	//finding out table values at particular index specified in ind, and saving it in outtable
	PDX_SEL_MXF32_T(outtable, vtable1, vtable1, ind, ind1);
	PDX_SEL_MXF32_T(outtable, vtable2, vtable2, ind, ind2);
	PDX_SEL_MXF32_T(outtable, vtable3, vtable3, ind, ind3);

	ind = ind + intone;
	ind10 = PDX_GE_MX32(ind, zero);
	ind11 = PDX_GE_MX32(seven, ind);
	ind1 = ind10 & ind11 & vbyte; //0=<ind<=7 
	ind20 = PDX_GE_MX32(fifteen, ind);
	ind21 = PDX_GT_MX32(ind, seven);
	ind2 = ind20 & ind21;  // 7<ind<=15
	ind3 = PDX_GT_MX32(ind, fifteen); // ind>15
	//finding out table values at particular index specified in ind, and saving it in outtable
	PDX_SEL_MXF32_T(outtable2, vtable1, vtable1, ind, ind1);
	PDX_SEL_MXF32_T(outtable2, vtable2, vtable2, ind, ind2);
	PDX_SEL_MXF32_T(outtable2, vtable3, vtable3, ind, ind3);
	out1 = (one-frac)*outtable; //(1 - frac) * table[ind]
	out2 = frac*outtable2; //frac * table[ind + 1]
	out3 = one*(out1+out2) ;
	out4 = PDX_SQRT_MXF32(x+pointone);
	out_main = PDX_DIV_MXF32(out3, out4);
	vbool_main1 = PDX_GE_MX32(twenty, ind); 
	vbool_main2 = PDX_GE_MX32(ind, intone); 
	vbool_main = vbool_main2 & vbool_main1; // 1<= ind <= 20
	out_main = PDX_MOV_MXF32_T(out_main, out_corner, vbool_main); //if vbool=1, out=out_main, else out_corner
	return out_main;
}

static void update_gains_critical_bands(SpeexPreprocessState2 * st, spx_word16_t Pframe)
{
	//Update gains in critical bands (MEL scale).
    int i;
    int N = st->ps_size;
    int M = st->nbands;

	xb_vecMxf32 *pps = (xb_vecMxf32 *)(st->ps + N);
	valign aps = PDX_LA_MXF32_PP(pps);
	xb_vecMxf32 *pprior = (xb_vecMxf32 *)(st->prior + N);
	valign aprior = PDX_LA_MXF32_PP(pprior);
	xb_vecMxf32 *ppost = (xb_vecMxf32 *)(st->post + N);
	valign apost = PDX_LA_MXF32_PP(ppost);
	xb_vecMxf32 *poldps = (xb_vecMxf32 *)(st->old_ps + N);
	valign aoldps = PDX_LA_MXF32_PP(poldps);
	xb_vecMxf32 *poldps2 = (xb_vecMxf32 *)(st->old_ps +N);
	valign aoldps2 = PDX_LA_MXF32_PP(poldps2);
	xb_vecMxf32 *pzeta = (xb_vecMxf32 *)(st->zeta + N);
	valign azeta = PDX_LA_MXF32_PP(pzeta);

	xb_vecMxf32 *pgain = (xb_vecMxf32 *)(st->gain + N);
	valign again = PDX_LA_MXF32_PP(pgain);
	xb_vecMxf32 *pgain2 = (xb_vecMxf32 *)(st->gain2 + N);
	valign again2 = PDX_LA_MXF32_PP(pgain2);

	xb_vecMxf32 vtheta, vMM, vps, vprior, vprior1, vP1, vq, vprior_ratio, vpost, vgain, vgain2, voldps,voldps2, vzeta, vout1, vout2, vout;
	xb_vecMxf32 one = 1.f; xb_vecMxf32 pointninteen = 0.199f; xb_vecMxf32 pointfifteen = 0.15f;
	xb_vecMxf32 pointtwo = 0.2f; xb_vecMxf32 pointeight = 0.8f;
	xb_vecMxf32 vframe = Pframe;

	// theta = (prior/(prior+1))*(1+post)
	//gain = MIN(1, (prior/(prior+1))*hypergeomgain)
	// Save old Bark power spectrum = 0.2*oldps[i] + 0.8*ps[i]*gain[i]^2
	//q = 1-frame*(0.199f+0.8f*(1/1+(0.15f/zeta)))
	//gain2 = 1/(1+ ((q/(1-q))*(1+prior)*exp(-theta)) )
    for (i = N; i < N + M; i+=8) 
	{
    	int byte = MIN(32, 4*(N+M - i));
    	PDX_LAV_MXF32_XP( vps, aps, pps, byte);
    	PDX_LAV_MXF32_XP( vprior, aprior, pprior, byte);
    	PDX_LAV_MXF32_XP( vpost, apost, ppost, byte);
    	PDX_LAV_MXF32_XP( voldps2, aoldps2, poldps2, byte);
    	PDX_LAV_MXF32_XP( vzeta, azeta, pzeta, byte);
    	vprior1 = vprior + one;
    	vprior_ratio = PDX_DIV_MXF32(vprior, vprior1);
    	vtheta = vprior_ratio * ( one + vpost);
    	vMM = hypergeom_gain_vector(vtheta, N+M - i); //calculate hypergeom gain for the entire vector values

    	xb_vecMxf32 val = vprior_ratio * vMM;
    	vgain = PDX_MIN_MXF32(one, val);
    	PDX_SAV_MXF32_XP( vgain, again, pgain, byte);
    	voldps = pointtwo * voldps2 + pointeight * vgain * vgain *vps;
    	PDX_SAV_MXF32_XP( voldps, aoldps, poldps, byte);

        vP1 = pointninteen + pointeight * PDX_DIV_MXF32(one, (one +  PDX_DIV_MXF32(pointfifteen, (one * vzeta)) ) );
        vq = one - vframe * vP1;
        vout1 = PDX_DIV_MXF32(vq, one - vq);
        xb_vecMxf32 vminustheta = -vtheta;
        xb_vecMxf32 vexp =  register_exp(vminustheta);
        vout2 = vout1 * (one + vprior)*vexp ;
        vout = one + vout2;
        vgain2 = PDX_DIV_MXF32(one, vout) ;
        PDX_SAV_MXF32_XP( vgain2, again2, pgain2, byte);
	}

    PDX_SAPOS_MXF32_FP(again2, pgain2);
    PDX_SAPOS_MXF32_FP(aoldps, poldps);
    PDX_SAPOS_MXF32_FP(again, pgain);
}
#endif

#ifdef OVERRIDE_ANR_UPDATE_GAINS_LINEAR
static void update_gains_linear(SpeexPreprocessState2 * st)
{
	//Update gains in linear spectral bands.
    int             i;
    int             N = st->ps_size;
    int             M = st->nbands;
    spx_word32_t   *ps = st->ps;

	xb_vecMxf32 *pprior = (xb_vecMxf32 *)(st->prior);
	valign aprior = PDX_LA_MXF32_PP(pprior);
	xb_vecMxf32 *ppost = (xb_vecMxf32 *)(st->post);
	valign apost = PDX_LA_MXF32_PP(ppost);

	xb_vecMxf32 *pgain = (xb_vecMxf32 *)(st->gain);
	valign again = PDX_LA_MXF32_PP(pgain);
	xb_vecMxf32 *pgainout = (xb_vecMxf32 *)(st->gain);
	valign againout = PDX_LA_MXF32_PP(pgainout);

	xb_vecMxf32 *pgain2 = (xb_vecMxf32 *)(st->gain2);
	valign again2 = PDX_LA_MXF32_PP(pgain2);
	xb_vecMxf32 *pgain2out = (xb_vecMxf32 *)(st->gain2);
	valign again2out = PDX_LA_MXF32_PP(pgain2out);

	xb_vecMxf32 *pps = (xb_vecMxf32 *)(st->ps);
	valign aps = PDX_LA_MXF32_PP(pps);
	xb_vecMxf32 *poldps = (xb_vecMxf32 *)(st->old_ps);
	valign aoldps = PDX_LA_MXF32_PP(poldps);
	xb_vecMxf32 *poldpsout = (xb_vecMxf32 *)(st->old_ps);
	valign aoldpsout = PDX_LA_MXF32_PP(poldpsout);
	xb_vecMxf32 *pgainfloor = (xb_vecMxf32 *)(st->gain_floor);
	valign againfloor = PDX_LA_MXF32_PP(pgainfloor);


	xb_vecMxf32 vprior, vpost, vprior1, vprior_ratio, vtheta, vMM, vtmp, vg, vp, vps;
    xb_vecMxf32 one = 1.0f; xb_vecMxf32 pointthree = 0.333f; xb_vecMxf32 pointtwo = 0.2f; xb_vecMxf32 pointeight = 0.8f;
    xb_vecMxf32 three = 3.0f;

	// theta = (prior/(prior+1))*(1+post)
	//gain = MIN(1, (prior/(prior+1))*hypergeomgain)
	// oldpsout = 0.2*oldps[i] + 0.8*ps[i]*gain[i]^2
	//gain2 = [p*sqrt(gain)+(1-p)*sqrt(gain _floor) ]^2
    xb_vecMxf32 vgain, vgainout, voldps, voldpsout, vgainfloor, vgain2, vgain2out;
    for (int i = 0; i < N; i+=8) 
	{
		int byte = MIN(32, 4*(N - i));
		PDX_LAV_MXF32_XP( vprior, aprior, pprior, byte);
		PDX_LAV_MXF32_XP( vpost, apost, ppost, byte);
		vprior1 = vprior + one;
		vprior_ratio = PDX_DIV_MXF32(vprior, vprior1);
		vtheta = vprior_ratio * ( one + vpost);
		vMM = hypergeom_gain_vector(vtheta, byte/4);

		PDX_LAV_MXF32_XP( vgain, again, pgain, byte);
		PDX_LAV_MXF32_XP(vgain2, again2, pgain2, byte);
		PDX_LAV_MXF32_XP( voldps, aoldps, poldps, byte);
		PDX_LAV_MXF32_XP( vps, aps, pps, byte);
		PDX_LAV_MXF32_XP( vgainfloor, againfloor, pgainfloor, byte);
		xb_vecMxf32 vvg = PDX_MIN_MXF32(one, vprior_ratio*vMM);
		vboolM vboolval = PDX_OLT_MXF32(vgain, pointthree*vvg);
		vgainout = PDX_MOV_MXF32_T(three*vgain, vvg, vboolval);
		voldpsout = pointtwo*voldps + pointeight*vps*vgainout*vgainout;

		vboolM vboolval2 = PDX_OLT_MXF32(vgainfloor, vgainout);
		vgainout = PDX_MOV_MXF32_T(vgainout, vgainfloor, vboolval2);

		xb_vecMxf32 vtmp = vgain2*PDX_SQRT_MXF32(vgainout) + (one-vgain2)*PDX_SQRT_MXF32(vgainfloor);
		vgain2out = vtmp*vtmp;

		PDX_SAV_MXF32_XP(vgainout, againout, pgainout, byte);
		PDX_SAV_MXF32_XP(voldpsout, aoldpsout, poldpsout, byte);
		PDX_SAV_MXF32_XP(vgain2out, again2out, pgain2out, byte);
    }
    PDX_SAPOS_MXF32_FP(againout, pgainout);
    PDX_SAPOS_MXF32_FP(aoldpsout, poldpsout);
    PDX_SAPOS_MXF32_FP(again2out, pgain2out);
}
#endif

#ifdef OVERRIDE_ANR_APPLY_SPEC_GAIN
static void apply_spectral_gain(SpeexPreprocessState2 * st)
{
	//Apply computed spectral gain.
    int             i;
    int             N = st->ps_size;

    xb_vecMxf32 *pft = (xb_vecMxf32 *)(st->ft+1);
    valign aft = PDX_LA_MXF32_PP(pft);
    xb_vecMxf32 *pftout = (xb_vecMxf32 *)(st->ft+1);
    valign aftout = PDX_LA_MXF32_PP(pftout);
    xb_vecMxf32 *pgain2 = (xb_vecMxf32 *)(st->gain2+1);
    valign again2 = PDX_LA_MXF32_PP(pgain2);
    xb_vecMxf32 vft1, vft2, vgain2;
    for (i = 1; i < N; i+=8) 
	{
    	int byte = 4*(N-i);
		PDX_LAV_MXF32_XP( vft1, aft, pft, 2*byte); //load first 8 values
		PDX_LAV_MXF32_XP( vft2, aft, pft, 2*byte-32); //load next 8-byte/4 values
		PDX_LAV_MXF32_XP( vgain2, again2, pgain2, byte);
		PDX_DSELI_MXF32(vft2, vft1, vft2, vft1, PDX_DSELI_32B_DEINTERLEAVE_1); //separate 2*i and 2*i-1
		vft1 = vft1*vgain2; //st->gain2[i], st->ft[2 * i]
		vft2 = vft2*vgain2; //st->gain2[i], st->ft[2 * i - 1]
		PDX_DSELI_MXF32(vft2, vft1, vft2, vft1, PDX_DSELI_32B_INTERLEAVE_1);
		PDX_SAV_MXF32_XP( vft1, aftout, pftout, 2*byte);
		PDX_SAV_MXF32_XP( vft2, aftout, pftout, 2*byte-32);
    }
    PDX_SAPOS_MXF32_FP(aftout, pftout);
    st->ft[0] = st->gain2[0] * st->ft[0];
    st->ft[2 * N - 1] = st->gain2[N - 1] * st->ft[2 * N - 1];
}

#endif


#ifdef OVERRIDE_ANR_UPDATE_NOISE_PROB
static void update_noise_prob(SpeexPreprocessState2 * st)
{
	// Update noise probabilities and smoothed power spectrum. 
	int             i;
	int             min_range;
	int             N = st->ps_size;


	if (st->nb_adapt == 1) 
	{
		for (i = 0; i < N; i++)
			st->Smin[i] = st->Stmp[i] = 0;
	}
	if (st->nb_adapt < 100)
		min_range = 15;
	else if (st->nb_adapt < 1000)
		min_range = 50;
	else if (st->nb_adapt < 10000)
		min_range = 150;
	else
		min_range = 300;

	xb_vecMxf32 vS, vSout, vSmin, vSminout, vStmp, vStmpout, vps0, vps1, vps2;
	st->S[0] = 0.8f*st->S[0] + 0.2f*st->ps[0];
	st->S[N - 1] = 0.8f*st->S[N - 1] + 0.2f* st->ps[N - 1];

	xb_vecMxf32 *pS = (xb_vecMxf32 *)(st->S + 1);
	valign aS = PDX_LA_MXF32_PP(pS);
	xb_vecMxf32 *pSout = (xb_vecMxf32 *)(st->S + 1);
	valign aSout = PDX_LA_MXF32_PP(pS);
	xb_vecMxf32 *pps0 = (xb_vecMxf32 *)(st->ps);
	valign aps0 = PDX_LA_MXF32_PP(pps0);
	xb_vecMxf32 *pps1 = (xb_vecMxf32 *)(st->ps + 1);
	valign aps1 = PDX_LA_MXF32_PP(pps1);
	xb_vecMxf32 *pps2 = (xb_vecMxf32 *)(st->ps + 2);
	valign aps2 = PDX_LA_MXF32_PP(pps2);
	xb_vecMxf32 pointeight = 0.8f; xb_vecMxf32 pointone = 0.1f;
	xb_vecMxf32 pointfive = 0.05f; xb_vecMxf32 pointfour = 0.4f;
	xb_vecMx32 zero = 0; xb_vecMx32 one = 1;
	// Sout = 0.8*S[i] + 0.5*ps[i] + 0.1*ps[i+1] + 0.5*ps[i+2]
	for (i = 1; i < N - 1; i+=8)
	{
		int byte = MIN(32, 4*(N-1 - i));
		PDX_LAV_MXF32_XP( vps0, aps0, pps0, byte);
		PDX_LAV_MXF32_XP( vps1, aps1, pps1, byte);
		PDX_LAV_MXF32_XP( vps2, aps2, pps2, byte);
		PDX_LAV_MXF32_XP( vS, aS, pS, byte);
		vSout = pointeight*vS + pointfive*vps0 + pointone*vps1 + pointfive*vps2;
		PDX_SAV_MXF32_XP( vSout, aSout, pSout, byte);
	}
	PDX_SAPOS_MXF32_FP(aSout, pSout);
	
	xb_vecMx32 *pupdateprob = (xb_vecMx32 *)(st->update_prob);
	valign aupdateprob = PDX_LA_MX32_PP(pupdateprob);
	xb_vecMx32 intupdateprob;  vboolM updateprob;

	xb_vecMxf32 *pSmin = (xb_vecMxf32 *)(st->Smin);
	valign aSmin = PDX_LA_MXF32_PP(pSmin);
	xb_vecMxf32 *pStmp = (xb_vecMxf32 *)(st->Stmp);
	valign aStmp = PDX_LA_MXF32_PP(pStmp);
	xb_vecMxf32 *pSminout = (xb_vecMxf32 *)(st->Smin);
	valign aSminout = PDX_LA_MXF32_PP(pSminout);
	xb_vecMxf32 *pStmpout = (xb_vecMxf32 *)(st->Stmp);
	valign aStmpout = PDX_LA_MXF32_PP(pStmpout);

	xb_vecMxf32 *pS1 = (xb_vecMxf32 *)(st->S);
	valign aS1 = PDX_LA_MXF32_PP(pS1);

	if (st->min_count > min_range) 
	{
		st->min_count = 0;
		for (i = 0; i < N; i+=8) 
		{
			int byte = MIN(32, 4*(N- i));
			PDX_LAV_MXF32_XP(vStmp, aStmp, pStmp, byte);
			PDX_LAV_MXF32_XP(vS, aS1, pS1, byte);
			vSminout = PDX_MIN_MXF32(vStmp, vS);
			vStmpout = vS;
			xb_vecMxf32 vcomp = pointfour*vS;
			updateprob = PDX_OLT_MXF32(vSminout, vcomp);
			intupdateprob = PDX_MOV_MX32_T(one, zero, updateprob);
			PDX_SAV_MX32_XP( intupdateprob, aupdateprob, pupdateprob, byte);
			PDX_SAV_MXF32_XP( vSminout, aSminout, pSminout, byte);
			PDX_SAV_MXF32_XP( vStmpout, aStmpout, pStmpout, byte);
		}
		PDX_SAPOS_MXF32_FP(aSminout, pSminout);
		PDX_SAPOS_MXF32_FP(aStmpout, pStmpout);
		PDX_SAPOS_MX32_FP(aupdateprob, pupdateprob);
	}
	else 
	{
		for (i = 0; i < N; i+=8) 
		{
			int byte = MIN(32, 4*(N- i));
			PDX_LAV_MXF32_XP(vStmp, aStmp, pStmp, byte);
			PDX_LAV_MXF32_XP(vSmin, aSmin, pSmin, byte);
			PDX_LAV_MXF32_XP(vS, aS1, pS1, byte);
			vSminout = PDX_MIN_MXF32(vSmin, vS);
			vStmpout =  PDX_MIN_MXF32(vStmp, vS);
			xb_vecMxf32 vcomp = pointfour*vS;
			updateprob = PDX_OLT_MXF32(vSminout, vcomp);
			intupdateprob = PDX_MOV_MX32_T(one, zero, updateprob);
			PDX_SAV_MX32_XP( intupdateprob, aupdateprob, pupdateprob, byte);
			PDX_SAV_MXF32_XP( vSminout, aSminout, pSminout, byte);
			PDX_SAV_MXF32_XP( vStmpout, aStmpout, pStmpout, byte);
		}
		PDX_SAPOS_MXF32_FP(aSminout, pSminout);
		PDX_SAPOS_MXF32_FP(aStmpout, pStmpout);
		PDX_SAPOS_MX32_FP(aupdateprob, pupdateprob);
	}
}
#endif

