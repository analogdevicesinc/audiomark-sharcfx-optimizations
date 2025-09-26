/* ----------------------------------------------------------------------
 * Title:        adi_sharcfx_cmplx_dot_prod_f32.c
 * Description:  Floating-point complex dot product
 * $Date:        22 April 2025
 * $Revision:    V.1.0
 *
 * Target Processor: SHARC-FX Processor
 * -------------------------------------------------------------------- */
/*
 * Copyright(c) 2025 Analog Devices, Inc. All Rights Reserved.
 * Copyright (C) 2010-2021 ARM Limited or its affiliates. All rights reserved.
 *
 * SPDX-License-Identifier: Apache-2.0
 *
 * Licensed under the Apache License, Version 2.0 (the License); you may
 * not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an AS IS BASIS, WITHOUT
 * WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */



#include <libs/DSP/Include/adi_sharcfx_math_types.h>
#include <vector.h>
#include <xtensa/sim.h>
#include <xtensa/tie/xt_pdxn.h>

/**
  @brief         Floating-point complex dot product.
  @param[in]     pSrcA       points to the first input vector
  @param[in]     pSrcB       points to the second input vector
  @param[in]     numSamples  number of samples in each vector
  @param[out]    realResult  real part of the result returned here
  @param[out]    imagResult  imaginary part of the result returned here
  @return        none
 */

#define MIN(X, Y) (((X) < (Y)) ? (X) : (Y))

void adi_sharcfx_cmplx_dot_prod_f32(
  const float32_t *restrict pSrcA,
  const float32_t *restrict pSrcB,
        uint32_t numSamples,
        float32_t * realResult,
        float32_t * imagResult)
{
	int N = numSamples;
	xb_vecMxf32 Are, Aim, Bre, Bim;
	xb_vecMxf32 Zre = 0, Zim = 0;
	float  real_sum, imag_sum;
	xb_vecMxf32 *restrict Ap = (xb_vecMxf32*)(pSrcA);
	xb_vecMxf32 *restrict Bp = (xb_vecMxf32*)(pSrcB);
	valign Aa = PDX_LA_MXF32_PP(Ap);
	valign Ba = PDX_LA_MXF32_PP(Bp);

	for ( int n = 0; n < 2*N; n+=16 )
	{
		int byte1 = MIN(32, 4*(2*N-n));
		int byte2 = MIN(32, 4*(2*N-n-8));
		PDX_LAV_MXF32_XP( Are, Aa, Ap, byte1 );
		PDX_LAV_MXF32_XP( Aim, Aa, Ap, byte2 );
		PDX_LAV_MXF32_XP( Bre, Ba, Bp, byte1 );
		PDX_LAV_MXF32_XP( Bim, Ba, Bp, byte2 );
		PDX_DSELI_MXF32(Aim, Are, Aim, Are, PDX_DSELI_32B_DEINTERLEAVE_1); //separating real and imaginary parts
		PDX_DSELI_MXF32(Bim, Bre, Bim, Bre, PDX_DSELI_32B_DEINTERLEAVE_1); //separating real and imaginary parts
		PDX_MULA_MXF32(Zre, Are, Bre); PDX_MULA_MXF32(Zre, -Aim, Bim); //re = A.re*B.re - A.im*B.im
		PDX_MULA_MXF32(Zim, Are, Bim); PDX_MULA_MXF32(Zim, Aim, Bre); //im = A.re*B.im + A.im*B.re
	}
	*realResult = PDX_RADD_MXF32(Zre); //Reduce add all the vector elements into one scaler
	*imagResult = PDX_RADD_MXF32(Zim);
}
