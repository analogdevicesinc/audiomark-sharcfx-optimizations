/* ----------------------------------------------------------------------
 * Title:        adi_sharcfx_q15_to_float.c
 * Description:  Converts the elements of the Q15 vector to floating-point vector
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
#include <xtensa/sim.h>
#include <xtensa/tie/xt_pdxn.h>


/**
  @brief         Converts the elements of the Q15 vector to floating-point vector.
  @param[in]     pSrc       points to the Q15 input vector
  @param[out]    pDst       points to the floating-point output vector
  @param[in]     blockSize  number of samples in each vector
  @return        none
 */


void adi_sharcfx_q15_to_float(
  const q15_t * pSrc,
        float32_t * pDst,
        uint32_t blockSize)
{

		xb_vecMxf32 * __restrict pvOut = (xb_vecMxf32 *)(pDst);
		valign pvOuta = PDX_LA_MXF32_PP(pvOut);
		xb_vecMx16 *pA = (xb_vecMx16 *)(pSrc);
		valign pAa = PDX_LA_MX16_PP(pA);
		xb_vecMx32 A; xb_vecMxf32 out, out_scale;
		xb_vecMxf32 factor =32768;

	    for (int32_t i = 0; i < blockSize; i+=8)
	    {
	    	PDX_LA32_MX16_XP(A, pAa, pA, 16);
	    	out =  PDX_FLOATF32_MX32(A,0); //immediate value divides input by 2^0
	    	out_scale = PDX_DIV_MXF32(out, factor);
	    	PDX_SAV_MXF32_XP(out_scale, pvOuta, pvOut, 32);
		}
	    PDX_SAPOS_MXF32_FP(pvOuta,pvOut);

}
