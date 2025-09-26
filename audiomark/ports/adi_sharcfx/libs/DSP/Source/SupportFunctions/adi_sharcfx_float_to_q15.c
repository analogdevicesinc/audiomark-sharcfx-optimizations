/* ----------------------------------------------------------------------
 * Title:        adi_sharcfx_float_to_q15.c
 * Description:  Converts the elements of the floating-point vector to Q15 vector
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
  @brief         Converts the elements of the floating-point vector to Q15 vector.
  @param[in]     pSrc       points to the floating-point input vector
  @param[out]    pDst       points to the Q15 output vector
  @param[in]     blockSize  number of samples in each vector
  @return        none
  */

void adi_sharcfx_float_to_q15(
  const float32_t * pSrc,
        q15_t * pDst,
        uint32_t blockSize)
{

	xb_vecMx16 * __restrict pvOut = (xb_vecMx16 *)(pDst);
	valign pvOuta = PDX_LA_MX16_PP(pvOut);
	xb_vecMxf32 *pA = (xb_vecMxf32 *)(pSrc);
	valign pAa = PDX_LA_MXF32_PP(pA);
	xb_vecMxf32 A, out_scale; xb_vecMx32 out;
	xb_vecMxf32 factor =32768;

    for (int32_t i = 0; i < blockSize; i+=8)
    {
    	PDX_LA_MXF32_XP(A, pAa, pA, 32);
    	out_scale = PDX_MUL_MXF32(A, factor);
    	out =  PDX_TRUNC32_MXF32(out_scale,0); //immediate value divides input by 2^0
    	PDX_SAV32_MX16_XP(out, pvOuta, pvOut, 32);
	}
    PDX_SAPOS_MX16_FP(pvOuta,pvOut);

}
