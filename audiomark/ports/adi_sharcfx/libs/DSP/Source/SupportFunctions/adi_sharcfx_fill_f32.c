/* ----------------------------------------------------------------------
 * Title:        adi_sharcfx_fill_f32.c
 * Description:  Fills a constant value into a floating-point vector
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
  @brief         Fills a constant value into a floating-point vector.
  @param[in]     value      input value to be filled
  @param[out]    pDst       points to output vector
  @param[in]     blockSize  number of samples in each vector
  @return        none
 */

void adi_sharcfx_fill_f32(
  float32_t value,
  float32_t * pDst,
  uint32_t blockSize)
{

	xb_vecMxf32 In = value;
	xb_vecMxf32 *pvOut = (xb_vecMxf32 *)(pDst);
	valign pvOuta = PDX_LA_MXF32_PP(pvOut);

    for (int i = 0; i < blockSize; i += 8)
    {
    	PDX_SAV_MXF32_XP(In, pvOuta, pvOut, 32);
   }
   PDX_SAPOS_MXF32_FP(pvOuta,pvOut);

}

