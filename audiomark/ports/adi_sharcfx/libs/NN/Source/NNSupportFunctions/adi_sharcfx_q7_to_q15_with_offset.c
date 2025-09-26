/*
 * Copyright(c) 2025 Analog Devices, Inc. All Rights Reserved.

 * SPDX-FileCopyrightText: Copyright 2010-2022 Arm Limited and/or its affiliates <open-source-office@arm.com>
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
/* ----------------------------------------------------------------------
 * Title:        adi_sharcfx_q7_to_q15_with_offset.c
 * Description:  Converts the elements of the Q7 vector to Q15 vector with an added offset
 * $Date:        22 April 2025
 * $Revision:    V.1.0
 * 
 * Target Processor: SHARC-FX Processor
 * -------------------------------------------------------------------- */


#include <libs/NN/Include/adi_sharcfx_nnsupportfunctions.h>


void adi_sharcfx_q7_to_q15_with_offset(const int8_t *src, int16_t *dst, uint32_t block_size, int16_t offset)
{
    xb_vec2Mx16 *restrict pvOut  = (xb_vec2Mx16 *)dst;
    valign pvOuta = PDX_LA_2MX16_PP(pvOut);

    xb_vec2Mx8 *restrict pSrc  = (xb_vec2Mx8 *)(src); //B 1st row
    valign pSrca = PDX_LA_2MX8_PP (pSrc);

    xb_vec2Mx16 off = offset; xb_vec2Mx16 out1;  xb_vec2Mx16 oSrc1;
    int loopend = (block_size/16)*16;
    for(int block_count = 0; block_count < loopend; block_count +=16)
	{
    	PDX_LV16_2MX8_IP(oSrc1, pSrc, 16);
    	out1 = PDX_ADDS_2MX16(off, oSrc1);
    	PDX_SA_2MX16_IP(out1, pvOuta, pvOut);
	}

    int inc = block_size -loopend;
	PDX_LAV16_2MX8_XP(oSrc1, pSrca, pSrc, inc);
	out1 = PDX_ADDS_2MX16(off, oSrc1);
	PDX_SAV_2MX16_XP(out1, pvOuta, pvOut, 2*inc);

    PDX_SAPOS_2MX16_FP(pvOuta,pvOut);

}