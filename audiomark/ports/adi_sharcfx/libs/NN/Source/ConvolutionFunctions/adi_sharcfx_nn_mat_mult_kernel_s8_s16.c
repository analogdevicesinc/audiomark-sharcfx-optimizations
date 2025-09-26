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
 * Title:        adi_sharcfx_nn_mat_mult_kernel_s8_s16.c
 * Description:  Matrix-multiplication function for convolution
 * $Date:        22th April 2025
 * $Revision:    V.1.0
 *
 * Target Processor:  SHARC-FX Processor
 *
 * -------------------------------------------------------------------- */

#include <sys/platform.h>
#include <stdio.h>

#include <stdio.h>
#include <math.h>
#include <float.h>
#include <complex.h>

#include <matrix.h>
#include <filter.h>
#include <libs/NN/Include/adi_sharcfx_nnfunctions.h>
#include <libs/NN/Include/adi_sharcfx_nnsupportfunctions.h>
#include <vector.h>

#define NVEC 32

#include "math_fixedpoint_vec.h"
/* Cross-platform data type definitions. */
#include "libdsp_types.h"


#pragma once
#include <inttypes.h>

#ifdef __XTENSA__
#include <xtensa/sim.h>
#include <xtensa/tie/xt_pdxn.h>
#endif

// Unaligned 16-way by 8b load
#define LOADU_2MX8(dst,ptr,inc) \
	ptr ## a = PDX_LA_2MX8_PP (ptr); \
    PDX_LAV16_2MX8_XP(dst, ptr ## a, ptr, inc)

#define REQUANT(i) \
        /* Multiply */ \
        vScaled = vDot ## i * vScale ## i; \
        /* Handle left shifts */ \
        vScaled1 = PDX_SLS_MX80 (vScaled, vShift ## i); \
        /* Get left shift outputs and cut off low bits for right shifts */ \
        vOut_ge_zero = PDX_PACKQSRV_MX80 (vScaled1, 0); \
        /* Put back for right shifts */ \
		vScaled2 = vOut_ge_zero; \
		/* Move up to Q32 and do the right shifts */ \
		vScaled3 = PDX_SLS_MX80 (vScaled2, v32_m_shift ## i); \
        /* Saturate, round, remove low 32 bits */ \
        vOut_lt_zero = PDX_PACKQSRV_MX80 (vScaled3, 1); \
        /* Select between outputs for shift >= 0 and < 0 */ \
        vOut = PDX_MOV_MX32_T (vOut_ge_zero, vOut_lt_zero, en_sh_ge_zero ## i); \
        /* Add yet another bias */ \
        vOut += vOutZeroPoint; \
        /* Saturate to 8 bits */ \
        vOut = PDX_MIN_MX32 (PDX_MAX_MX32 (vOut, activation_min), activation_max)


/*
 * Matrix-multiplication function for convolution with per-channel requantization.
 *
 * Refer header file for details.
 *
 */
int8_t *adi_sharcfx_nn_mat_mult_kernel_s8_s16(const int8_t *input_a,
                                      const int16_t *input_b,
                                      const uint16_t output_ch,
                                      const int32_t *out_shift,
                                      const int32_t *out_mult,
                                      const int32_t out_offset,
                                      const int16_t activation_min,
                                      const int16_t activation_max,
                                      const uint16_t num_col_a,
                                      const int32_t *const output_bias,
                                      int8_t *out_0)
{
	xb_vec2Mx16 ina0, ina1, ina2, ina3, inb0, inb1;
	xb_vecMx80 vScaled, vScaled1, vScaled2, vScaled3;
	xb_vecMx32 vOut_ge_zero, vOut_lt_zero;
	xb_vecMx32 v32 = 32, vOne = 1;
	xb_vecMx32 vDot0,vShift0,vBias0,vScale0;
	xb_vecMx32 vDot1,vShift1,vBias1,vScale1;
    xb_vecMx32 vOutZeroPoint = out_offset;
    xb_vecMx32 a00, a01, a10, a11, a20, a21, a30, a31;
    xb_vecMx32 vOut;

    //first output row
    xb_vecMx8*restrict pvOut0  = (xb_vecMx8 *)out_0;
    valign pvOuta0 = PDX_LA_MX8_PP(pvOut0);
    //second output row
    xb_vecMx8*restrict pvOut1  = (xb_vecMx8 *)(out_0 + output_ch );
    valign pvOuta1 = PDX_LA_MX8_PP(pvOut1);
    // pointers to Bias, Shift, and Scale for each A row/output column
	xb_vecMx32 *restrict pBias  = (xb_vecMx32 *)(output_bias);
	valign pBiasa = PDX_LA_MX32_PP (pBias);
	xb_vecMx32 *restrict pScale  = (xb_vecMx32 *) (out_mult);
	valign pScalea = PDX_LA_MX32_PP (pScale);
	xb_vecMx32 *restrict pShift  = (xb_vecMx32 *) (out_shift);
	valign pShifta = PDX_LA_MX32_PP (pShift);

    const int8_t *a_row_ptr = input_a; // A first row
    const int16_t *b_row_ptr = input_b; // B first row

    for (int32_t a_rows_idx = 0; a_rows_idx < output_ch; a_rows_idx +=16) // A rows from 0-64 channels
    {
    	xb_int32 add_out[32] = {0}; //getting 32 output elements at a time 16 for 1st output row, 16 for 2nd output row
    	int row_count = 0;
    	//loop till all the a_rows_idx selected are complete
    	while(row_count< MIN(16, output_ch -  a_rows_idx))
        {
			xb_vec2Mx16 *inbp0  = (xb_vec2Mx16 *)(b_row_ptr); //B 1st row
			xb_vec2Mx16 *inbp1  = (xb_vec2Mx16 *)(b_row_ptr + num_col_a);  //B 2nd row
			valign bp0a = PDX_LA_2MX16_PP (inbp0);
			valign bp1a = PDX_LA_2MX16_PP (inbp1);

			xb_vec2Mx8 *ap0 = (xb_vec2Mx8 *)(a_row_ptr ); //A 1st row
			xb_vec2Mx8 *ap1 = (xb_vec2Mx8 *)(a_row_ptr + 1*num_col_a); //A 2nd row
			valign ap0a = PDX_LA_2MX8_PP (ap0);
			valign ap1a = PDX_LA_2MX8_PP (ap1);

    		xb_vec2Mx40 a0 = 0; xb_vec2Mx40 a1 = 0; xb_vec2Mx40 a2 = 0; xb_vec2Mx40 a3 = 0; //initialize accumulator results
    		// iterating over all columns of a = 40, 16 at a time
    		for (int32_t a_cols_idx = 0; a_cols_idx < num_col_a; a_cols_idx+= 16)
			{
    			int inc = MIN(16, num_col_a - a_cols_idx );
				PDX_LA_2MX16_IP(inb0,bp0a,inbp0); // B[0][0-15]
				PDX_LA_2MX16_IP(inb1,bp1a,inbp1);//  B[1][0-15]
				PDX_LAV16_2MX8_XP(ina0, ap0a, ap0, inc); // A[a_row_idx][0-15]
				PDX_LAV16_2MX8_XP(ina1, ap1a, ap1, inc);// A[a_row_idx next row][0-15]
				// accumulate the multiplication for
				//A[a_row_idx0][all], A[a_row_idx next row][all] with B[0][all], B[1][all]
				PDX_MULAW_2MX16(a0,ina0,inb0);
				PDX_MULAW_2MX16(a1,ina1,inb0);
				PDX_MULAW_2MX16(a2,ina0,inb1);
				PDX_MULAW_2MX16(a3,ina1,inb1);
            }
			// convert 40 bit to two 8 lane 32 bit vectors
			PDX_CVT32D_2MX40(a00, a01, a0);
			PDX_CVT32D_2MX40(a10, a11, a1);
			PDX_CVT32D_2MX40(a20, a21, a2);
			PDX_CVT32D_2MX40(a30, a31, a3);
			//adding all the multiplied elements together to get the vector multiplication result
			//first output row
			add_out[row_count + 0]  += PDX_RADDS_MX32(a00) + PDX_RADDS_MX32(a01); //8-way 32-bit Signed Add operation, reduction to scalar with saturation to int32. PDX_RADD_2MX16
			add_out[row_count + 1]  += PDX_RADDS_MX32(a10) + PDX_RADDS_MX32(a11);
			//second output row
			add_out[16 + row_count + 0]  += PDX_RADDS_MX32(a20) + PDX_RADDS_MX32(a21);
			add_out[16 + row_count + 1]  += PDX_RADDS_MX32(a30) + PDX_RADDS_MX32(a31);

			a_row_ptr +=2*num_col_a; //update to next two rows
			row_count += 2; //increase row count till 16 rows
        }
		//loading two 8 way 32 bit vector for bias
		PDX_LA_MX32_IP (vBias0, pBiasa, pBias);
		PDX_LA_MX32_IP (vBias1, pBiasa, pBias);
		//loading two 8 way 32 bit vector for scale
		PDX_LA_MX32_IP (vScale0, pScalea, pScale);
		PDX_LA_MX32_IP (vScale1, pScalea, pScale);
		//loading two 8 way 32 bit vector for shift
		PDX_LA_MX32_IP (vShift0, pShifta, pShift);
		PDX_LA_MX32_IP (vShift1, pShifta, pShift);
		//loading two 8 boolean vectors to separate positive and negative shifts
		vboolM en_sh_ge_zero0 = PDX_GE_MX32 (vShift0, 0);
		vboolM en_sh_ge_zero1 = PDX_GE_MX32 (vShift1, 0);
		//add 32bits to shift value
		xb_vecMx32 v32_m_shift0 = v32 + vShift0;
		xb_vecMx32 v32_m_shift1 = v32 + vShift1;
		// Add one to convert Q31 to Q32, use 1 if shift < 0
		vShift0 = PDX_MOV_MX32_T (vShift0 + vOne, 1, en_sh_ge_zero0);
		vShift1 = PDX_MOV_MX32_T (vShift1 + vOne, 1, en_sh_ge_zero1);
		//Calculate number of columns to process as nbytestowrite0
		int nBytesToWrite0 = MIN(16, output_ch - a_rows_idx);
		int nBytesToWrite1 = nBytesToWrite0 - PDX_M;
		//vector pointer to the addition result obtained
		xb_int32 *add_ptr = add_out;
		xb_vecMx32 *vAdd0 = (xb_vecMx32 *)(add_ptr);
		valign wAdd0 =  PDX_LA_MX32_PP(vAdd0);
		//loading first 8 way 32 bit vector output first row
		PDX_LA_MX32_IP (vDot0, wAdd0, vAdd0);
		vDot0 += vBias0;
		REQUANT(0);
		PDX_SAV32_MX8_XP (vOut, pvOuta0, pvOut0, nBytesToWrite0);
		//loading first 8 way 32 bit vector output first row
		PDX_LA_MX32_IP (vDot1, wAdd0, vAdd0);
		vDot1 += vBias1;
		REQUANT(1);
		PDX_SAV32_MX8_XP (vOut, pvOuta0, pvOut0, nBytesToWrite1);

		//loading first 8 way 32 bit vector output second row
		PDX_LA_MX32_IP (vDot0, wAdd0, vAdd0);
		vDot0 += vBias0;
		REQUANT(0);
		PDX_SAV32_MX8_XP (vOut, pvOuta1, pvOut1, nBytesToWrite0);
		//loading first 8 way 32 bit vector output second row
		PDX_LA_MX32_IP (vDot1, wAdd0, vAdd0);
		vDot1 += vBias1;
		REQUANT(1);
		PDX_SAV32_MX8_XP (vOut, pvOuta1, pvOut1, nBytesToWrite1);

     } //all A rows complete

	PDX_SAPOS_MX8_FP(pvOuta0,pvOut0); // Flush out last bytes
	PDX_SAPOS_MX8_FP(pvOuta1,pvOut1); // Flush out last bytes
    out_0 += 2*output_ch; //return the new output pointer

    return out_0;
}
