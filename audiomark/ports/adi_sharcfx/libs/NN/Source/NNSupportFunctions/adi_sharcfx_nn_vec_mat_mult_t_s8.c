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
 * Title:        adi_sharcfx_nn_vec_mat_mult_t_s8
 * Description:  s8 vector by matrix (transposed) multiplication
 * $Date:        22 April 2025
 * $Revision:    V.1.0
 *
 * Target Processor: SHARC-FX Processor
 *
 * -------------------------------------------------------------------- */


#include "libdsp_types.h"

#include <inttypes.h>
#include <libs/NN/Include/adi_sharcfx_nnsupportfunctions.h>

#ifdef __XTENSA__
#include <xtensa/sim.h>
#include <xtensa/tie/xt_pdxn.h>
#endif

#include <stdlib.h>  // For malloc and free
#include <string.h>  // For memset

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



adi_sharcfx_nn_status adi_sharcfx_nn_vec_mat_mult_t_s8(const int8_t *lhs,
                                             const int8_t *rhs,
                                             const int32_t *bias,
                                             int8_t *dst,
                                             const int32_t lhs_offset,
                                             const int32_t dst_offset,
                                             const int32_t dst_multiplier,
                                             const int32_t dst_shift,
                                             const int32_t rhs_cols,
                                             const int32_t rhs_rows,
                                             const int32_t activation_min,
                                             const int32_t activation_max,
                                             const int32_t address_offset)
{

    uint32_t numRows = rhs_rows;
    uint32_t numCols = rhs_cols;
    const int8_t *restrict prhs = rhs;
    const xb_vecMx8* __restrict pvOut  = (const xb_vecMx8 *)(dst);
    valign pvOuta = PDX_LA_MX8_PP(pvOut); // Aligning the pointer

	xb_vecMx80 vScaled, vScaled1, vScaled2, vScaled3;
	xb_vecMx32 vOut_ge_zero, vOut_lt_zero;
	xb_vecMx32 v32 = 32, vOne = 1;
	xb_vecMx32 vDot0,vShift0,vBias0,vScale0;
	xb_vecMx32 vDot1,vShift1,vBias1,vScale1;
    xb_vec2Mx16 rhs0, rhs1, rhs2, rhs3, lhs0;
    xb_vecMx32 vOut;
    xb_vec2Mx16 vInZeroPoint = lhs_offset;
    xb_vecMx32 vOutZeroPoint = dst_offset;
    xb_vecMx32 a00, a01, a10, a11, a20, a21, a30, a31;

	xb_vecMx32 *restrict pBias  = (xb_vecMx32 *)(bias);
	valign pBiasa = PDX_LA_MX32_PP (pBias);

    for (int32_t rhs_rows_idx = 0; rhs_rows_idx < numRows; rhs_rows_idx +=8) // A rows from 0-64 channels
    {
    	xb_int32 *add_out = (xb_int32 *)malloc(8*address_offset * sizeof(int));
    	memset(add_out, 0, 8*address_offset * sizeof(int));
    	int row_count = 0;
    	//loop till all the a_rows_idx selected are complete
    	while(row_count< MIN(8, numRows -  rhs_rows_idx))
        {
    		xb_vec2Mx8 *plhs = (xb_vec2Mx8 *)(lhs);
			valign plhsa = PDX_LA_2MX8_PP (plhs);//alligning vector

			xb_vec2Mx8 *prhs0 = (xb_vec2Mx8 *)(prhs);
			valign prhs0a = PDX_LA_2MX8_PP (prhs0);//alligning vector
			xb_vec2Mx8 *prhs1 = (xb_vec2Mx8 *)(prhs +  1*numCols);
			valign prhs1a = PDX_LA_2MX8_PP (prhs1);//alligning vector
			xb_vec2Mx8 *prhs2 = (xb_vec2Mx8 *)(prhs +  2*numCols);
			valign prhs2a = PDX_LA_2MX8_PP (prhs2);//alligning vector
			xb_vec2Mx8 *prhs3 = (xb_vec2Mx8 *)(prhs +  3*numCols);
			valign prhs3a = PDX_LA_2MX8_PP (prhs3);//alligning vector


			xb_vec2Mx40 a0 = 0; xb_vec2Mx40 a1 = 0; xb_vec2Mx40 a2 = 0; xb_vec2Mx40 a3 = 0; //initialize accumulator results
    		// iterating over all columns of A, 8 at a time
    		for (int32_t a_cols_idx = 0; a_cols_idx < numCols; a_cols_idx += 16)
			{
    			PDX_LA16_2MX8_XP(rhs0, prhs0a, prhs0, 16);  //LHS[0][0-15]
    			PDX_LA16_2MX8_XP(rhs1, prhs1a, prhs1, 16);  //LHS[0][0-15]
    			PDX_LA16_2MX8_XP(rhs2, prhs2a, prhs2, 16);  //LHS[0][0-15]
    			PDX_LA16_2MX8_XP(rhs3, prhs3a, prhs3, 16);  //LHS[0][0-15]
    			PDX_LA16_2MX8_XP(lhs0, plhsa, plhs, 16);  //LHS[0][0-15]
    			lhs0 += vInZeroPoint;
    			PDX_MULAW_2MX16(a0,lhs0,rhs0);
				PDX_MULAW_2MX16(a1,lhs0,rhs1);
				PDX_MULAW_2MX16(a2,lhs0,rhs2);
				PDX_MULAW_2MX16(a3,lhs0,rhs3);

            }
			//adding all the multiplied elements together to get the vector multiplication result
			//output row
    		PDX_CVT32D_2MX40(a00, a01, a0);
			PDX_CVT32D_2MX40(a10, a11, a1);
			PDX_CVT32D_2MX40(a20, a21, a2);
			PDX_CVT32D_2MX40(a30, a31, a3);
			//adding all the multiplied elements together to get the vector multiplication result

			//first output row
			add_out[row_count + 0*address_offset ]  += PDX_RADDS_MX32(a00) + PDX_RADDS_MX32(a01); //8-way 32-bit Signed Add operation, reduction to scalar with saturation to int32. PDX_RADD_2MX16
			add_out[row_count + 1*address_offset]  += PDX_RADDS_MX32(a10) + PDX_RADDS_MX32(a11);
			add_out[row_count + 2*address_offset]  += PDX_RADDS_MX32(a20) + PDX_RADDS_MX32(a21);
			add_out[row_count + 3*address_offset]  += PDX_RADDS_MX32(a30) + PDX_RADDS_MX32(a31);

			prhs +=4*numCols; //update to next two rows
			row_count += 4; //increase row count till 16 rows
        }

		// Bias, Shift, and Scale for RHS matrix
		xb_vecMx32 vScale0  = dst_multiplier;
		xb_vecMx32 vShift0  = dst_shift;
    	//loading two 8 way 32 bit vector for bias
		PDX_LA_MX32_IP (vBias0, pBiasa, pBias);
		//loading two 8 boolean vectors to separate positive and negative shifts
		vboolM en_sh_ge_zero0 = PDX_GE_MX32 (vShift0, 0);
		//add 32bits to shift value
		xb_vecMx32 v32_m_shift0 = v32 + vShift0;
		// Add one to convert Q31 to Q32, use 1 if shift < 0
		vShift0 = PDX_MOV_MX32_T (vShift0 + vOne, 1, en_sh_ge_zero0);
		//Calculate number of columns to process as nbytestowrite0
		int nBytesToWrite0 = MIN(8, numRows - rhs_rows_idx);
		//vector pointer to the addition result obtained

		xb_int32 *add_ptr = add_out;
		xb_vecMx32 *vAdd0 = (xb_vecMx32 *)(add_ptr);
		valign wAdd0 =  PDX_LA_MX32_PP(vAdd0);
		//loading first 8 way 32 bit vector output first row
		PDX_LA_MX32_IP (vDot0, wAdd0, vAdd0);
		vDot0 += vBias0;
		REQUANT(0);
		PDX_SAV32_MX8_XP(vOut, pvOuta, pvOut, nBytesToWrite0);
		free(add_out);
    }

	PDX_SAPOS_MX8_FP(pvOuta,pvOut);

	return ADI_SHARCFX_NN_SUCCESS;
}

