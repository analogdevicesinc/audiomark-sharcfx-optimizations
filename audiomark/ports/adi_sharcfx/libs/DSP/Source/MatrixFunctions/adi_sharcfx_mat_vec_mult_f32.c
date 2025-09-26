/* ----------------------------------------------------------------------
 * Title:        adi_sharcfx_mat_vec_mult_f32.c
 * Description:  Floating-point matrix and vector multiplication
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


#include <stdio.h>
#include <math.h>
#include <float.h>
#include <complex.h>

#include <matrix.h>
#include <filter.h>
#include <vector.h>

#define NVEC 32

#include "math_fixedpoint_vec.h"
/* Cross-platform data type definitions. */
#include "libdsp_types.h"

#include <inttypes.h>

#ifdef __XTENSA__
#include <xtensa/sim.h>
#include <xtensa/tie/xt_pdxn.h>
#endif

#include <libs/DSP/Include/dsp/matrix_functions.h>

/**
 * @brief Floating-point matrix and vector multiplication.
 * @param[in]       *pSrcMat points to the input matrix structure
 * @param[in]       *pVec points to input vector
 * @param[out]      *pDst points to output vector
 */


#define MIN(X, Y) (((X) < (Y)) ? (X) : (Y))
#define MAX(X, Y) (((X) > (Y)) ? (X) : (Y))

void adi_sharcfx_mat_vec_mult_f32(const adi_sharcfx_matrix_instance_f32 *pSrcMat, const float32_t *pVec, float32_t *pDst)
{
    uint32_t numRows = pSrcMat->numRows;
    uint32_t numCols = pSrcMat->numCols;
    const float32_t *pSrcA = pSrcMat->pData;

	xb_vecMxf32* __restrict pvOut  = (xb_vecMxf32 *)(pDst);
    valign pvOuta = PDX_LA_MXF32_PP(pvOut);
    xb_vecMxf32 inA0, inA1, inA2, inA3, inB;
    xb_vecMxf32 out, vDot0;

    for (int32_t a_rows_idx = 0; a_rows_idx < numRows; a_rows_idx +=8) // A rows from 0-64 channels
    {
    	xtfloat add_out[8] = {0}; //getting 32 output elements at a time 16 for 1st output row, 16 for 2nd output row
    	int row_count = 0;
    	//loop till all the a_rows_idx selected are complete
    	int bytes = MIN(8, numRows -  a_rows_idx);
		while(row_count< bytes)
        {
    		xb_vecMxf32 *pB = (xb_vecMxf32 *)(pVec);
    		valign pBa = PDX_LA_MXF32_PP(pB);

    	    xb_vecMxf32 *pA0 = (xb_vecMxf32 *)(pSrcA + 0*numCols);
    	    valign pA0a = PDX_LA_MXF32_PP(pA0);
    	    xb_vecMxf32 *pA1 = (xb_vecMxf32 *)(pSrcA + 1*numCols);
    	    valign pA1a = PDX_LA_MXF32_PP(pA1);
    	    xb_vecMxf32 *pA2 = (xb_vecMxf32 *)(pSrcA + 2*numCols);
    	    valign pA2a = PDX_LA_MXF32_PP(pA1);
    	    xb_vecMxf32 *pA3 = (xb_vecMxf32 *)(pSrcA + 3*numCols);
    	    valign pA3a = PDX_LA_MXF32_PP(pA3);

    	    xb_vecMxf32 a0 = 0; xb_vecMxf32 a1 = 0; xb_vecMxf32 a2 = 0; xb_vecMxf32 a3 = 0; //initialize accumulator results
    		// iterating over all columns of A, 8 at a time
    		for (int32_t a_cols_idx = 0; a_cols_idx < numCols; a_cols_idx += 8)
			{
    			PDX_LA_MXF32_XP(inA0, pA0a, pA0, 32);
    			PDX_LA_MXF32_XP(inA1, pA1a, pA1, 32);
    			PDX_LA_MXF32_XP(inA2, pA2a, pA2, 32);
    			PDX_LA_MXF32_XP(inA3, pA3a, pA3, 32);
    			PDX_LA_MXF32_XP(inB, pBa, pB,32);

    			PDX_MULA_MXF32(a0,inA0,inB);
    			PDX_MULA_MXF32(a1,inA1,inB);
    			PDX_MULA_MXF32(a2,inA2,inB);
    			PDX_MULA_MXF32(a3,inA3,inB);
            }
			//adding all the multiplied elements together to get the vector multiplication result
			add_out[row_count + 0]  += PDX_RADD_MXF32(a0) ; //8-way 32-bit Signed Add operation, reduction to scalar with saturation to int32. PDX_RADD_2MX16
			add_out[row_count + 1]  += PDX_RADD_MXF32(a1) ;
			add_out[row_count + 2]  += PDX_RADD_MXF32(a2) ;
			add_out[row_count + 3]  += PDX_RADD_MXF32(a3) ;
			pSrcA +=4*numCols; //update to next two rows
			row_count += 4; //increase row count till 16 rows
        }
    	xtfloat *add_ptr = add_out;
    	xb_vecMxf32 *vAdd0 = (xb_vecMxf32 *)(add_ptr);
    	valign wAdd0 =  PDX_LA_MXF32_PP(vAdd0);
    	PDX_LA_MXF32_IP (out, wAdd0, vAdd0);
    	PDX_SAV_MXF32_XP(out, pvOuta, pvOut, 4*bytes);
	}
	PDX_SAPOS_MXF32_FP(pvOuta,pvOut);
}
