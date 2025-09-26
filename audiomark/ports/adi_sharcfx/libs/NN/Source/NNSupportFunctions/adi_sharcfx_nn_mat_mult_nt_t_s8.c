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
 * Title:        adi_sharcfx_nn_mat_mult_s8_nt_t_s8
 * Description:  Matrix multiplication support function with the right-hand-side (rhs) matrix transposed
 * $Date:        22 April 2025
 * $Revision:    V.1.0
 * 
 * Target Processor:  SHARC-FX Processor
 * -------------------------------------------------------------------- */


#include <sys/platform.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include <complex.h>
#include <matrix.h>
#include <filter.h>
#include <libs/NN/Include/adi_sharcfx_nnsupportfunctions.h>
#include <vector.h>

#define NVEC 32

#include "math_fixedpoint_vec.h"
/* Cross-platform data type definitions. */
#include "libdsp_types.h"

//#pragma once
#include <inttypes.h>

#ifdef __XTENSA__
#include <xtensa/sim.h>
#include <xtensa/tie/xt_pdxn.h>
#endif

#define MIN(X, Y) (((X) < (Y)) ? (X) : (Y))
#define MAX(X, Y) (((X) > (Y)) ? (X) : (Y))

/*
 * transpose() function to transpose a matrix of signed 8-bit integers.
 */
void transpose(const int8_t *matrix, int8_t *transposeMatrix, int rows, int cols) {

	// Iterate through each column of the original matrix
    for (int i = 0; i < cols; i++) {
    	// Iterate through each row of the original matrix
        for (int j = 0; j < rows; j++) {
        	// Access the element at (j, i) in the original matrix and store it at (i, j) in the transposed matrix
            *(transposeMatrix + i * rows + j) = *(matrix + j * cols + i);
        }
    }
}


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
  * adi_sharcfx_nn_mat_mult_nt_t_s8() function performs matrix multiplication between two matrices, lhs and rhs, where rhs is transposed (denoted as rhs').
  * The function first computes the transpose of rhs matrix, resulting in rhs', and then multiplies matrix lhs with rhs'.
  * The result is requantized and stored in dst.
  */
adi_sharcfx_nn_status adi_sharcfx_nn_mat_mult_nt_t_s8(const int8_t *lhs,					/*[in] input LHS matrix of size [lhs_rows][rhs_cols] */
                                            const int8_t *rhs,					/*[in] input RHS matrix of size [rhs_rows][rhs_cols] */
                                            const int32_t *bias,				/*[in] bias buffer of size [rhs_rows] */
                                            int8_t *dst,						/*[out] output matrix of size [lhs_rows][rhs_cols] */
                                            const int32_t *dst_multipliers,		/*[in] multiplier buffer of size [rhs_rows] */
                                            const int32_t *dst_shifts,			/*[in] shift buffer of size [rhs_rows] */
                                            const int32_t lhs_rows,				/*[in] number of LHS rows */
                                            const int32_t rhs_rows,				/*[in] number of RHS rows */
                                            const int32_t rhs_cols,				/*[in] number of RHS cols */
                                            const int32_t lhs_offset,			/*[in] Offset to be applied to the LHS input value */
                                            const int32_t dst_offset,			/*[in] Offset to be applied the output result */
                                            const int32_t activation_min,		/*[in] Minimum value to clamp down the output. Range : int8 */
                                            const int32_t activation_max,		/*[in] Maximum value to clamp up the output. Range : int8 */
                                            const int32_t rhs_cols_offset)		/*[in] Column offset between subsequent lhs_rows */
{

    if (rhs_cols_offset < rhs_cols)
    {
        return ADI_SHARCFX_NN_ARG_ERROR;
    }
    xb_vec2Mx16 lhs0, lhs1;		/*Declaring 16 way 16 bit vector to store 2 LHS rows */
	xb_vec2Mx16 rhs0, rhs1, rhs2, rhs3;		/*Declaring 16 way 16 bit vector to store 4 RHS columns */
	xb_vec2Mx16 vInZP16 = lhs_offset;				/*Declaring 16 way 16 bit vector with InputZeroPoint */
	xb_vecMx32 vOutZeroPoint = dst_offset;			/*Declaring 8 way 32 bit vector with OutputZeroPoint */
	xb_vecMx80 vScaled, vScaled1, vScaled2, vScaled3; /*Declaring 16 way 16 bit vector to store 4 scale values in requantization */
	xb_vecMx32 vOut_ge_zero, vOut_lt_zero;			/*For shift value */
	xb_vecMx32 v32 = 32, vOne = 1; 					/*For shift value */
	vboolM en_sh_ge_zero0, en_sh_ge_zero1, en_sh_ge_zero2, en_sh_ge_zero3; /*For shift value */
	xb_vecMx32 v32_m_shift0, v32_m_shift1, v32_m_shift2, v32_m_shift3; /*For shift value */
	xb_vecMx32  vDot0, vDot1, vDot2, vDot3; 				/*Declaring 8 way 32 bit vector to store Dot product result */
	xb_vecMx32 vShift0, vShift1, vShift2, vShift3; /*Declaring 8 way 32 bit vector to store Shift values */
	xb_vecMx32 vBias0, vBias1, vBias2, vBias3;		/*Declaring 8 way 32 bit vector to store Bias values */
	xb_vecMx32 vScale0, vScale1, vScale2, vScale3;	/*Declaring 8 way 32 bit vector to store Multiplier values */
	xb_vecMx32 vOut;				/*Declaring 8 way 32 bit vector to store Output value */
	int nBytesToWrite0, nBytesToWrite1, nBytesToWrite2, nBytesToWrite3; /*For keeping count of output values to be written */

    // to store transposed RHS matrix
    int8_t transposeRHS[rhs_rows*rhs_cols];
    // transpose RHS matrix
    transpose(rhs, transposeRHS, rhs_rows, rhs_cols);


	int8_t *lhs_row_ptr = lhs;  //start with first LHS row, first LHS column
	//iterating over 2 LHS rows at a time
	for (int32_t lhs_rows_idx = 0; lhs_rows_idx <  (lhs_rows / 2) * 2; lhs_rows_idx += 2)
	{
		int nBytesLeft = rhs_cols; //reset nBytesleft for every new LHS row ---------------------------------- rhs_cols ? rhs_rows
		const int8_t *rhs_row_ptr = transposeRHS; //start with first RHS row, first RHS column
		xb_vecMx8*  pvOut0  = (xb_vecMx8 *)(dst + (0 + lhs_rows_idx)* rhs_cols );//point to first output row
		xb_vecMx8*  pvOut1  = (xb_vecMx8 *)(dst + (1 + lhs_rows_idx)* rhs_cols );//point to 2nd output row
		valign pvOuta0 = PDX_LA_MX8_PP(pvOut0); // Aligning the pointer
		valign pvOuta1 = PDX_LA_MX8_PP(pvOut1); // Aligning the pointer

		// pointers to Bias, Shift, and Scale for RHS matrix, which is reset after every 2 LHS rows
		xb_vecMx32 *pBias  = (xb_vecMx32 *)(bias);
		valign pBiasa = PDX_LA_MX32_PP (pBias);
		xb_vecMx32 *pScale  = (xb_vecMx32 *) (dst_multipliers);
		valign pScalea = PDX_LA_MX32_PP (pScale);
		xb_vecMx32 *pShift  = (xb_vecMx32 *) (dst_shifts);
		valign pShifta = PDX_LA_MX32_PP (pShift);

		// iterating over 32 RHS columns at a time
        for (int32_t rhs_cols_idx = 0; rhs_cols_idx < rhs_cols; rhs_cols_idx += 32)
        {
        	xb_vec2Mx40 a0 = 0; xb_vec2Mx40 a1 = 0; xb_vec2Mx40 a2 = 0; xb_vec2Mx40 a3 = 0; //initialize accumulator results

        	//iterating over each RHS row at a time
			for (int32_t rhs_rows_idx = 0; rhs_rows_idx < rhs_rows; rhs_rows_idx++)
			{
				int8_t *lhsp0  = (lhs_row_ptr + 0*rhs_cols_offset + rhs_rows_idx); //point to 1st LHS row
				int8_t *lhsp1  = (lhs_row_ptr + 1*rhs_cols_offset + rhs_rows_idx);//point to 2nd LHS row
				//1-way 8-bit Signed element scalar load from memory converting to 16-bit scalar, and replicating to all 16 lanes of a vector register
				PDX_LSR16_8_IP( lhs0,lhsp0, 1 );
				PDX_LSR16_8_IP( lhs1,lhsp1, 1 );
				//adding zeropoint (input offset) value to the LHS value
				lhs0 += vInZP16;
				lhs1 += vInZP16;

				//number of columns to load to RHS values (max 16 way)
				xb_vec2Mx8 *rhsp0 = (xb_vec2Mx8 *)(rhs_row_ptr + 0 + rhs_cols_idx + rhs_rows_idx*rhs_cols);
				valign rhsp0a = PDX_LA_2MX8_PP (rhsp0);//alligning vector
			    PDX_LA16_2MX8_XP(rhs0, rhsp0a, rhsp0, 16);  //RHS[0][0-15]
			    PDX_LA16_2MX8_XP(rhs1, rhsp0a, rhsp0, rhs_cols);  //RHS[0][16-32]

				// accumulate the multiplication for
				//LHS[0][i], LHS[1][i], with RHS[i][0-15], RHS[i][16-32]
				PDX_MULAW_2MX16(a0,lhs0,rhs0);
				PDX_MULAW_2MX16(a1,lhs0,rhs1);
				PDX_MULAW_2MX16(a2,lhs1,rhs0);
				PDX_MULAW_2MX16(a3,lhs1,rhs1);
			}

			//loading two 8 way 32 bit vector for bias
			PDX_LA_MX32_IP (vBias0, pBiasa, pBias);
			PDX_LA_MX32_IP (vBias1, pBiasa, pBias);
			PDX_LA_MX32_IP (vBias2, pBiasa, pBias);
			PDX_LA_MX32_IP (vBias3, pBiasa, pBias);
			//loading twp 8 way 32 bit vector for multiply
			PDX_LA_MX32_IP (vScale0, pScalea, pScale);
			PDX_LA_MX32_IP (vScale1, pScalea, pScale);
			PDX_LA_MX32_IP (vScale2, pScalea, pScale);
			PDX_LA_MX32_IP (vScale3, pScalea, pScale);
			//loading two 8 way 32 bit vector for shift
			PDX_LA_MX32_IP (vShift0, pShifta, pShift);
			PDX_LA_MX32_IP (vShift1, pShifta, pShift);
			PDX_LA_MX32_IP (vShift2, pShifta, pShift);
			PDX_LA_MX32_IP (vShift3, pShifta, pShift);
			//loading two 8 boolean vectors to separate positive and negative shifts
			en_sh_ge_zero0 = PDX_GE_MX32 (vShift0, 0);
			en_sh_ge_zero1 = PDX_GE_MX32 (vShift1, 0);
			en_sh_ge_zero2 = PDX_GE_MX32 (vShift2, 0);
			en_sh_ge_zero3 = PDX_GE_MX32 (vShift3, 0);
			//add 32bits to shift value
			v32_m_shift0 = v32 + vShift0;
			v32_m_shift1 = v32 + vShift1;
			v32_m_shift2 = v32 + vShift2;
			v32_m_shift3 = v32 + vShift3;
			// Add one to convert Q31 to Q32, use 1 if shift < 0
			vShift0 = PDX_MOV_MX32_T (vShift0 + vOne, 1, en_sh_ge_zero0);
			vShift1 = PDX_MOV_MX32_T (vShift1 + vOne, 1, en_sh_ge_zero1);
			vShift2 = PDX_MOV_MX32_T (vShift2 + vOne, 1, en_sh_ge_zero2);
			vShift3 = PDX_MOV_MX32_T (vShift3 + vOne, 1, en_sh_ge_zero3);

			// Continue iterating until processing of nBytesToWrite0 (32 columns max for 2 rows) is complete.
			//Calculate number of columns to process as nbytestowrite0
			nBytesToWrite0 = (rhs_cols-rhs_cols_idx);
			nBytesToWrite1 = nBytesToWrite0 -PDX_M;
			nBytesToWrite2 = nBytesToWrite1 -PDX_M;
			nBytesToWrite3 = nBytesToWrite2 -PDX_M;
			//for 1st LHS row* RHS row = 1st output row (8,8 columns)
			PDX_CVT32D_2MX40(vDot1, vDot0, a0 );//storing 16 way 40 bit accumulator result into two 8 way 32bit vectors
			PDX_CVT32D_2MX40(vDot3, vDot2, a1 );//storing 16 way 40 bit accumulator result into two 8 way 32bit vectors
			vDot0 += vBias0; //add bias
			REQUANT(0);
			PDX_SAV32_MX8_XP (vOut, pvOuta0, pvOut0, nBytesToWrite0);
	        vDot1 += vBias1;
			REQUANT(1);
			PDX_SAV32_MX8_XP (vOut, pvOuta0, pvOut0, nBytesToWrite1);
			vDot2 += vBias2; //add bias
			REQUANT(2);
			PDX_SAV32_MX8_XP (vOut, pvOuta0, pvOut0, nBytesToWrite2);
	        vDot3 += vBias3;
			REQUANT(3);
			PDX_SAV32_MX8_XP (vOut, pvOuta0, pvOut0, nBytesToWrite3);
			//for 2nd LHS row* RHS row = 2nd output row (8,8 columns)
			PDX_CVT32D_2MX40(vDot1, vDot0, a2 );//storing 16 way 40 bit accumulator result into two 8 way 32bit vectors
			PDX_CVT32D_2MX40(vDot3, vDot2, a3 );//storing 16 way 40 bit accumulator result into two 8 way 32bit vectors
			vDot0 += vBias0; //add bias
			REQUANT(0);
			PDX_SAV32_MX8_XP (vOut, pvOuta1, pvOut1, nBytesToWrite0);
	        vDot1 += vBias1;
			REQUANT(1);
			PDX_SAV32_MX8_XP (vOut, pvOuta1, pvOut1, nBytesToWrite1);
			vDot2 += vBias2;
			REQUANT(2);
			PDX_SAV32_MX8_XP (vOut, pvOuta1, pvOut1, nBytesToWrite2);
			vDot3 += vBias3;
			REQUANT(3);
			PDX_SAV32_MX8_XP (vOut, pvOuta1, pvOut1, nBytesToWrite3);

        } // end RHS columns loop (which was processing 32 columns at a time)
		//Flush out all the bytes for the 2 output rows
        PDX_SAPOS_MX8_FP(pvOuta0,pvOut0);
        PDX_SAPOS_MX8_FP(pvOuta1,pvOut1);
        //increment to next 2 LHS rows
        lhs_row_ptr -= 2*rhs_cols;
        lhs_row_ptr += 2*2*rhs_cols_offset;
     }

	//if LHS rows are not multiple of 2
	if(lhs_rows % 2 !=0)
	{
		// pointer for the last output row
		xb_vecMx8*  pvOut  = (xb_vecMx8 *)(dst + (lhs_rows / 2) * 2* rhs_cols );
		valign pvOuta = PDX_LA_MX8_PP(pvOut); // Aligning the pointer
		valign rhsp0a, rhsp1a, rhsp2a, rhsp3a, lhsp0a;

		const int8_t *rhs_row_ptr = transposeRHS;

		// pointers to Bias, Shift, and Scale for RHS matrix, which is reset after every 1 LHS row
		xb_vecMx32 *pBias  = (xb_vecMx32 *)(bias);
		valign pBiasa = PDX_LA_MX32_PP (pBias);
		xb_vecMx32 *pScale  = (xb_vecMx32 *) (dst_multipliers);
		valign pScalea = PDX_LA_MX32_PP (pScale);
		xb_vecMx32 *pShift  = (xb_vecMx32 *) (dst_shifts);
		valign pShifta = PDX_LA_MX32_PP (pShift);

		for (int32_t rhs_cols_idx = 0; rhs_cols_idx < rhs_cols; rhs_cols_idx += 64)
		{
			xb_vec2Mx40 a0 = 0; xb_vec2Mx40 a1 = 0; xb_vec2Mx40 a2 = 0; xb_vec2Mx40 a3 = 0; //initialize accumulator results
			//iterating over each RHS row at a time
			for (int32_t rhs_rows_idx = 0; rhs_rows_idx < rhs_rows; rhs_rows_idx++)
			{
				int8_t *lhsp0  = (lhs_row_ptr + rhs_rows_idx); //point to 1st LHS row
				PDX_LSR16_8_IP( lhs0,lhsp0, 1 ); //1-way 8-bit Signed element scalar load
				lhs0 += vInZP16; //adding zeropoint (input offset) value to the LHS value

				xb_vec2Mx8 *rhsp0 = (xb_vec2Mx8 *)(rhs_row_ptr + 0 + rhs_cols_idx + rhs_rows_idx*rhs_cols); //RHS[0][0-15]
				xb_vec2Mx8 *rhsp1 = (xb_vec2Mx8 *)(rhs_row_ptr + 16 + rhs_cols_idx + rhs_rows_idx*rhs_cols); //RHS[0][16-31]
				xb_vec2Mx8 *rhsp2 = (xb_vec2Mx8 *)(rhs_row_ptr + 32 + rhs_cols_idx + rhs_rows_idx*rhs_cols); //RHS[0][32-47]
				xb_vec2Mx8 *rhsp3 = (xb_vec2Mx8 *)(rhs_row_ptr + 48 + rhs_cols_idx + rhs_rows_idx*rhs_cols); //RHS[0][48-64]
				//number of columns to load to RHS values (max 16 way)
				LOADU_2MX8(rhs0,rhsp0,MIN(16, rhs_cols - rhs_cols_idx ));
				LOADU_2MX8(rhs1,rhsp1,MIN(16, rhs_cols - rhs_cols_idx - 16));
				LOADU_2MX8(rhs2,rhsp2,MIN(16, rhs_cols - rhs_cols_idx - 32));
				LOADU_2MX8(rhs3,rhsp3,MIN(16, rhs_cols - rhs_cols_idx - 48));
				// accumulate the multiplication for
				//LHS[0][i] with RHS[i][0-15], RHS[i][16-31], RHS[i][32-47], RHS[i][48-64]
				PDX_MULAW_2MX16(a0,rhs0,lhs0);
				PDX_MULAW_2MX16(a1,rhs1,lhs0);
				PDX_MULAW_2MX16(a2,rhs2,lhs0);
				PDX_MULAW_2MX16(a3,rhs3,lhs0);
			}

			PDX_CVT32D_2MX40(vDot1, vDot0, a0 );  //storing 16 way 40 bit accumulator result into two 8 way 32bit vectors
			PDX_CVT32D_2MX40(vDot3, vDot2, a1 );  //storing 16 way 40 bit accumulator result into two 8 way 32bit vectors

			//loading four 8 way 32 bit vector for bias
			PDX_LA_MX32_IP (vBias0, pBiasa, pBias);
			PDX_LA_MX32_IP (vBias1, pBiasa, pBias);
			PDX_LA_MX32_IP (vBias2, pBiasa, pBias);
			PDX_LA_MX32_IP (vBias3, pBiasa, pBias);
			//loading four 8 way 32 bit vector for multiply
			PDX_LA_MX32_IP (vScale0, pScalea, pScale);
			PDX_LA_MX32_IP (vScale1, pScalea, pScale);
			PDX_LA_MX32_IP (vScale2, pScalea, pScale);
			PDX_LA_MX32_IP (vScale3, pScalea, pScale);
			//loading four 8 way 32 bit vector for shift
			PDX_LA_MX32_IP (vShift0, pShifta, pShift);
			PDX_LA_MX32_IP (vShift1, pShifta, pShift);
			PDX_LA_MX32_IP (vShift2, pShifta, pShift);
			PDX_LA_MX32_IP (vShift3, pShifta, pShift);
			//loading four 8 boolean vectors to separate positive and negative shifts
			en_sh_ge_zero0 = PDX_GE_MX32 (vShift0, 0);
			en_sh_ge_zero1 = PDX_GE_MX32 (vShift1, 0);
			en_sh_ge_zero2 = PDX_GE_MX32 (vShift2, 0);
			en_sh_ge_zero3 = PDX_GE_MX32 (vShift3, 0);
			//add 32bits to shift value
			v32_m_shift0 = v32 + vShift0;
			v32_m_shift1 = v32 + vShift1;
			v32_m_shift2 = v32 + vShift2;
			v32_m_shift3 = v32 + vShift3;
			// Add one to convert Q31 to Q32, use 1 if shift < 0
			vShift0 = PDX_MOV_MX32_T (vShift0 + vOne, 1, en_sh_ge_zero0);
			vShift1 = PDX_MOV_MX32_T (vShift1 + vOne, 1, en_sh_ge_zero1);
			vShift2 = PDX_MOV_MX32_T (vShift2 + vOne, 1, en_sh_ge_zero2);
			vShift3 = PDX_MOV_MX32_T (vShift3 + vOne, 1, en_sh_ge_zero3);
			//Calculate number of columns to process as nbytestowrite
			nBytesToWrite0 = MIN(32, rhs_cols - rhs_cols_idx);
			nBytesToWrite1 = nBytesToWrite0 - PDX_M;
			nBytesToWrite2 = nBytesToWrite1 - PDX_M;
			nBytesToWrite3 = nBytesToWrite2 - PDX_M;
			//for first 8 columns
			vDot0 += vBias0;
			REQUANT(0);
			PDX_SAV32_MX8_XP (vOut, pvOuta, pvOut, nBytesToWrite0);
			//for 2nd 8 columns
			vDot1 += vBias1;
			REQUANT(1);
			PDX_SAV32_MX8_XP (vOut, pvOuta, pvOut, nBytesToWrite1);
			//for 3rd 8 columns
			vDot2 += vBias2;
			REQUANT(2);
			PDX_SAV32_MX8_XP (vOut, pvOuta, pvOut, nBytesToWrite2);
			//for 4th 8 columns
			vDot3 += vBias3;
			REQUANT(3);
			PDX_SAV32_MX8_XP (vOut, pvOuta, pvOut, nBytesToWrite3);

			if(rhs_cols-rhs_cols_idx>32){

				PDX_CVT32D_2MX40(vDot1, vDot0, a2 );  //storing 16 way 40 bit accumulator result into two 8 way 32bit vectors
				PDX_CVT32D_2MX40(vDot3, vDot2, a3 );  //storing 16 way 40 bit accumulator result into two 8 way 32bit vectors

				//loading four 8 way 32 bit vector for bias
				PDX_LA_MX32_IP (vBias0, pBiasa, pBias);
				PDX_LA_MX32_IP (vBias1, pBiasa, pBias);
				PDX_LA_MX32_IP (vBias2, pBiasa, pBias);
				PDX_LA_MX32_IP (vBias3, pBiasa, pBias);
				//loading four 8 way 32 bit vector for multiply
				PDX_LA_MX32_IP (vScale0, pScalea, pScale);
				PDX_LA_MX32_IP (vScale1, pScalea, pScale);
				PDX_LA_MX32_IP (vScale2, pScalea, pScale);
				PDX_LA_MX32_IP (vScale3, pScalea, pScale);
				//loading four 8 way 32 bit vector for shift
				PDX_LA_MX32_IP (vShift0, pShifta, pShift);
				PDX_LA_MX32_IP (vShift1, pShifta, pShift);
				PDX_LA_MX32_IP (vShift2, pShifta, pShift);
				PDX_LA_MX32_IP (vShift3, pShifta, pShift);
				//loading four 8 boolean vectors to separate positive and negative shifts
				en_sh_ge_zero0 = PDX_GE_MX32 (vShift0, 0);
				en_sh_ge_zero1 = PDX_GE_MX32 (vShift1, 0);
				en_sh_ge_zero2 = PDX_GE_MX32 (vShift2, 0);
				en_sh_ge_zero3 = PDX_GE_MX32 (vShift3, 0);
				//add 32bits to shift value
				v32_m_shift0 = v32 + vShift0;
				v32_m_shift1 = v32 + vShift1;
				v32_m_shift2 = v32 + vShift2;
				v32_m_shift3 = v32 + vShift3;
				// Add one to convert Q31 to Q32, use 1 if shift < 0
				vShift0 = PDX_MOV_MX32_T (vShift0 + vOne, 1, en_sh_ge_zero0);
				vShift1 = PDX_MOV_MX32_T (vShift1 + vOne, 1, en_sh_ge_zero1);
				vShift2 = PDX_MOV_MX32_T (vShift2 + vOne, 1, en_sh_ge_zero2);
				vShift3 = PDX_MOV_MX32_T (vShift3 + vOne, 1, en_sh_ge_zero3);

				nBytesToWrite0 = MIN(32, rhs_cols - rhs_cols_idx);
				nBytesToWrite1 = nBytesToWrite0 - PDX_M;
				nBytesToWrite2 = nBytesToWrite1 - PDX_M;
				nBytesToWrite3 = nBytesToWrite2 - PDX_M;
				//for 5th 8 columns
				vDot0 += vBias0;
				REQUANT(0);
				PDX_SAV32_MX8_XP (vOut, pvOuta, pvOut, nBytesToWrite0);
				//for 6th 8 columns
				vDot1 += vBias1;
				REQUANT(1);
				PDX_SAV32_MX8_XP (vOut, pvOuta, pvOut, nBytesToWrite1);
				//for 7th 8 columns
				vDot2 += vBias2;
				REQUANT(2);
				PDX_SAV32_MX8_XP (vOut, pvOuta, pvOut, nBytesToWrite2);
				//for 8th 8 columns
				vDot3 += vBias3;
				REQUANT(3);
				PDX_SAV32_MX8_XP (vOut, pvOuta, pvOut, nBytesToWrite3);
			}
		} // end RHS columns loop (which was processing 64 columns at a time)
	  PDX_SAPOS_MX8_FP(pvOuta,pvOut); //Flush out last bytes

	}

   return ADI_SHARCFX_NN_SUCCESS;

 }
