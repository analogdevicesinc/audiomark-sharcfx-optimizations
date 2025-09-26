
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
 * Title:        adi_sharcfx_convolve_s8.c
 * Description:  s8 version of convolution using symmetric quantization.
 * $Date:        22th April 2025
 * $Revision:    V.1.0
 *
 * Target Processor:  SHARC-FX Processor
 * -------------------------------------------------------------------- */

 #include <libs/NN/Include/adi_sharcfx_nnfunctions.h>
 #include <libs/NN/Include/adi_sharcfx_nnsupportfunctions.h>

 #include <libs/NN/Include/adi_sharcfx_nnfunctions.h>
 #include <libs/NN/Include/adi_sharcfx_nnsupportfunctions.h>


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
 /**
  *  @ingroup Public
  */

 /**
  * @addtogroup NNConv
  * @{
  */

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
         vOut = PDX_MIN_MX32 (PDX_MAX_MX32 (vOut, out_activation_min), out_activation_max)

 /*
  * Basic s8 convolution function.
  *
  * Refer header file for details. Optimal use case for the DSP/MVE implementation is when input and output channels
  * are multiples of 4 or atleast greater than 4.
  *
  */
 adi_sharcfx_nn_status adi_sharcfx_convolve_s8(const nn_context *ctx,
                                     const nn_conv_params *conv_params,
                                     const nn_per_channel_quant_params *quant_params,
                                     const nn_dims *input_dims,
                                     const int8_t *input_data,
                                     const nn_dims *filter_dims,
                                     const int8_t *filter_data,
                                     const nn_dims *bias_dims,
                                     const int32_t *bias_data,
                                     const nn_dims *output_dims,
                                     int8_t *output_data)
 {
	(void)bias_dims;

	if (ctx->buf == NULL && adi_sharcfx_convolve_s8_get_buffer_size(input_dims, filter_dims) > 0)
	{
		return ADI_SHARCFX_NN_ARG_ERROR;
	}
	int16_t *buffer_a = (int16_t *)ctx->buf;

	const int32_t input_batches = input_dims->n;
	const uint16_t input_x = input_dims->w;
	const uint16_t input_y = input_dims->h;
	const uint16_t input_ch = input_dims->c;
	const uint16_t kernel_x = filter_dims->w;
	const uint16_t kernel_y = filter_dims->h;
	const uint16_t output_x = output_dims->w;
	const uint16_t output_y = output_dims->h;
	const uint16_t output_ch = output_dims->c;

	const uint16_t pad_x = conv_params->padding.w;
	const uint16_t pad_y = conv_params->padding.h;
	const uint16_t stride_x = conv_params->stride.w;
	const uint16_t stride_y = conv_params->stride.h;

	const int32_t input_offset = conv_params->input_offset;
	const int32_t out_offset = conv_params->output_offset;
	const int32_t out_activation_min = conv_params->activation.min;
	const int32_t out_activation_max = conv_params->activation.max;
	int32_t *output_mult = quant_params->multiplier;
	int32_t *output_shift = quant_params->shift;


	int8_t transpose_filterdata[(input_ch * kernel_y * kernel_x)*output_ch];
	// transpose RHS matrix
	transpose(filter_data, transpose_filterdata, output_ch, input_ch * kernel_y * kernel_x);

	const uint16_t dilation_x = conv_params->dilation.w;
	const uint16_t dilation_y = conv_params->dilation.h;

	if(dilation_x == 1 && dilation_y==1) // processing all x values together
	{
		for ( int i_batch = 0; i_batch < input_batches; i_batch++)
		{
			int32_t i_out_y, i_out_x, i_ker_y, i_ker_x;
			int16_t *two_column_buf = buffer_a; //Generate two columns from the input tensor a GEMM computation
			int8_t *out = output_data;

			/* This part implements the im2col function */
			for (i_out_y = 0; i_out_y < output_y; i_out_y++)
			{
				for (i_out_x = 0; i_out_x < output_x; i_out_x++)
				{
					const int32_t base_idx_y = stride_y * i_out_y - pad_y;
					const int32_t base_idx_x = stride_x * i_out_x - pad_x;
					const size_t block_size = kernel_x * input_ch * sizeof(int8_t);
					memset(two_column_buf, 0, 2*block_size*kernel_y);

					for (i_ker_y = 0; i_ker_y < kernel_y; i_ker_y++)
					{
						const int32_t k_y = base_idx_y + dilation_y * i_ker_y;
						if (k_y >= 0 && k_y < input_y) // Check if the row is completely out of bounds
						{
							int32_t valid_x_start = MAX(0, -base_idx_x);
							int32_t valid_x_end = MIN(kernel_x, input_x - base_idx_x);
							// Convert valid middle block
							if (valid_x_end > valid_x_start)
							{
								const int8_t *src = input_data + (k_y * input_x + base_idx_x + valid_x_start) * input_ch;
								int16_t *dst = two_column_buf + valid_x_start * input_ch;
								size_t block_len = (valid_x_end - valid_x_start) * input_ch;

								xb_vec2Mx16 *restrict pvOut  = (xb_vec2Mx16 *)(dst);
								valign pvOuta = PDX_LA_2MX16_PP(pvOut);
								xb_vec2Mx8 *restrict pSrc  = (xb_vec2Mx8 *)(src); //B 1st row
								valign pSrca = PDX_LA_2MX8_PP (pSrc);
								xb_vec2Mx16 off = input_offset; xb_vec2Mx16 out1; xb_vec2Mx16 oSrc1;
								int loopend = (block_len/16)*16;
								for(int block_count = 0; block_count < loopend; block_count +=16)
								{
									PDX_LV16_2MX8_IP(oSrc1, pSrc, 16);
									out1 = PDX_ADDS_2MX16(off, oSrc1);
									PDX_SA_2MX16_IP(out1, pvOuta, pvOut);
								}
								int inc = block_len -loopend;
								PDX_LAV16_2MX8_XP(oSrc1, pSrca, pSrc, inc);
								out1 = PDX_ADDS_2MX16(off, oSrc1);
								PDX_SAV_2MX16_XP(out1, pvOuta, pvOut, 2*inc);
								PDX_SAPOS_2MX16_FP(pvOuta,pvOut);
							}
						}
						two_column_buf += block_size / sizeof(int8_t);  // Move the pointer by one row's worth of data
					}
					//Computation is filed for every 2 columns
					if (two_column_buf == buffer_a + 2 * input_ch * kernel_y * kernel_x)
					{
						xb_vec2Mx16 ina0, ina1, ina2, ina3, inb0, inb1;
						xb_vecMx80 vScaled, vScaled1, vScaled2, vScaled3;
						xb_vecMx32 vOut_ge_zero, vOut_lt_zero;
						xb_vecMx32 v32 = 32, vOne = 1;
						xb_vecMx32 vDot0,vShift0,vBias0,vScale0;
						xb_vecMx32 vDot1,vShift1,vBias1,vScale1;
						xb_vecMx32 vDot2,vShift2,vBias2,vScale2;
						xb_vecMx32 vDot3,vShift3,vBias3,vScale3;
						xb_vecMx32 vOutZeroPoint = out_offset;
						xb_vecMx32 a00, a01, a10, a11, a20, a21, a30, a31;
						xb_vecMx32 vOut;
						int num_col_a = input_ch * kernel_y * kernel_x;
						//first output row
						xb_vecMx8*restrict pvOut0  = (xb_vecMx8 *)out;
						valign pvOuta0 = PDX_LA_MX8_PP(pvOut0);
						//second output row
						xb_vecMx8*restrict pvOut1  = (xb_vecMx8 *)(out + output_ch );
						valign pvOuta1 = PDX_LA_MX8_PP(pvOut1);
						// pointers to Bias, Shift, and Scale for each A row/output column
						xb_vecMx32 *restrict pBias  = (xb_vecMx32 *)(bias_data);
						valign pBiasa = PDX_LA_MX32_PP (pBias);
						xb_vecMx32 *restrict pScale  = (xb_vecMx32 *) (output_mult);
						valign pScalea = PDX_LA_MX32_PP (pScale);
						xb_vecMx32 *restrict pShift  = (xb_vecMx32 *) (output_shift);
						valign pShifta = PDX_LA_MX32_PP (pShift);

						const int8_t *a_row_ptr = transpose_filterdata; // A first row
						const int16_t *b_row_ptr = buffer_a; // B first row
						for (int32_t a_rows_idx = 0; a_rows_idx < output_ch; a_rows_idx += 32) // A rows from 0-64 channels
						{
							xb_vec2Mx40 a0 = 0; xb_vec2Mx40 a1 = 0; xb_vec2Mx40 a2 = 0; xb_vec2Mx40 a3 = 0; //initialize accumulator results
							// iterating over all columns of a = 40, 16 at a time
							for (int32_t a_cols_idx = 0; a_cols_idx < num_col_a; a_cols_idx++)
							{
								int16_t *inbp0  = (b_row_ptr + a_cols_idx); //point to 1st LHS row
								int16_t *inbp1  = (b_row_ptr + a_cols_idx +  num_col_a);//point to 2nd LHS row
								//1-way 8-bit Signed element scalar load from memory converting to 16-bit scalar, and replicating to all 16 lanes of a vector register
								PDX_LSR_16_XP( inb0,inbp0,1); //B[0][0]
								PDX_LSR_16_XP( inb1,inbp1,1);//B[1][0]
								xb_vec2Mx8 *ap0 = (xb_vec2Mx8 *)(a_row_ptr + (a_rows_idx) + a_cols_idx*output_ch ); //A [0-15]
								valign ap0a = PDX_LA_2MX8_PP (ap0);
								PDX_LA16_2MX8_XP(ina0, ap0a, ap0, 16); // A[a_row_idx][0-15]
								PDX_LA16_2MX8_XP(ina1, ap0a, ap0, output_ch);// A[a_row_idx][16-32]
								// accumulate the multiplication for
								//A[all][0-15], A[all][16-32] with B[0][0], B[1][0]
								PDX_MULAW_2MX16(a0,ina0,inb0); //out1
								PDX_MULAW_2MX16(a1,ina1,inb0); //out1
								PDX_MULAW_2MX16(a2,ina0,inb1); // out2
								PDX_MULAW_2MX16(a3,ina1,inb1); //out2
							}
							// convert 40 bit to two 8 lane 32 bit vectors
							PDX_CVT32D_2MX40(a01, a00, a0);
							PDX_CVT32D_2MX40(a11, a10, a1);
							PDX_CVT32D_2MX40(a21, a20, a2);
							PDX_CVT32D_2MX40(a31, a30, a3);
							//loading two 8 way 32 bit vector for bias
							PDX_LA_MX32_IP (vBias0, pBiasa, pBias);
							PDX_LA_MX32_IP (vBias1, pBiasa, pBias);
							PDX_LA_MX32_IP (vBias2, pBiasa, pBias);
							PDX_LA_MX32_IP (vBias3, pBiasa, pBias);
							//loading two 8 way 32 bit vector for scale
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
							vboolM en_sh_ge_zero0 = PDX_GE_MX32 (vShift0, 0);
							vboolM en_sh_ge_zero1 = PDX_GE_MX32 (vShift1, 0);
							vboolM en_sh_ge_zero2 = PDX_GE_MX32 (vShift2, 0);
							vboolM en_sh_ge_zero3 = PDX_GE_MX32 (vShift3, 0);
							//add 32bits to shift value
							xb_vecMx32 v32_m_shift0 = v32 + vShift0;
							xb_vecMx32 v32_m_shift1 = v32 + vShift1;
							xb_vecMx32 v32_m_shift2 = v32 + vShift2;
							xb_vecMx32 v32_m_shift3 = v32 + vShift3;
							// Add one to convert Q31 to Q32, use 1 if shift < 0
							vShift0 = PDX_MOV_MX32_T (vShift0 + vOne, 1, en_sh_ge_zero0);
							vShift1 = PDX_MOV_MX32_T (vShift1 + vOne, 1, en_sh_ge_zero1);
							vShift2 = PDX_MOV_MX32_T (vShift2 + vOne, 1, en_sh_ge_zero2);
							vShift3 = PDX_MOV_MX32_T (vShift3 + vOne, 1, en_sh_ge_zero3);

							//Calculate number of columns to process as nbytestowrite0
							int nBytesToWrite0 = MIN(32, output_ch - a_rows_idx);
							int nBytesToWrite1 = nBytesToWrite0 - PDX_M;
							int nBytesToWrite2 = nBytesToWrite1 - PDX_M;
							int nBytesToWrite3 = nBytesToWrite2 - PDX_M;

							vDot0 = a00 + vBias0;
							REQUANT(0);
							PDX_SAV32_MX8_XP (vOut, pvOuta0, pvOut0, nBytesToWrite0);
							vDot1 = a01 + vBias1;
							REQUANT(1);
							PDX_SAV32_MX8_XP (vOut, pvOuta0, pvOut0, nBytesToWrite1);
							vDot2 = a10 + vBias2;
							REQUANT(2);
							PDX_SAV32_MX8_XP (vOut, pvOuta0, pvOut0, nBytesToWrite2);
							vDot3 = a11 + vBias3;
							REQUANT(3);
							PDX_SAV32_MX8_XP (vOut, pvOuta0, pvOut0, nBytesToWrite3);

							vDot0 = a20 + vBias0;
							REQUANT(0);
							PDX_SAV32_MX8_XP (vOut, pvOuta1, pvOut1, nBytesToWrite0);
							vDot1 = a21 + vBias1;
							REQUANT(1);
							PDX_SAV32_MX8_XP (vOut, pvOuta1, pvOut1, nBytesToWrite1);
							vDot2 = a30 + vBias2;
							REQUANT(2);
							PDX_SAV32_MX8_XP (vOut, pvOuta1, pvOut1, nBytesToWrite2);
							vDot3 = a31 + vBias3;
							REQUANT(3);
							PDX_SAV32_MX8_XP (vOut, pvOuta1, pvOut1, nBytesToWrite3);
						} //all A rows complete

						PDX_SAPOS_MX8_FP(pvOuta0,pvOut0); // Flush out last bytes
						PDX_SAPOS_MX8_FP(pvOuta1,pvOut1); // Flush out last bytes
						out += 2*output_ch; //return the new output pointer
						two_column_buf = buffer_a; // counter reset
					}
				}
			}
			//left-over because odd number of output pixels
			if (two_column_buf != buffer_a)
			{
				xb_vec2Mx16 ina0, ina1, ina2, ina3, inb0, inb1;
				xb_vecMx80 vScaled, vScaled1, vScaled2, vScaled3;
				xb_vecMx32 vOut_ge_zero, vOut_lt_zero;
				xb_vecMx32 v32 = 32, vOne = 1;
				xb_vecMx32 vDot0,vShift0,vBias0,vScale0;
				xb_vecMx32 vDot1,vShift1,vBias1,vScale1;
				xb_vecMx32 vDot2,vShift2,vBias2,vScale2;
				xb_vecMx32 vDot3,vShift3,vBias3,vScale3;
				xb_vecMx32 vOutZeroPoint = out_offset;
				xb_vecMx32 a00, a01, a10, a11, a20, a21, a30, a31;
				xb_vecMx32 vOut;

				int num_col_a = input_ch * kernel_y * kernel_x;
				//first output row
				xb_vecMx8*restrict pvOut0  = (xb_vecMx8 *)out;
				valign pvOuta0 = PDX_LA_MX8_PP(pvOut0);
				// pointers to Bias, Shift, and Scale for each A row/output column
				xb_vecMx32 *restrict pBias  = (xb_vecMx32 *)(bias_data);
				valign pBiasa = PDX_LA_MX32_PP (pBias);
				xb_vecMx32 *restrict pScale  = (xb_vecMx32 *) (output_mult);
				valign pScalea = PDX_LA_MX32_PP (pScale);
				xb_vecMx32 *restrict pShift  = (xb_vecMx32 *) (output_shift);
				valign pShifta = PDX_LA_MX32_PP (pShift);

				const int8_t *a_row_ptr = transpose_filterdata; // A first row
				const int16_t *b_row_ptr = buffer_a; // B first row

				for (int32_t a_rows_idx = 0; a_rows_idx < output_ch; a_rows_idx += 32) // A rows from 0-64 channels
				{
						xb_vec2Mx40 a0 = 0; xb_vec2Mx40 a1 = 0; xb_vec2Mx40 a2 = 0; xb_vec2Mx40 a3 = 0; //initialize accumulator results
						// iterating over all columns of a = 40, 16 at a time
						for (int32_t a_cols_idx = 0; a_cols_idx < num_col_a; a_cols_idx++)
						{
							int16_t *inbp0  = (b_row_ptr + a_cols_idx); //point to 1st LHS row
							//1-way 8-bit Signed element scalar load from memory converting to 16-bit scalar, and replicating to all 16 lanes of a vector register
							PDX_LSR_16_XP( inb0,inbp0,1);
							xb_vec2Mx8 *ap0 = (xb_vec2Mx8 *)(a_row_ptr + (a_rows_idx) + a_cols_idx*output_ch ); //A [0-15]
							valign ap0a = PDX_LA_2MX8_PP (ap0);
							PDX_LA16_2MX8_XP(ina0, ap0a, ap0, 16); // A[a_row_idx][0-15]
							PDX_LA16_2MX8_XP(ina1, ap0a, ap0, output_ch);// A[a_row_idx][16-32]
							// accumulate the multiplication for
							//A[all][0-15], A[all][16-32] with B[0][0], B[1][0]
							PDX_MULAW_2MX16(a0,ina0,inb0); //out1
							PDX_MULAW_2MX16(a1,ina1,inb0); //out1
						}
						// convert 40 bit to two 8 lane 32 bit vectors
						PDX_CVT32D_2MX40(a01, a00, a0);
						PDX_CVT32D_2MX40(a11, a10, a1);

					//loading two 8 way 32 bit vector for bias
					PDX_LA_MX32_IP (vBias0, pBiasa, pBias);
					PDX_LA_MX32_IP (vBias1, pBiasa, pBias);
					PDX_LA_MX32_IP (vBias2, pBiasa, pBias);
					PDX_LA_MX32_IP (vBias3, pBiasa, pBias);
					//loading two 8 way 32 bit vector for scale
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
					vboolM en_sh_ge_zero0 = PDX_GE_MX32 (vShift0, 0);
					vboolM en_sh_ge_zero1 = PDX_GE_MX32 (vShift1, 0);
					vboolM en_sh_ge_zero2 = PDX_GE_MX32 (vShift2, 0);
					vboolM en_sh_ge_zero3 = PDX_GE_MX32 (vShift3, 0);
					//add 32bits to shift value
					xb_vecMx32 v32_m_shift0 = v32 + vShift0;
					xb_vecMx32 v32_m_shift1 = v32 + vShift1;
					xb_vecMx32 v32_m_shift2 = v32 + vShift2;
					xb_vecMx32 v32_m_shift3 = v32 + vShift3;
					// Add one to convert Q31 to Q32, use 1 if shift < 0
					vShift0 = PDX_MOV_MX32_T (vShift0 + vOne, 1, en_sh_ge_zero0);
					vShift1 = PDX_MOV_MX32_T (vShift1 + vOne, 1, en_sh_ge_zero1);
					vShift2 = PDX_MOV_MX32_T (vShift2 + vOne, 1, en_sh_ge_zero2);
					vShift3 = PDX_MOV_MX32_T (vShift3 + vOne, 1, en_sh_ge_zero3);

					//Calculate number of columns to process as nbytestowrite0
					int nBytesToWrite0 = MIN(32, output_ch - a_rows_idx);
					int nBytesToWrite1 = nBytesToWrite0 - PDX_M;
					int nBytesToWrite2 = nBytesToWrite1 - PDX_M;
					int nBytesToWrite3 = nBytesToWrite2 - PDX_M;

					vDot0 = a00 + vBias0;
					REQUANT(0);
					PDX_SAV32_MX8_XP (vOut, pvOuta0, pvOut0, nBytesToWrite0);
					vDot1 = a01 + vBias1;
					REQUANT(1);
					PDX_SAV32_MX8_XP (vOut, pvOuta0, pvOut0, nBytesToWrite1);
					vDot2 = a10 + vBias2;
					REQUANT(2);
					PDX_SAV32_MX8_XP (vOut, pvOuta0, pvOut0, nBytesToWrite2);
					vDot3 = a11 + vBias3;
					REQUANT(3);
					PDX_SAV32_MX8_XP (vOut, pvOuta0, pvOut0, nBytesToWrite3);
				} //all A rows complete
				PDX_SAPOS_MX8_FP(pvOuta0,pvOut0); // Flush out last bytes
				out += output_ch; //return the new output pointer
			}
			input_data += (input_x * input_y * input_ch); // Advance to the next batch
			output_data += (output_x * output_y * output_ch);
		}
	} //dilation == 1
	else
	{
		for ( int i_batch = 0; i_batch < input_batches; i_batch++)
		{
			int32_t i_out_y, i_out_x, i_ker_y, i_ker_x;	
			int16_t *two_column_buf = buffer_a; //Generate two columns from the input tensor a GEMM computation
			int8_t *out = output_data;
			//This part implements the im2col function
			for (i_out_y = 0; i_out_y < output_y; i_out_y++) 
			{
				for (i_out_x = 0; i_out_x < output_x; i_out_x++)
				{
					const int32_t base_idx_y = stride_y * i_out_y - pad_y;
					const int32_t base_idx_x = stride_x * i_out_x - pad_x;
					const size_t block_size = kernel_x * input_ch * sizeof(int8_t);
					for (i_ker_y = 0; i_ker_y < kernel_y; i_ker_y++)
					{
						const int32_t k_y = base_idx_y + dilation_y * i_ker_y;
						if (k_y < 0 || k_y >= input_y)// Check if the row is completely out of bounds
						{
							memset(two_column_buf, 0, block_size*sizeof(int16_t)/sizeof(int8_t)); // Zero out the entire block if out of bounds
						}
						else
						{
							xb_vec2Mx16 *restrict pvOut  = (xb_vec2Mx16 *)two_column_buf;
							valign pvOuta = PDX_LA_2MX16_PP(pvOut);
							xb_vec2Mx8 *restrict pSrc  = (xb_vec2Mx8 *)(input_data + (k_y * input_x + base_idx_x) * input_ch); //B 1st row
							valign pSrca = PDX_LA_2MX8_PP (pSrc);
							xb_vec2Mx16 off = input_offset; xb_vec2Mx16 out1; xb_vec2Mx16 out2; xb_vec2Mx16 oSrc1;  xb_vec2Mx16 oSrc2;
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
							// Now zero out specific positions within the copied row if needed
							for (i_ker_x = 0; i_ker_x < kernel_x; i_ker_x++)
							{
								const int32_t k_x = base_idx_x + dilation_x * i_ker_x;
								if (k_x < 0 || k_x >= input_x)
								{
									memset(two_column_buf + i_ker_x * input_ch, 0, sizeof(int16_t) * input_ch); // Zero out the pixel data at the out-of-bounds position
								}
							}
						}
						two_column_buf += block_size / sizeof(int8_t);  // Move the pointer by one row's worth of data
					}
					// Computation is filed for every 2 columns
					if (two_column_buf == buffer_a + 2 * input_ch * kernel_y * kernel_x)
					{
						xb_vec2Mx16 ina0, ina1, ina2, ina3, inb0, inb1;
						xb_vecMx80 vScaled, vScaled1, vScaled2, vScaled3;
						xb_vecMx32 vOut_ge_zero, vOut_lt_zero;
						xb_vecMx32 v32 = 32, vOne = 1;
						xb_vecMx32 vDot0,vShift0,vBias0,vScale0;
						xb_vecMx32 vDot1,vShift1,vBias1,vScale1;
						xb_vecMx32 vDot2,vShift2,vBias2,vScale2;
						xb_vecMx32 vDot3,vShift3,vBias3,vScale3;
						xb_vecMx32 vOutZeroPoint = out_offset;
						xb_vecMx32 a00, a01, a10, a11, a20, a21, a30, a31;
						xb_vecMx32 vOut;
						int num_col_a = input_ch * kernel_y * kernel_x;
						//first output row
						xb_vecMx8*restrict pvOut0  = (xb_vecMx8 *)out;
						valign pvOuta0 = PDX_LA_MX8_PP(pvOut0);
						//second output row
						xb_vecMx8*restrict pvOut1  = (xb_vecMx8 *)(out + output_ch );
						valign pvOuta1 = PDX_LA_MX8_PP(pvOut1);
						// pointers to Bias, Shift, and Scale for each A row/output column
						xb_vecMx32 *restrict pBias  = (xb_vecMx32 *)(bias_data);
						valign pBiasa = PDX_LA_MX32_PP (pBias);
						xb_vecMx32 *restrict pScale  = (xb_vecMx32 *) (output_mult);
						valign pScalea = PDX_LA_MX32_PP (pScale);
						xb_vecMx32 *restrict pShift  = (xb_vecMx32 *) (output_shift);
						valign pShifta = PDX_LA_MX32_PP (pShift);

						const int8_t *a_row_ptr = transpose_filterdata; // A first row
						const int16_t *b_row_ptr = buffer_a; // B first row

						for (int32_t a_rows_idx = 0; a_rows_idx < output_ch; a_rows_idx += 32) // A rows from 0-64 channels
						{
							xb_vec2Mx40 a0 = 0; xb_vec2Mx40 a1 = 0; xb_vec2Mx40 a2 = 0; xb_vec2Mx40 a3 = 0; //initialize accumulator results
							// iterating over all columns of a = 40, 16 at a time
							for (int32_t a_cols_idx = 0; a_cols_idx < num_col_a; a_cols_idx++)
							{
								int16_t *inbp0  = (b_row_ptr + a_cols_idx); //point to 1st LHS row
								int16_t *inbp1  = (b_row_ptr + a_cols_idx +  num_col_a);//point to 2nd LHS row
								//1-way 8-bit Signed element scalar load from memory converting to 16-bit scalar, and replicating to all 16 lanes of a vector register
								PDX_LSR_16_XP( inb0,inbp0,1);
								PDX_LSR_16_XP( inb1,inbp1,1);
								xb_vec2Mx8 *ap0 = (xb_vec2Mx8 *)(a_row_ptr + (a_rows_idx) + a_cols_idx*output_ch ); //A [0-15]
								valign ap0a = PDX_LA_2MX8_PP (ap0);
								PDX_LA16_2MX8_XP(ina0, ap0a, ap0, 16); // A[a_row_idx][0-15]
								PDX_LA16_2MX8_XP(ina1, ap0a, ap0, output_ch);// A[a_row_idx][16-32]
								// accumulate the multiplication for
								//A[all][0-15], A[all][16-32] with B[0][0], B[1][0]
								PDX_MULAW_2MX16(a0,ina0,inb0); //out1
								PDX_MULAW_2MX16(a1,ina1,inb0); //out1
								PDX_MULAW_2MX16(a2,ina0,inb1); // out2
								PDX_MULAW_2MX16(a3,ina1,inb1); //out2
							}
							// convert 40 bit to two 8 lane 32 bit vectors
							PDX_CVT32D_2MX40(a01, a00, a0);
							PDX_CVT32D_2MX40(a11, a10, a1);
							PDX_CVT32D_2MX40(a21, a20, a2);
							PDX_CVT32D_2MX40(a31, a30, a3);
							//loading two 8 way 32 bit vector for bias
							PDX_LA_MX32_IP (vBias0, pBiasa, pBias);
							PDX_LA_MX32_IP (vBias1, pBiasa, pBias);
							PDX_LA_MX32_IP (vBias2, pBiasa, pBias);
							PDX_LA_MX32_IP (vBias3, pBiasa, pBias);
							//loading two 8 way 32 bit vector for scale
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
							vboolM en_sh_ge_zero0 = PDX_GE_MX32 (vShift0, 0);
							vboolM en_sh_ge_zero1 = PDX_GE_MX32 (vShift1, 0);
							vboolM en_sh_ge_zero2 = PDX_GE_MX32 (vShift2, 0);
							vboolM en_sh_ge_zero3 = PDX_GE_MX32 (vShift3, 0);
							//add 32bits to shift value
							xb_vecMx32 v32_m_shift0 = v32 + vShift0;
							xb_vecMx32 v32_m_shift1 = v32 + vShift1;
							xb_vecMx32 v32_m_shift2 = v32 + vShift2;
							xb_vecMx32 v32_m_shift3 = v32 + vShift3;
							// Add one to convert Q31 to Q32, use 1 if shift < 0
							vShift0 = PDX_MOV_MX32_T (vShift0 + vOne, 1, en_sh_ge_zero0);
							vShift1 = PDX_MOV_MX32_T (vShift1 + vOne, 1, en_sh_ge_zero1);
							vShift2 = PDX_MOV_MX32_T (vShift2 + vOne, 1, en_sh_ge_zero2);
							vShift3 = PDX_MOV_MX32_T (vShift3 + vOne, 1, en_sh_ge_zero3);
							//Calculate number of columns to process as nbytestowrite0
							int nBytesToWrite0 = MIN(32, output_ch - a_rows_idx);
							int nBytesToWrite1 = nBytesToWrite0 - PDX_M;
							int nBytesToWrite2 = nBytesToWrite1 - PDX_M;
							int nBytesToWrite3 = nBytesToWrite2 - PDX_M;
							vDot0 = a00 + vBias0;
							REQUANT(0);
							PDX_SAV32_MX8_XP (vOut, pvOuta0, pvOut0, nBytesToWrite0);
							vDot1 = a01 + vBias1;
							REQUANT(1);
							PDX_SAV32_MX8_XP (vOut, pvOuta0, pvOut0, nBytesToWrite1);
							vDot2 = a10 + vBias2;
							REQUANT(2);
							PDX_SAV32_MX8_XP (vOut, pvOuta0, pvOut0, nBytesToWrite2);
							vDot3 = a11 + vBias3;
							REQUANT(3);
							PDX_SAV32_MX8_XP (vOut, pvOuta0, pvOut0, nBytesToWrite3);

							vDot0 = a20 + vBias0;
							REQUANT(0);
							PDX_SAV32_MX8_XP (vOut, pvOuta1, pvOut1, nBytesToWrite0);
							vDot1 = a21 + vBias1;
							REQUANT(1);
							PDX_SAV32_MX8_XP (vOut, pvOuta1, pvOut1, nBytesToWrite1);
							vDot2 = a30 + vBias2;
							REQUANT(2);
							PDX_SAV32_MX8_XP (vOut, pvOuta1, pvOut1, nBytesToWrite2);
							vDot3 = a31 + vBias3;
							REQUANT(3);
							PDX_SAV32_MX8_XP (vOut, pvOuta1, pvOut1, nBytesToWrite3);
						} //all A rows complete
						PDX_SAPOS_MX8_FP(pvOuta0,pvOut0); // Flush out last bytes
						PDX_SAPOS_MX8_FP(pvOuta1,pvOut1); // Flush out last bytes
						out += 2*output_ch; //return the new output pointer
						// counter reset
						two_column_buf = buffer_a;
					}
				}
			}

			// left-over because odd number of output pixels
			if (two_column_buf != buffer_a)
			{
				xb_vec2Mx16 ina0, ina1, ina2, ina3, inb0, inb1;
				xb_vecMx80 vScaled, vScaled1, vScaled2, vScaled3;
				xb_vecMx32 vOut_ge_zero, vOut_lt_zero;
				xb_vecMx32 v32 = 32, vOne = 1;
				xb_vecMx32 vDot0,vShift0,vBias0,vScale0;
				xb_vecMx32 vDot1,vShift1,vBias1,vScale1;
				xb_vecMx32 vDot2,vShift2,vBias2,vScale2;
				xb_vecMx32 vDot3,vShift3,vBias3,vScale3;
				xb_vecMx32 vOutZeroPoint = out_offset;
				xb_vecMx32 a00, a01, a10, a11, a20, a21, a30, a31;
				xb_vecMx32 vOut;

				int num_col_a = input_ch * kernel_y * kernel_x;
				xb_vecMx8*restrict pvOut0  = (xb_vecMx8 *)out;
				valign pvOuta0 = PDX_LA_MX8_PP(pvOut0);
				// pointers to Bias, Shift, and Scale for each A row/output column
				xb_vecMx32 *restrict pBias  = (xb_vecMx32 *)(bias_data);
				valign pBiasa = PDX_LA_MX32_PP (pBias);
				xb_vecMx32 *restrict pScale  = (xb_vecMx32 *) (output_mult);
				valign pScalea = PDX_LA_MX32_PP (pScale);
				xb_vecMx32 *restrict pShift  = (xb_vecMx32 *) (output_shift);
				valign pShifta = PDX_LA_MX32_PP (pShift);
				const int8_t *a_row_ptr = transpose_filterdata; // A first row
				const int16_t *b_row_ptr = buffer_a; // B first row
				for (int32_t a_rows_idx = 0; a_rows_idx < output_ch; a_rows_idx += 32) // A rows from 0-64 channels
				{
						xb_vec2Mx40 a0 = 0; xb_vec2Mx40 a1 = 0; xb_vec2Mx40 a2 = 0; xb_vec2Mx40 a3 = 0; //initialize accumulator results
						// iterating over all columns of a = 40, 16 at a time
						for (int32_t a_cols_idx = 0; a_cols_idx < num_col_a; a_cols_idx++)
						{
							int16_t *inbp0  = (b_row_ptr + a_cols_idx); //point to 1st LHS row
							//1-way 8-bit Signed element scalar load from memory converting to 16-bit scalar, and replicating to all 16 lanes of a vector register
							PDX_LSR_16_XP( inb0,inbp0,1);
							xb_vec2Mx8 *ap0 = (xb_vec2Mx8 *)(a_row_ptr + (a_rows_idx) + a_cols_idx*output_ch ); //A [0-15]
							valign ap0a = PDX_LA_2MX8_PP (ap0);
							PDX_LA16_2MX8_XP(ina0, ap0a, ap0, 16); // A[a_row_idx][0-15]
							PDX_LA16_2MX8_XP(ina1, ap0a, ap0, output_ch);// A[a_row_idx][16-32]
							// accumulate the multiplication for
							//A[all][0-15], A[all][16-32] with B[0][0], B[1][0]
							PDX_MULAW_2MX16(a0,ina0,inb0); //out1
							PDX_MULAW_2MX16(a1,ina1,inb0); //out1
						}
						// convert 40 bit to two 8 lane 32 bit vectors
						PDX_CVT32D_2MX40(a01, a00, a0);
						PDX_CVT32D_2MX40(a11, a10, a1);
					//loading two 8 way 32 bit vector for bias
					PDX_LA_MX32_IP (vBias0, pBiasa, pBias);
					PDX_LA_MX32_IP (vBias1, pBiasa, pBias);
					PDX_LA_MX32_IP (vBias2, pBiasa, pBias);
					PDX_LA_MX32_IP (vBias3, pBiasa, pBias);
					//loading two 8 way 32 bit vector for scale
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
					vboolM en_sh_ge_zero0 = PDX_GE_MX32 (vShift0, 0);
					vboolM en_sh_ge_zero1 = PDX_GE_MX32 (vShift1, 0);
					vboolM en_sh_ge_zero2 = PDX_GE_MX32 (vShift2, 0);
					vboolM en_sh_ge_zero3 = PDX_GE_MX32 (vShift3, 0);
					//add 32bits to shift value
					xb_vecMx32 v32_m_shift0 = v32 + vShift0;
					xb_vecMx32 v32_m_shift1 = v32 + vShift1;
					xb_vecMx32 v32_m_shift2 = v32 + vShift2;
					xb_vecMx32 v32_m_shift3 = v32 + vShift3;
					// Add one to convert Q31 to Q32, use 1 if shift < 0
					vShift0 = PDX_MOV_MX32_T (vShift0 + vOne, 1, en_sh_ge_zero0);
					vShift1 = PDX_MOV_MX32_T (vShift1 + vOne, 1, en_sh_ge_zero1);
					vShift2 = PDX_MOV_MX32_T (vShift2 + vOne, 1, en_sh_ge_zero2);
					vShift3 = PDX_MOV_MX32_T (vShift3 + vOne, 1, en_sh_ge_zero3);
					//Calculate number of columns to process as nbytestowrite0
					int nBytesToWrite0 = MIN(32, output_ch - a_rows_idx);
					int nBytesToWrite1 = nBytesToWrite0 - PDX_M;
					int nBytesToWrite2 = nBytesToWrite1 - PDX_M;
					int nBytesToWrite3 = nBytesToWrite2 - PDX_M;
					vDot0 = a00 + vBias0;
					REQUANT(0);
					PDX_SAV32_MX8_XP (vOut, pvOuta0, pvOut0, nBytesToWrite0);
					vDot1 = a01 + vBias1;
					REQUANT(1);
					PDX_SAV32_MX8_XP (vOut, pvOuta0, pvOut0, nBytesToWrite1);
					vDot2 = a10 + vBias2;
					REQUANT(2);
					PDX_SAV32_MX8_XP (vOut, pvOuta0, pvOut0, nBytesToWrite2);
					vDot3 = a11 + vBias3;
					REQUANT(3);
					PDX_SAV32_MX8_XP (vOut, pvOuta0, pvOut0, nBytesToWrite3);
				} //all A rows complete
				PDX_SAPOS_MX8_FP(pvOuta0,pvOut0); // Flush out last bytes
				out += output_ch; //return the new output pointer
			}
			/// Advance to the next batch 
			input_data += (input_x * input_y * input_ch);
			output_data += (output_x * output_y * output_ch);
		}
	}
	// Return to application 
	return ADI_SHARCFX_NN_SUCCESS;
 }
 
 int32_t adi_sharcfx_convolve_s8_get_buffer_size(const nn_dims *input_dims, const nn_dims *filter_dims)
 {

     return (2 * input_dims->c * filter_dims->w * filter_dims->h) * (int32_t)sizeof(int16_t);

 }

 /**
  * @} end of NNConv group
  */

