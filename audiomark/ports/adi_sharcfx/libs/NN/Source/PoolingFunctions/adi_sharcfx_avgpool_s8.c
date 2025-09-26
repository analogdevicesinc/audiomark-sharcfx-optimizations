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
 * Title:        adi_sharcfx_avgpool_s8.c
 * Description:  Pooling function implementations
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
#include <libs/NN/Include/adi_sharcfx_nnfunctions.h>
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

/*
 * s8 average pooling function
 *
 * Refer to header file for details.
 *
 */
adi_sharcfx_nn_status adi_sharcfx_avgpool_s8(const nn_context *ctx,
                                   const nn_pool_params *pool_params,
                                   const nn_dims *input_dims,
                                   const int8_t *src,
                                   const nn_dims *filter_dims,
                                   const nn_dims *output_dims,
                                   int8_t *dst)
{
	(void)ctx;
	    const int32_t input_y = input_dims->h;
	    const int32_t input_x = input_dims->w;
	    const int32_t output_y = output_dims->h;
	    const int32_t output_x = output_dims->w;
	    const int32_t stride_y = pool_params->stride.h;
	    const int32_t stride_x = pool_params->stride.w;
	    const int32_t kernel_y = filter_dims->h;
	    const int32_t kernel_x = filter_dims->w;
	    const int32_t pad_y = pool_params->padding.h;
	    const int32_t pad_x = pool_params->padding.w;
	    const int32_t act_min = pool_params->activation.min;
	    const int32_t act_max = pool_params->activation.max;
	    const int32_t ch_src = input_dims->c;

	    if (ctx->buf == NULL && adi_sharcfx_avgpool_s8_get_buffer_size(output_dims->w, input_dims->c))
	    {
	        return ADI_SHARCFX_NN_ARG_ERROR;
	    }

	    int32_t *buffer = (int32_t *)ctx->buf;
	    (void)buffer;

	    xb_vecMx32 sum_ch00, sum_ch01, sum_ch10, sum_ch11;
	    xb_vec2Mx16 in_vec0, in_vec1, in_vec2, in_vec3;
	    xb_vecMx32 vDot0, vDot1, rem;
	    xb_vec2Mx16 in_zero;

	    for (int i_y = 0; i_y < output_y; i_y++)
	    {
	        for (int i_x = 0; i_x < output_x; i_x++)
	        {
	        	//kernel X and Y axis starting and end points
	            const int32_t k_y_start = MAX(0, i_y * stride_y - pad_y);
	            const int32_t k_y_end = MIN(i_y * stride_y - pad_y + kernel_y, input_y);
	            const int32_t k_x_start = MAX(0, i_x * stride_x - pad_x);
	            const int32_t k_x_end = MIN(i_x * stride_x - pad_x + kernel_x, input_x);

	            //input and output
	            const int8_t *src_in = src;
	            xb_vecMx8* __restrict pvOut  = (xb_vecMx8 *)(dst + ch_src * (i_x + i_y * output_x));
	            valign pvOuta = PDX_LA_MX8_PP(pvOut);

	            //loop over every 32 channels
		        for (int ch_idx = 0; ch_idx < ch_src; ch_idx+=32) // src_in+32
		        {
		        	int32_t count = 0; // counter for averaging at the end
		        	xb_vec2Mx40 out_ch0 = 0; xb_vec2Mx40 out_ch1 = 0; xb_vec2Mx40 out_ch2 = 0; xb_vec2Mx40 out_ch3 = 0;
		        	xb_vec2Mx16 in_zero = 0; //adding zero value for add and accumulate  "PDX_ADDAW_2MX16"
		        	xb_vecMx32 out0 = 0; xb_vecMx32 out1 = 0; xb_vecMx32 out2 = 0; xb_vecMx32 out3 = 0;

	                for (int k_y = k_y_start; k_y < k_y_end; k_y++)
	                {
	                    for (int k_x = k_x_start; k_x < k_x_end; k_x++)
	                    {
	                        xb_vec2Mx8 *in_ptr = (xb_vec2Mx8 *)(src_in + (ch_src * (k_x + k_y * input_x)) ); //src_in modified at end
	                        valign in_ptra = PDX_LA_2MX8_PP (in_ptr);
	                        PDX_LA16_2MX8_XP(in_vec0, in_ptra, in_ptr, 16);  //src_in[0][0-15]
	                        PDX_LA16_2MX8_XP(in_vec1, in_ptra, in_ptr, ch_src);  //src_in[0][16-32]

	                        // Add and accumulate the src_in values to out_ch0, out_ch1.
	                        PDX_ADDAW_2MX16(out_ch0, in_vec0, in_zero); //out_ch0 += in_vec0 + in_zero
	                        PDX_ADDAW_2MX16(out_ch1, in_vec1, in_zero); //out_ch1 += in_vec1 + in_zero
	                        count++;
	                    }
	                }

	                if (count == 0)
					{
						return ADI_SHARCFX_NN_ARG_ERROR;
					}

	                //storing 16 way 40 bit accumulator result into two 8 way 32bit vectors
	                PDX_CVT32D_2MX40(sum_ch01, sum_ch00, out_ch0 );
	                PDX_CVT32D_2MX40(sum_ch11, sum_ch10, out_ch1 );

					// performing this condition
	                //sum = sum > 0 ? (sum + count / 2) / count : (sum - count / 2) / count;
	                xb_vecMx32 half_count = count / 2;
	                //greater than and less than condition for 32 channels 00, 01, 10, 11
					vboolM en_sum_ge_zero00 = PDX_GE_MX32 (sum_ch00, 0);
					vboolM en_sum_lt_zero00 = PDX_LT_MX32 (sum_ch00, 0);
					vboolM en_sum_ge_zero01 = PDX_GE_MX32 (sum_ch01, 0);
					vboolM en_sum_lt_zero01 = PDX_LT_MX32 (sum_ch01, 0);
					vboolM en_sum_ge_zero10 = PDX_GE_MX32 (sum_ch10, 0);
					vboolM en_sum_lt_zero10 = PDX_LT_MX32 (sum_ch10, 0);
					vboolM en_sum_ge_zero11 = PDX_GE_MX32 (sum_ch11, 0);
					vboolM en_sum_lt_zero11 = PDX_LT_MX32 (sum_ch11, 0);

					PDX_ADDS_MX32_T(out0, sum_ch00, half_count, en_sum_ge_zero00 ); //(sum + count / 2) / count
					PDX_SUBS_MX32_T(out0, sum_ch00, half_count, en_sum_lt_zero00 ); //(sum - count / 2) / count
					PDX_ADDS_MX32_T(out1, sum_ch01, half_count, en_sum_ge_zero01 );
					PDX_SUBS_MX32_T(out1, sum_ch01, half_count, en_sum_lt_zero01 );
					PDX_ADDS_MX32_T(out2, sum_ch10, half_count, en_sum_ge_zero10 );
					PDX_SUBS_MX32_T(out2, sum_ch10, half_count, en_sum_lt_zero10 );
					PDX_ADDS_MX32_T(out3, sum_ch11, half_count, en_sum_ge_zero11 );
					PDX_SUBS_MX32_T(out3, sum_ch11, half_count, en_sum_lt_zero11 );

					//dividing the sum by count value to get the average result
	                PDX_DIV_MX32(out0, rem, out0, count);
	                PDX_DIV_MX32(out1, rem, out1, count);
	                PDX_DIV_MX32(out2, rem, out2, count);
	                PDX_DIV_MX32(out3, rem, out3, count);
	                //saturating to 8 bits
					out0 = PDX_MIN_MX32 (PDX_MAX_MX32 (out0, act_min), act_max);
					out1 = PDX_MIN_MX32 (PDX_MAX_MX32 (out1, act_min), act_max);
					out2 = PDX_MIN_MX32 (PDX_MAX_MX32 (out2, act_min), act_max);
					out3 = PDX_MIN_MX32 (PDX_MAX_MX32 (out3, act_min), act_max);
					//save value at the end
					PDX_SAV32_MX8_XP (out0, pvOuta, pvOut, ch_src-ch_idx);
					PDX_SAV32_MX8_XP (out1, pvOuta, pvOut, ch_src-ch_idx - 8);
					PDX_SAV32_MX8_XP (out2, pvOuta, pvOut, ch_src-ch_idx - 8*2);
					PDX_SAV32_MX8_XP (out3, pvOuta, pvOut, ch_src-ch_idx - 8*3);

					PDX_SAPOS_MX8_FP(pvOuta,pvOut); //Flush out all the bytes

		        	src_in += 32; // update to next 32 channels

		        }
	        }
	    }

	    return ADI_SHARCFX_NN_SUCCESS;
	}



int32_t adi_sharcfx_avgpool_s8_get_buffer_size(const int output_x, const int ch_src)
{
    (void)output_x;

    (void)ch_src;
    return 0;

}




