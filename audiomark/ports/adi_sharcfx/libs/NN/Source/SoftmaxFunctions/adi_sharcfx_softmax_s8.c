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
 * Title:        adi_sharcfx_softmax_s8.c
 * Description:  S8 softmax function
 * $Date:        22 April 2025
 * $Revision:    V.1.0
 * 
 * Target Processor: SHARC-FX Processor
 * -------------------------------------------------------------------- */


#include <libs/NN/Include/adi_sharcfx_nnfunctions.h>
#include <libs/NN/Include/adi_sharcfx_nnsupportfunctions.h>


#define ACCUM_BITS 12

void adi_sharcfx_softmax_s8(const int8_t *input,
                    const int32_t num_rows,
                    const int32_t row_size,
                    const int32_t mult,
                    const int32_t shift,
                    const int32_t diff_min,
                    int8_t *output)
{

    const int32_t mask = (1 << shift);

    int32_t col = 0;
    int32_t row_idx;

    for (row_idx = 0; row_idx < num_rows; ++row_idx)
    {

        xb_int8 max_reduced;
        xb_vec4Mx8 max_vector, first_col, second_col;

        xb_vec4Mx8* pinput;
        pinput = (xb_vec4Mx8 *)(input);
        valign pinputa = PDX_LA_4MX8_PP(pinput);
        PDX_LA_4MX8_XP(first_col,pinputa, pinput, 32);

        if(row_size>32){
			for (col = 32; col < row_size; col+=32)
			{
				PDX_LA_4MX8_XP(second_col, pinputa, pinput, 32);
				max_vector = PDX_MAX_4MX8(first_col, second_col);
				first_col = second_col;
			}
        }
        else{
        	max_vector = first_col;
        }
        max_reduced = PDX_RMAX_4MX8(max_vector);

        int8_t max = max_reduced;
        int32_t diff = 0;
        int32_t sum = 0;int32_t multiplyval;

        for (col = 0; col < row_size; ++col)
        {
            diff = input[col] - max;
            if (diff >= diff_min)
            {
            	multiplyval = MUL_SAT(diff *mask, mult);
                sum += DIV_POW2(EXP_ON_NEG(multiplyval), ACCUM_BITS);

            }
        }


        const int32_t headroom = __CLZ(sum);
        const int32_t shifted_scale = ONE_OVER1((sum > 0 ? sum << headroom : 0) - (1 << 31));
        int32_t bits_over_unit;

            int8_t *output_s8 = (int8_t *)output + row_idx * row_size;

            bits_over_unit = ACCUM_BITS - headroom + 23;
            int in, out;
            for (col = 0; col < row_size; ++col)
            {
                diff = input[col] - max;
                 in = input[col];
                if (diff >= diff_min)
                {
                    const int32_t res =
                        DIV_POW2(MUL_SAT(shifted_scale, EXP_ON_NEG(MUL_SAT(diff * mask, mult))), bits_over_unit) +
                        NN_Q7_MIN;
                    output_s8[col] = (int8_t)CLAMP(res, (int32_t)NN_Q7_MAX, (int32_t)NN_Q7_MIN);
                }
                else
                {
                    output_s8[col] = NN_Q7_MIN;
                }
                 out =  output_s8[col] ;
            }

        input += row_size;
    }

}
