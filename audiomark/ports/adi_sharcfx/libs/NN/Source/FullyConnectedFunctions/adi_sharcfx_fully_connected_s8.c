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
 * Title:        adi_sharcfx_fully_connected_s8
 * Description:  Fully connected function compatible with TF Lite.
 *
 * $Date:        22 April 2025
 * $Revision:    V.1.0
 * 
 * Target Processor: SHARC-FX Processor
 * -------------------------------------------------------------------- */



#include <libs/NN/Include/adi_sharcfx_nnfunctions.h>
#include <libs/NN/Include/adi_sharcfx_nnsupportfunctions.h>


/*
 * S8 basic fully-connected and matrix multiplication layer function for TensorFlow Lite
 *
 * Refer header file for details.
 *
 */

adi_sharcfx_nn_status adi_sharcfx_fully_connected_s8(const nn_context *ctx,
                                           const nn_fc_params *fc_params,
                                           const nn_per_tensor_quant_params *quant_params,
                                           const nn_dims *input_dims,
                                           const int8_t *input,
                                           const nn_dims *filter_dims,
                                           const int8_t *kernel,
                                           const nn_dims *bias_dims,
                                           const int32_t *bias,
                                           const nn_dims *output_dims,
                                           int8_t *output)
{
    (void)bias_dims;
    (void)ctx;
    (void)fc_params->filter_offset;

    int32_t batch_cnt = input_dims->n;

    while (batch_cnt)
    {
        adi_sharcfx_nn_vec_mat_mult_t_s8(input,
                                 kernel,
                                 bias,
                                 output,
                                 fc_params->input_offset,
                                 fc_params->output_offset,
                                 quant_params->multiplier,
                                 quant_params->shift,
                                 filter_dims->n, /* col_dim or accum_depth */
                                 output_dims->c, /* row_dim or output_depth */
                                 fc_params->activation.min,
                                 fc_params->activation.max,
                                 1L);
        input += filter_dims->n;
        output += output_dims->c;
        batch_cnt--;
    }
    return (ADI_SHARCFX_NN_SUCCESS);
}

int32_t adi_sharcfx_fully_connected_s8_get_buffer_size(const nn_dims *filter_dims)
{
    (void)filter_dims;
    return 0;
}

