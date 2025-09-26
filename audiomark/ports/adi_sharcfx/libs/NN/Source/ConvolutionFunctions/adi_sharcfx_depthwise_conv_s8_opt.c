
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
 * Title:        adi_sharcfx_depthwise_conv_s8_opt.c
 * Description:  Optimized s8 depthwise separable convolution function for
 *               channel multiplier of 1.
 * $Date:        22 April 2025
 * $Revision:     V.1.0
 *
 * Target Processor:  SHARC-FX Processor
 *
 * -------------------------------------------------------------------- */
#include <libs/NN/Include/adi_sharcfx_nnfunctions.h>
#include <libs/NN/Include/adi_sharcfx_nnsupportfunctions.h>

/**
 *  @ingroup Public
 */

/**
 * @addtogroup NNConv
 * @{
 */

/*
 * Optimized s8 depthwise convolution function with constraint that in_channel equals out_channel
 *
 *  Refer prototype header file for details.
 *
 */

adi_sharcfx_nn_status adi_sharcfx_depthwise_conv_s8_opt(const nn_context *ctx,
                                              const nn_dw_conv_params *dw_conv_params,
                                              const nn_per_channel_quant_params *quant_params,
                                              const nn_dims *input_dims,
                                              const int8_t *input,
                                              const nn_dims *filter_dims,
                                              const int8_t *kernel,
                                              const nn_dims *bias_dims,
                                              const int32_t *bias,
                                              const nn_dims *output_dims,
                                              int8_t *output)
{

    const int32_t input_ch = input_dims->c;
    const int32_t output_ch = output_dims->c;

    /* Check depth multiplier is 1 */
    if (input_ch != output_ch)
    {
        return ADI_SHARCFX_NN_ARG_ERROR;
    }

    if (ctx->buf == NULL && adi_sharcfx_depthwise_conv_s8_opt_get_buffer_size(input_dims, filter_dims) > 0)
    {
        return ADI_SHARCFX_NN_ARG_ERROR;
    }
    return adi_sharcfx_depthwise_conv_s8(ctx,
                                 dw_conv_params,
                                 quant_params,
                                 input_dims,
                                 input,
                                 filter_dims,
                                 kernel,
                                 bias_dims,
                                 bias,
                                 output_dims,
                                 output);

    //Return to application
    return ADI_SHARCFX_NN_SUCCESS;
}

int32_t adi_sharcfx_depthwise_conv_s8_opt_get_buffer_size(const nn_dims *input_dims, const nn_dims *filter_dims)
{
    (void)input_dims;
    (void)filter_dims;
    return 0;

}

/**
 * @} end of NNConv group
 */
