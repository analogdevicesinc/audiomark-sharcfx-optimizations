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
 * Title:        adi_sharcfx_depthwise_conv_wrapper_s8.c
 * Description:  Wrapper API to select appropriate depthwise conv API based
 *               on dimensions.
 * $Date:        22 April 2025
 * $Revision:     V.1.0
 *
 * Target Processor:  SHARC-FX Processor
 * -------------------------------------------------------------------- */

#include <libs/NN/Include/adi_sharcfx_nnfunctions.h>


/*
 *  s8 Depthwise conv wrapper function
 *
 *  Refer header file for details.
 *
 */
adi_sharcfx_nn_status adi_sharcfx_depthwise_conv_wrapper_s8(const nn_context *ctx,
                                                  const nn_dw_conv_params *dw_conv_params,
                                                  const nn_per_channel_quant_params *quant_params,
                                                  const nn_dims *input_dims,
                                                  const int8_t *input,
                                                  const nn_dims *filter_dims,
                                                  const int8_t *filter,
                                                  const nn_dims *bias_dims,
                                                  const int32_t *bias,
                                                  const nn_dims *output_dims,
                                                  int8_t *output)
{
    adi_sharcfx_nn_status status = ADI_SHARCFX_NN_SUCCESS;
    if (1 == dw_conv_params->ch_mult && input_dims->n == 1 && dw_conv_params->dilation.w == 1 &&
        dw_conv_params->dilation.h == 1)
    {

        if ((filter_dims->w == 3) && (filter_dims->h == 3) && (dw_conv_params->padding.h <= 1) &&
            (dw_conv_params->padding.w <= 1))
        {
            status = adi_sharcfx_depthwise_conv_3x3_s8(ctx,
                                               dw_conv_params,
                                               quant_params,
                                               input_dims,
                                               input,
                                               filter_dims,
                                               filter,
                                               bias_dims,
                                               bias,
                                               output_dims,
                                               output);
        }

    }
    else
    {
        status = adi_sharcfx_depthwise_conv_s8(ctx,
                                       dw_conv_params,
                                       quant_params,
                                       input_dims,
                                       input,
                                       filter_dims,
                                       filter,
                                       bias_dims,
                                       bias,
                                       output_dims,
                                       output);
    }

    /* Return to application */
    return status;
}

int32_t adi_sharcfx_depthwise_conv_wrapper_s8_get_buffer_size(const nn_dw_conv_params *dw_conv_params,
                                                      const nn_dims *input_dims,
                                                      const nn_dims *filter_dims,
                                                      const nn_dims *output_dims)
{
    (void)dw_conv_params;
    int32_t size = 0;

    if (input_dims->c == output_dims->c && input_dims->n == 1 && dw_conv_params->dilation.w == 1 &&
        dw_conv_params->dilation.h == 1)
    {
        size = adi_sharcfx_depthwise_conv_s8_opt_get_buffer_size(input_dims, filter_dims);
    }

    return size;
}

