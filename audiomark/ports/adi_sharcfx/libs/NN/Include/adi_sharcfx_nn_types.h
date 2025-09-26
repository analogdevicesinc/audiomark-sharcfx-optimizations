 /*
 * Copyright(c) 2025 Analog Devices, Inc. All Rights Reserved

 * SPDX-FileCopyrightText: Copyright 2020-2022 Arm Limited and/or its affiliates <open-source-office@arm.com>
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
 * Title:        adi_sharcfx_nn_types.h
 * Description:  Public header file to contain the NN structs for the
 *               TensorFlowLite micro compliant functions
 *
 * $Date:        22 April 2025
 * $Revision:    V.1.0
 *
 * Target Processor:  SHARC-FX Processor
 * -------------------------------------------------------------------- */

#ifndef _ADI_SHARCFX_NN_TYPES_H
#define _ADI_SHARCFX_NN_TYPES_H

#include <stdint.h>

/** Enum for specifying activation function types */
typedef enum
{
    ADI_SHARCFX_SIGMOID = 0, /**< Sigmoid activation function */
    ADI_SHARCFX_TANH = 1,    /**< Tanh activation function */
} nn_activation_type;

/** NN object to contain the width and height of a tile */
typedef struct
{
    int32_t w; /**< Width */
    int32_t h; /**< Height */
} nn_tile;

/** NN object used for the function context. */
typedef struct
{
    void *buf;    /**< Pointer to a buffer needed for the optimization */
    int32_t size; /**< Buffer size */
} nn_context;

/** NN object to contain the dimensions of the tensors */
typedef struct
{
    int32_t n; /**< Generic dimension to contain either the batch size or output channels.
                     Please refer to the function documentation for more information */
    int32_t h; /**< Height */
    int32_t w; /**< Width */
    int32_t c; /**< Input channels */
} nn_dims;


/** NN object for the per-channel quantization parameters */
typedef struct
{
    int32_t *multiplier; /**< Multiplier values */
    int32_t *shift;      /**< Shift values */
} nn_per_channel_quant_params;

/** NN object for the per-tensor quantization parameters */
typedef struct
{
    int32_t multiplier; /**< Multiplier value */
    int32_t shift;      /**< Shift value */
} nn_per_tensor_quant_params;

/** NN object for the quantized Relu activation */
typedef struct
{
    int32_t min; /**< Min value used to clamp the result */
    int32_t max; /**< Max value used to clamp the result */
} nn_activation;

/** NN object for the convolution layer parameters */
typedef struct
{
    int32_t input_offset;  /**< Zero value for the input tensor */
    int32_t output_offset; /**< Zero value for the output tensor */
    nn_tile stride;
    nn_tile padding;
    nn_tile dilation;
    nn_activation activation;
} nn_conv_params;

/** NN object for Depthwise convolution layer parameters */
typedef struct
{
    int32_t input_offset;  /**< Zero value for the input tensor */
    int32_t output_offset; /**< Zero value for the output tensor */
    int32_t ch_mult;       /**< Channel Multiplier. ch_mult * in_ch = out_ch */
    nn_tile stride;
    nn_tile padding;
    nn_tile dilation;
    nn_activation activation;
} nn_dw_conv_params;
/** NN object for pooling layer parameters */
typedef struct
{
    nn_tile stride;
    nn_tile padding;
    nn_activation activation;
} nn_pool_params;

/** NN object for Fully Connected layer parameters */
typedef struct
{
    int32_t input_offset;  /**< Zero value for the input tensor */
    int32_t filter_offset; /**< Zero value for the filter tensor. Not used */
    int32_t output_offset; /**< Zero value for the output tensor */
    nn_activation activation;
} nn_fc_params;


/** NN object for quantization parameters */
typedef struct
{
    int32_t multiplier; /**< Multiplier value */
    int32_t shift;      /**< Shift value */
} nn_scaling;

/** NN norm layer coefficients */
typedef struct
{
    int16_t *input_weight;
    int16_t *forget_weight;
    int16_t *cell_weight;
    int16_t *output_weight;
} nn_layer_norm;



#endif // _ADI_SHARCFX_NN_TYPES_H
