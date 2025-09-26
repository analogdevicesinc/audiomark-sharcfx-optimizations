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
 * Title:        adi_sharcfx_nnfunctions.h
 * Description:  Public header file for NN Library
 *
 * $Date:        22 April 2025
 * $Revision:    V.1.0
 *
 * Target Processor:  SHARC-FX Processor
 * -------------------------------------------------------------------- */


/**
 * @defgroup Public Public
 * A collection of functions to perform basic operations for neural network layers. Functions with a _s8 suffix support
 * TensorFlow Lite framework.
 */

#ifndef _ADI_SHARCFX_NNFUNCTIONS_H
#define _ADI_SHARCFX_NNFUNCTIONS_H

#include "../../NN/Include/adi_sharcfx_nn_math_types.h"
#include "../../NN/Include/adi_sharcfx_nn_types.h"

#define USE_INTRINSIC

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @defgroup NNConv Convolution Functions
 *
 * Collection of convolution, depthwise convolution functions and their variants.
 *
 * The convolution is implemented in 2 steps: im2col and General Matrix Multiplication(GEMM)
 *
 * im2col is a process of converting each patch of image data into
 * a column. After im2col, the convolution is computed as matrix-matrix
 * multiplication.
 *
 * To reduce the memory footprint, the im2col is performed partially.
 * Each iteration, only a few column (i.e., patches) are generated followed
 * by GEMM.
 *
 */

/**
 * @brief s8 convolution layer wrapper function with the main purpose to call the optimal kernel available in
 *        NN  to perform the convolution.
 *
 * @param[in, out] ctx            Function context that contains the additional buffer if required by the function.
 *                                adi_sharcfx_convolve_wrapper_s8_get_buffer_size will return the buffer_size if required.
 *                                The caller is expected to clear the buffer ,if applicable, for security reasons.
 * @param[in]      conv_params    Convolution parameters (e.g. strides, dilations, pads,...).
 *                                Range of conv_params->input_offset  : [-127, 128]
 *                                Range of conv_params->output_offset : [-128, 127]
 * @param[in]      quant_params   Per-channel quantization info.
 *                                It contains the multiplier and shift values to be applied to each output channel
 * @param[in]      input_dims     Input (activation) tensor dimensions. Format: [N, H, W, C_IN]
 * @param[in]      input_data     Input (activation) data pointer. Data type: int8
 * @param[in]      filter_dims    Filter tensor dimensions. Format: [C_OUT, HK, WK, C_IN] where HK and WK are the
 *                                spatial filter dimensions
 * @param[in]      filter_data    Filter data pointer. Data type: int8
 * @param[in]      bias_dims      Bias tensor dimensions. Format: [C_OUT]
 * @param[in]      bias_data      Bias data pointer. Data type: int32
 * @param[in]      output_dims    Output tensor dimensions. Format: [N, H, W, C_OUT]
 * @param[out]     output_data    Output data pointer. Data type: int8
 *
 * @return     The function returns either
 *                  <code>ADI_SHARCFX_NN_ARG_ERROR</code> if argument constraints fail. or,
 *                  <code>ADI_SHARCFX_NN_SUCCESS</code> on successful completion.
 *
 */
adi_sharcfx_nn_status adi_sharcfx_convolve_wrapper_s8(const nn_context *ctx,
                                            const nn_conv_params *conv_params,
                                            const nn_per_channel_quant_params *quant_params,
                                            const nn_dims *input_dims,
                                            const int8_t *input_data,
                                            const nn_dims *filter_dims,
                                            const int8_t *filter_data,
                                            const nn_dims *bias_dims,
                                            const int32_t *bias_data,
                                            const nn_dims *output_dims,
                                            int8_t *output_data);

/**
 * @brief Get the required buffer size for adi_sharcfx_convolve_wrapper_s8
 *
 * @param[in]      conv_params    Convolution parameters (e.g. strides, dilations, pads,...).
 *                                Range of conv_params->input_offset  : [-127, 128]
 *                                Range of conv_params->output_offset : [-128, 127]
 * @param[in]      input_dims     Input (activation) dimensions. Format: [N, H, W, C_IN]
 * @param[in]      filter_dims    Filter dimensions. Format: [C_OUT, HK, WK, C_IN] where HK and WK are the spatial
 *                                filter dimensions
 * @param[in]      output_dims    Output tensor dimensions. Format: [N, H, W, C_OUT]
 *
 * @return         The function returns  required buffer size(bytes)
 *
 */
int32_t adi_sharcfx_convolve_wrapper_s8_get_buffer_size(const nn_conv_params *conv_params,
                                                const nn_dims *input_dims,
                                                const nn_dims *filter_dims,
                                                const nn_dims *output_dims);

/**
 * @brief Basic s8 convolution function
 * @param[in, out] ctx            Function context that contains the additional buffer if required by the function.
 *                                adi_sharcfx_convolve_s8_get_buffer_size will return the buffer_size if required.
 *                                The caller is expected to clear the buffer ,if applicable, for security reasons.
 * @param[in]      conv_params    Convolution parameters (e.g. strides, dilations, pads,...).
 *                                Range of conv_params->input_offset  : [-127, 128]
 *                                Range of conv_params->output_offset : [-128, 127]
 * @param[in]      quant_params   Per-channel quantization info.
 *                                It contains the multiplier and shift values to be applied to each output channel
 * @param[in]      input_dims     Input (activation) tensor dimensions. Format: [N, H, W, C_IN]
 * @param[in]      input_data     Input (activation) data pointer. Data type: int8
 * @param[in]      filter_dims    Filter tensor dimensions. Format: [C_OUT, HK, WK, C_IN] where HK and WK are the
 *                                spatial filter dimensions
 * @param[in]      filter_data    Filter data pointer. Data type: int8
 * @param[in]      bias_dims      Bias tensor dimensions. Format: [C_OUT]
 * @param[in]      bias_data      Optional bias data pointer. Data type: int32
 * @param[in]      output_dims    Output tensor dimensions. Format: [N, H, W, C_OUT]
 * @param[out]     output_data    Output data pointer. Data type: int8

 * @return     The function returns <code>ADI_SHARCFX_NN_SUCCESS</code>
 *
 * @details
 *    1. Supported framework: TensorFlow Lite micro
 *    2. Additional memory is required for optimization. Refer to argument 'ctx' for details.
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
                                    int8_t *output_data);

/**
 * @brief Get the required buffer size for s8 convolution function
 *
 * @param[in]       input_dims            Input (activation) tensor dimensions. Format: [N, H, W, C_IN]
 * @param[in]       filter_dims           Filter tensor dimensions. Format: [C_OUT, HK, WK, C_IN] where HK and WK
 * are the spatial filter dimensions
 * @return          The function returns  required buffer size(bytes)
 *
 */
int32_t adi_sharcfx_convolve_s8_get_buffer_size(const nn_dims *input_dims, const nn_dims *filter_dims);



/**
 * @brief Fast s8 version for 1x1 convolution (non-square shape)
 *
 * @param[in, out] ctx           Function context that contains the additional buffer if required by the function.
 *                               adi_sharcfx_convolve_1x1_s8_fast_get_buffer_size will return the buffer_size if required.
 *                               The caller is expected to clear the buffer ,if applicable, for security reasons.
 * @param[in]      conv_params   Convolution parameters (e.g. strides, dilations, pads,...).
 *                               Range of conv_params->input_offset  : [-127, 128]
 *                               Range of conv_params->output_offset : [-128, 127]
 * @param[in]      quant_params  Per-channel quantization info.
 *                               It contains the multiplier and shift values to be applied to each output channel
 * @param[in]      input_dims    Input (activation) tensor dimensions. Format: [N, H, W, C_IN]
 * @param[in]      input_data    Input (activation) data pointer. Data type: int8
 * @param[in]      filter_dims   Filter tensor dimensions. Format: [C_OUT, 1, 1, C_IN]
 * @param[in]      filter_data   Filter data pointer. Data type: int8
 * @param[in]      bias_dims     Bias tensor dimensions. Format: [C_OUT]
 * @param[in]      bias_data     Optional bias data pointer. Data type: int32
 * @param[in]      output_dims   Output tensor dimensions. Format: [N, H, W, C_OUT]
 * @param[out]     output_data   Output data pointer. Data type: int8
 *
 * @return     The function returns either
 *                  <code>ADI_SHARCFX_NN_ARG_ERROR</code> if argument constraints fail. or,
 *                  <code>ADI_SHARCFX_NN_SUCCESS</code> on successful completion.
 *
 * @details
 *   - Supported framework : TensorFlow Lite Micro
 *   - The following constrains on the arguments apply
 *      -# input_dims->c is a multiple of 4
 *      -# conv_params->padding.w = conv_params->padding.h = 0
 *      -# conv_params->stride.w = conv_params->stride.h = 1
 *
 */
adi_sharcfx_nn_status adi_sharcfx_convolve_1x1_s8_fast(const nn_context *ctx,
                                             const nn_conv_params *conv_params,
                                             const nn_per_channel_quant_params *quant_params,
                                             const nn_dims *input_dims,
                                             const int8_t *input_data,
                                             const nn_dims *filter_dims,
                                             const int8_t *filter_data,
                                             const nn_dims *bias_dims,
                                             const int32_t *bias_data,
                                             const nn_dims *output_dims,
                                             int8_t *output_data);

/**
 * @brief Get the required buffer size for adi_sharcfx_convolve_1x1_s8_fast
 *
 * @param[in]       input_dims            Input (activation) dimensions
 * @return          The function returns the required buffer size in bytes
 *
 */
int32_t adi_sharcfx_convolve_1x1_s8_fast_get_buffer_size(const nn_dims *input_dims);

/**
 * @brief s8 version for 1x1 convolution with support for non-unity stride values
 *
 * @param[in, out] ctx           Function context that contains the additional buffer if required by the function.
 *                               None is required by this function.
 * @param[in]      conv_params   Convolution parameters (e.g. strides, dilations, pads,...).
 *                               Range of conv_params->input_offset  : [-127, 128]
 *                               Range of conv_params->output_offset : [-128, 127]
 * @param[in]      quant_params  Per-channel quantization info.
 *                               It contains the multiplier and shift values to be applied to each output channel
 * @param[in]      input_dims    Input (activation) tensor dimensions. Format: [N, H, W, C_IN]
 * @param[in]      input_data    Input (activation) data pointer. Data type: int8
 * @param[in]      filter_dims   Filter tensor dimensions. Format: [C_OUT, 1, 1, C_IN]
 * @param[in]      filter_data   Filter data pointer. Data type: int8
 * @param[in]      bias_dims     Bias tensor dimensions. Format: [C_OUT]
 * @param[in]      bias_data     Optional bias data pointer. Data type: int32
 * @param[in]      output_dims   Output tensor dimensions. Format: [N, H, W, C_OUT]
 * @param[out]     output_data   Output data pointer. Data type: int8
 *
 * @return     The function returns either
 *                  <code>ADI_SHARCFX_NN_ARG_ERROR</code> if argument constraints fail. or,
 *                  <code>ADI_SHARCFX_NN_SUCCESS</code> on successful completion.
 * @details
 *   - Supported framework : TensorFlow Lite Micro
 *   - The following constrains on the arguments apply
 *      -# conv_params->padding.w = conv_params->padding.h = 0
 *
 */
adi_sharcfx_nn_status adi_sharcfx_convolve_1x1_s8(const nn_context *ctx,
                                        const nn_conv_params *conv_params,
                                        const nn_per_channel_quant_params *quant_params,
                                        const nn_dims *input_dims,
                                        const int8_t *input_data,
                                        const nn_dims *filter_dims,
                                        const int8_t *filter_data,
                                        const nn_dims *bias_dims,
                                        const int32_t *bias_data,
                                        const nn_dims *output_dims,
                                        int8_t *output_data);

/**
 * @brief 1xn convolution
 *
 * @param[in, out] ctx           Function context that contains the additional buffer if required by the function.
 *                               adi_sharcfx_convolve_1_x_n_s8_get_buffer_size will return the buffer_size if required
 *                               The caller is expected to clear the buffer ,if applicable, for security reasons.
 * @param[in]      conv_params   Convolution parameters (e.g. strides, dilations, pads,...).
 *                               Range of conv_params->input_offset  : [-127, 128]
 *                               Range of conv_params->output_offset : [-128, 127]
 * @param[in]      quant_params  Per-channel quantization info.
 *                               It contains the multiplier and shift values to be applied to each output channel
 * @param[in]      input_dims    Input (activation) tensor dimensions. Format: [N, H, W, C_IN]
 * @param[in]      input_data    Input (activation) data pointer. Data type: int8
 * @param[in]      filter_dims   Filter tensor dimensions. Format: [C_OUT, 1, WK, C_IN] where WK is the horizontal
 *                               spatial filter dimension
 * @param[in]      filter_data   Filter data pointer. Data type: int8
 * @param[in]      bias_dims     Bias tensor dimensions. Format: [C_OUT]
 * @param[in]      bias_data     Optional bias data pointer. Data type: int32
 * @param[in]      output_dims   Output tensor dimensions. Format: [N, H, W, C_OUT]
 * @param[out]     output_data   Output data pointer. Data type: int8
 *
 * @return     The function returns either
 *                  <code>ADI_SHARCFX_NN_ARG_ERROR</code> if argument constraints fail. or,
 *                  <code>ADI_SHARCFX_NN_SUCCESS</code> on successful completion.
 *
 * @details
 *   - Supported framework : TensorFlow Lite Micro
 *   - The following constrains on the arguments apply
 *      -# input_dims->n equals 1
 *      -# ouput_dims->w is a multiple of 4
 *      -# Explicit constraints(since it is for 1xN convolution)
 *      -## input_dims->h equals 1
 *      -## output_dims->h equals 1
 *      -## filter_dims->h equals 1
 *@todo  Remove constraint on output_dims->w to make the function generic.
 *
 */
adi_sharcfx_nn_status adi_sharcfx_convolve_1_x_n_s8(const nn_context *ctx,
                                          const nn_conv_params *conv_params,
                                          const nn_per_channel_quant_params *quant_params,
                                          const nn_dims *input_dims,
                                          const int8_t *input_data,
                                          const nn_dims *filter_dims,
                                          const int8_t *filter_data,
                                          const nn_dims *bias_dims,
                                          const int32_t *bias_data,
                                          const nn_dims *output_dims,
                                          int8_t *output_data);

/**
 * @brief Get the required additional buffer size for 1xn convolution
 *
 * @param[in]       input_dims            Input (activation) tensor dimensions. Format: [N, H, W, C_IN]
 * @param[in]       filter_dims           Filter tensor dimensions. Format: [C_OUT, 1, WK, C_IN] where WK is the
 *                                        horizontal spatial filter dimension
 * @return          The function returns  required buffer size(bytes)
 *
 */
int32_t adi_sharcfx_convolve_1_x_n_s8_get_buffer_size(const nn_dims *input_dims, const nn_dims *filter_dims);

/**
 * @brief Wrapper function to pick the right optimized s8 depthwise convolution function
 *
 * @param[in, out] ctx             Function context (e.g. temporary buffer). Check the function
 *                                 definition file to see if an additional buffer is required.
 *                                 Optional function {API}_get_buffer_size() provides the buffer
 *                                 size if required.
 *                                 The caller is expected to clear the buffer ,if applicable, for security reasons.
 * @param[in]      dw_conv_params  Depthwise convolution parameters (e.g. strides, dilations, pads,...)
 *                                 dw_conv_params->dilation is not used.
 *                                 Range of dw_conv_params->input_offset : [-127, 128]
 *                                 Range of dw_conv_params->output_offset : [-128, 127]
 * @param[in]      quant_params    Per-channel quantization info.
 *                                 It contains the multiplier and shift values to be applied to each
 *                                 output channel
 * @param[in]      input_dims      Input (activation) tensor dimensions. Format: [H, W, C_IN]
 *                                 Batch argument N is not used and assumed to be 1.
 * @param[in]      input_data      Input (activation) data pointer. Data type: int8
 * @param[in]      filter_dims     Filter tensor dimensions. Format: [1, H, W, C_OUT]
 * @param[in]      filter_data     Filter data pointer. Data type: int8
 * @param[in]      bias_dims       Bias tensor dimensions. Format: [C_OUT]
 * @param[in]      bias_data       Bias data pointer. Data type: int32
 * @param[in]      output_dims     Output tensor dimensions. Format: [1, H, W, C_OUT]
 * @param[in, out] output_data     Output data pointer. Data type: int8
 * @return     The function returns
 *                <code>ADI_SHARCFX_NN_SUCCESS</code>   -  Successful completion.
 *
 * @details
 *    - Supported framework: TensorFlow Lite
 *    - Picks one of the the following functions
 *        -# adi_sharcfx_depthwise_conv_s8()
 *        -# adi_sharcfx_depthwise_conv_3x3_s8() - Cortex-M CPUs with DSP extension only
 *        -# adi_sharcfx_depthwise_conv_s8_opt()
 *    - Check details of adi_sharcfx_depthwise_conv_s8_opt() for potential data that can be accessed outside of the
 * boundary.
 */
adi_sharcfx_nn_status adi_sharcfx_depthwise_conv_wrapper_s8(const nn_context *ctx,
                                                  const nn_dw_conv_params *dw_conv_params,
                                                  const nn_per_channel_quant_params *quant_params,
                                                  const nn_dims *input_dims,
                                                  const int8_t *input_data,
                                                  const nn_dims *filter_dims,
                                                  const int8_t *filter_data,
                                                  const nn_dims *bias_dims,
                                                  const int32_t *bias_data,
                                                  const nn_dims *output_dims,
                                                  int8_t *output_data);

/**
 * @brief Get size of additional buffer required by adi_sharcfx_depthwise_conv_wrapper_s8()
 *
 * @param[in]      dw_conv_params  Depthwise convolution parameters (e.g. strides, dilations, pads,...)
 *                                 Range of dw_conv_params->input_offset : [-127, 128]
 *                                 Range of dw_conv_params->input_offset : [-128, 127]
 * @param[in]      input_dims      Input (activation) tensor dimensions. Format: [H, W, C_IN]
 *                                 Batch argument N is not used and assumed to be 1.
 * @param[in]      filter_dims     Filter tensor dimensions. Format: [1, H, W, C_OUT]
 * @param[in]      output_dims     Output tensor dimensions. Format: [1, H, W, C_OUT]
 * @return                         Size of additional memory required for optimizations in bytes.
 *
 */
int32_t adi_sharcfx_depthwise_conv_wrapper_s8_get_buffer_size(const nn_dw_conv_params *dw_conv_params,
                                                      const nn_dims *input_dims,
                                                      const nn_dims *filter_dims,
                                                      const nn_dims *output_dims);

/**
 * @brief Basic s8 depthwise convolution function that doesn't have any constraints on the input dimensions.
 *
 * @param[in, out] ctx             Function context (e.g. temporary buffer). Check the function
 *                                 definition file to see if an additional buffer is required.
 *                                 Optional function {API}_get_buffer_size() provides the buffer
 *                                 size if an additional buffer is required exists if additional memory is.
 *                                 The caller is expected to clear the buffer ,if applicable, for security reasons.
 * @param[in]      dw_conv_params  Depthwise convolution parameters (e.g. strides, dilations, pads,...)
 *                                 dw_conv_params->dilation is not used.
 *                                 Range of dw_conv_params->input_offset : [-127, 128]
 *                                 Range of dw_conv_params->input_offset : [-128, 127]
 * @param[in]      quant_params    Per-channel quantization info.
 *                                 It contains the multiplier and shift values to be applied to each
 *                                 output channel
 * @param[in]      input_dims      Input (activation) tensor dimensions. Format: [N, H, W, C_IN]
 *                                 Batch argument N is not used.
 * @param[in]      input_data      Input (activation) data pointer. Data type: int8
 * @param[in]      filter_dims     Filter tensor dimensions. Format: [1, H, W, C_OUT]
 * @param[in]      filter_data     Filter data pointer. Data type: int8
 * @param[in]      bias_dims       Bias tensor dimensions. Format: [C_OUT]
 * @param[in]      bias_data       Bias data pointer. Data type: int32
 * @param[in]      output_dims     Output tensor dimensions. Format: [N, H, W, C_OUT]
 * @param[in, out] output_data     Output data pointer. Data type: int8
 * @return     The function returns <code>ADI_SHARCFX_NN_SUCCESS</code>
 *
 * @details
 *    - Supported framework: TensorFlow Lite
 */
adi_sharcfx_nn_status adi_sharcfx_depthwise_conv_s8(const nn_context *ctx,
                                          const nn_dw_conv_params *dw_conv_params,
                                          const nn_per_channel_quant_params *quant_params,
                                          const nn_dims *input_dims,
                                          const int8_t *input_data,
                                          const nn_dims *filter_dims,
                                          const int8_t *filter_data,
                                          const nn_dims *bias_dims,
                                          const int32_t *bias_data,
                                          const nn_dims *output_dims,
                                          int8_t *output_data);




adi_sharcfx_nn_status adi_sharcfx_depthwise_conv_3x3_s8(const nn_context *ctx,
                                              const nn_dw_conv_params *dw_conv_params,
                                              const nn_per_channel_quant_params *quant_params,
                                              const nn_dims *input_dims,
                                              const int8_t *input_data,
                                              const nn_dims *filter_dims,
                                              const int8_t *filter_data,
                                              const nn_dims *bias_dims,
                                              const int32_t *bias_data,
                                              const nn_dims *output_dims,
                                              int8_t *output_data);

/**
 * @brief Optimized s8 depthwise convolution function with constraint that in_channel equals out_channel.
 *        Refer adi_sharcfx_depthwise_conv_s8() for function argument details.
 *
 * @return     The function returns one of the following
 *                <code>ADI_SHARCFX_NN_ARG_ERROR</code> - input channel != output channel or
 *                                                      ch_mult != 1
 *                <code>ADI_SHARCFX_NN_SUCCESS</code> - Successful operation
 *
 * @note       If number of channels is not a multiple of 4, upto 3 elements outside the boundary will be read out
 *             for the following if MVE optimizations(Arm Helium Technology) are used.
 *               - Output shift
 *               - Output multiplier
 *               - Output bias
 *               - kernel
 * @details
 *    - Supported framework: TensorFlow Lite
 *    - The following constrains on the arguments apply
 *        -# Number of input channel equals number of output channels or ch_mult equals 1
 *    - Reccomended when number of channels is 4 or greater.
 *
 */
adi_sharcfx_nn_status adi_sharcfx_depthwise_conv_s8_opt(const nn_context *ctx,
                                              const nn_dw_conv_params *dw_conv_params,
                                              const nn_per_channel_quant_params *quant_params,
                                              const nn_dims *input_dims,
                                              const int8_t *input_data,
                                              const nn_dims *filter_dims,
                                              const int8_t *filter_data,
                                              const nn_dims *bias_dims,
                                              const int32_t *bias_data,
                                              const nn_dims *output_dims,
                                              int8_t *output_data);

/**
 * @brief Get the required buffer size for optimized s8 depthwise convolution
 * function with constraint that in_channel equals out_channel.
 * @param[in]       input_dims   Input (activation) tensor dimensions. Format: [1, H, W, C_IN]
 *                               Batch argument N is not used.
 * @param[in]       filter_dims  Filter tensor dimensions. Format: [1, H, W, C_OUT]
 * @return          The function returns  required buffer size in bytes
 *
 */
int32_t adi_sharcfx_depthwise_conv_s8_opt_get_buffer_size(const nn_dims *input_dims, const nn_dims *filter_dims);

/**
 * @defgroup FC Fully-connected Layer Functions
 *
 * Collection of fully-connected and matrix multiplication functions.
 *
 * Fully-connected layer is basically a matrix-vector multiplication
 * with bias. The matrix is the weights and the input/output vectors
 * are the activation values. Supported {weight, activation} precisions
 * include {8-bit, 8-bit} and {8-bit, 16-bit}
 *
 *
 */

/**
 * @brief Basic s8 Fully Connected function.
 *
 * @param[in, out] ctx           Function context (e.g. temporary buffer). Check the function
 *                               definition file to see if an additional buffer is required.
 *                               Optional function {API}_get_buffer_size() provides the buffer
 *                               size if an additional buffer is required.
 *                               The caller is expected to clear the buffer ,if applicable, for security reasons.
 * @param[in]      fc_params     Fully Connected layer parameters.
 *                               Range of fc_params->input_offset  : [-127, 128]
 *                               fc_params->filter_offset : 0
 *                               Range of fc_params->output_offset : [-128, 127]
 * @param[in]      quant_params  Per-tensor quantization info.
 *                               It contains the multiplier and shift values to be applied to the output tensor.
 * @param[in]      input_dims    Input (activation) tensor dimensions. Format: [N, H, W, C_IN]
 *                               Input dimension is taken as Nx(H * W * C_IN)
 * @param[in]      input_data    Input (activation) data pointer. Data type: int8
 * @param[in]      filter_dims   Two dimensional filter dimensions. Format: [N, C]
 *                               N : accumulation depth and equals (H * W * C_IN) from input_dims
 *                               C : output depth and equals C_OUT in output_dims
 *                               H & W : Not used
 * @param[in]      filter_data   Filter data pointer. Data type: int8
 * @param[in]      bias_dims     Bias tensor dimensions. Format: [C_OUT]
 *                               N, H, W : Not used
 * @param[in]      bias_data     Bias data pointer. Data type: int32
 * @param[in]      output_dims   Output tensor dimensions. Format: [N, C_OUT]
 *                               N : Batches
 *                               C_OUT : Output depth
 *                               H & W : Not used.
 * @param[in, out] output_data    Output data pointer. Data type: int8
 * @return     The function returns <code>ADI_SHARCFX_NN_SUCCESS</code>
 *
 * @details
 *    - Supported framework: TensorFlow Lite
 */
adi_sharcfx_nn_status adi_sharcfx_fully_connected_s8(const nn_context *ctx,
                                           const nn_fc_params *fc_params,
                                           const nn_per_tensor_quant_params *quant_params,
                                           const nn_dims *input_dims,
                                           const int8_t *input_data,
                                           const nn_dims *filter_dims,
                                           const int8_t *filter_data,
                                           const nn_dims *bias_dims,
                                           const int32_t *bias_data,
                                           const nn_dims *output_dims,
                                           int8_t *output_data);

/**
 * @brief Get the required buffer size for S8 basic fully-connected and
 * matrix multiplication layer function for TF Lite
 * @param[in]      filter_dims             dimension of filter
 * @return         The function returns    required buffer size in bytes
 *
 */
int32_t adi_sharcfx_fully_connected_s8_get_buffer_size(const nn_dims *filter_dims);


/**
 * @brief s8 average pooling function.
 *
 * @param[in, out] ctx          Function context (e.g. temporary buffer). Check the function
 *                              definition file to see if an additional buffer is required.
 *                              Optional function {API}_get_buffer_size() provides the buffer
 *                              size if an additional buffer is required.
 *                              The caller is expected to clear the buffer ,if applicable, for security reasons.
 * @param[in]      pool_params  Pooling parameters
 * @param[in]      input_dims   Input (activation) tensor dimensions. Format: [H, W, C_IN]
 *                              Argument 'N' is not used.
 * @param[in]      input_data   Input (activation) data pointer. Data type: int8
 * @param[in]      filter_dims  Filter tensor dimensions. Format: [H, W]
 *                              Argument N and C are not used.
 * @param[in]      output_dims  Output tensor dimensions. Format: [H, W, C_OUT]
 *                              Argument N is not used.
 *                              C_OUT equals C_IN.
 * @param[in, out] output_data Output data pointer. Data type: int8
 * @return                     The function returns
 *                             <code>ADI_SHARCFX_NN_SUCCESS</code> - Successful operation
 *
 * @details
 *    - Supported Framework: TensorFlow Lite
 *
 */
adi_sharcfx_nn_status adi_sharcfx_avgpool_s8(const nn_context *ctx,
                                   const nn_pool_params *pool_params,
                                   const nn_dims *input_dims,
                                   const int8_t *input_data,
                                   const nn_dims *filter_dims,
                                   const nn_dims *output_dims,
                                   int8_t *output_data);

/**
 * @brief Get the required buffer size for S8 average pooling function
 * @param[in]       dim_dst_width         output tensor dimension
 * @param[in]       ch_src                number of input tensor channels
 * @return          The function returns  required buffer size in bytes
 *
 */
int32_t adi_sharcfx_avgpool_s8_get_buffer_size(const int dim_dst_width, const int ch_src);


/**
 * @defgroup Softmax Softmax Functions
 *
 *
 */

/**
 * @brief S8 softmax function
 * @param[in]  input     Pointer to the input tensor
 * @param[in]  num_rows  Number of rows in the input tensor
 * @param[in]  row_size  Number of elements in each input row
 * @param[in]  mult      Input quantization multiplier
 * @param[in]  shift     Input quantization shift within the range [0, 31]
 * @param[in]  diff_min  Minimum difference with max in row. Used to check if
 *                       the quantized exponential operation can be performed
 * @param[out] output    Pointer to the output tensor
 *
 * @note Supported framework: TensorFlow Lite micro (bit-accurate)
 *
 */
void adi_sharcfx_softmax_s8(const int8_t *input,
                    const int32_t num_rows,
                    const int32_t row_size,
                    const int32_t mult,
                    const int32_t shift,
                    const int32_t diff_min,
                    int8_t *output);




#ifdef __cplusplus
}
#endif

#endif
