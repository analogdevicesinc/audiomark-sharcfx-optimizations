/******************************************************************************
 * @file     fast_math_functions.h
 * @brief    Public header file for DSP Library
 * $Date:        22 April 2025
 * $Revision:    V.1.0
 *
 * Target Processor: SHARC-FX Processor
 ******************************************************************************/
/*
 * Copyright(c) 2025 Analog Devices, Inc. All Rights Reserved.
 * Copyright (c) 2010-2020 Arm Limited or its affiliates. All rights reserved.
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

 
#ifndef _FAST_MATH_FUNCTIONS_H_
#define _FAST_MATH_FUNCTIONS_H_

#include "../../../DSP/Include/adi_sharcfx_math_types.h"
#include "../../../DSP/Include/dsp/basic_math_functions.h"


#ifdef   __cplusplus
extern "C"
{
#endif



/**
  @brief         Floating-point vector of log values.
  @param[in]     pSrc       points to the input vector
  @param[out]    pDst       points to the output vector
  @param[in]     blockSize  number of samples in each vector
  @return        none
 */
  void adi_sharcfx_vlog_f32(
  const float32_t * pSrc,
        float32_t * pDst,
        uint32_t blockSize);



/**
  @brief         Floating-point vector of exp values.
  @param[in]     pSrc       points to the input vector
  @param[out]    pDst       points to the output vector
  @param[in]     blockSize  number of samples in each vector
  @return        none
 */
  void adi_sharcfx_vexp_f32(
  const float32_t * pSrc,
        float32_t * pDst,
        uint32_t blockSize);


#define FAST_MATH_TABLE_SIZE  512

#ifdef   __cplusplus
}
#endif

#endif /* ifndef _FAST_MATH_FUNCTIONS_H_ */
