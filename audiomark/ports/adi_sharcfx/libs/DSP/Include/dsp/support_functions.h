/******************************************************************************
 * @file     support_functions.h
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


 
#ifndef _SUPPORT_FUNCTIONS_H_
#define _SUPPORT_FUNCTIONS_H_

#include "../../../DSP/Include/adi_sharcfx_math_types.h"

#ifdef   __cplusplus
extern "C"
{
#endif

  /**
   * @brief Converts the elements of the floating-point vector to Q15 vector.
   * @param[in]  pSrc       points to the floating-point input vector
   * @param[out] pDst       points to the Q15 output vector
   * @param[in]  blockSize  length of the input vector
   */
  void adi_sharcfx_float_to_q15(
  const float32_t * pSrc,
        q15_t * pDst,
        uint32_t blockSize);


  /**
   * @brief  Converts the elements of the Q15 vector to floating-point vector.
   * @param[in]  pSrc       is input pointer
   * @param[out] pDst       is output pointer
   * @param[in]  blockSize  is the number of samples to process
   */
  void adi_sharcfx_q15_to_float(
  const q15_t * pSrc,
        float32_t * pDst,
        uint32_t blockSize);

 
  /**
   * @brief  Copies the elements of a floating-point vector.
   * @param[in]  pSrc       input pointer
   * @param[out] pDst       output pointer
   * @param[in]  blockSize  number of samples to process
   */
  void adi_sharcfx_copy_f32(
  const float32_t * pSrc,
        float32_t * pDst,
        uint32_t blockSize);

 
 
  /**
   * @brief  Fills a constant value into a floating-point vector.
   * @param[in]  value      input value to be filled
   * @param[out] pDst       output pointer
   * @param[in]  blockSize  number of samples to process
   */
  void adi_sharcfx_fill_f32(
        float32_t value,
        float32_t * pDst,
        uint32_t blockSize);



#ifdef   __cplusplus
}
#endif

#endif /* ifndef _SUPPORT_FUNCTIONS_H_ */
