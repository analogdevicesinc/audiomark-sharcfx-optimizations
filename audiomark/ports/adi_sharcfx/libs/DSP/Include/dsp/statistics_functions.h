/******************************************************************************
 * @file     statistics_functions.h
 * @brief    Public header file for DSP Library
 * $Date:         22 April 2025
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

 
#ifndef _STATISTICS_FUNCTIONS_H_
#define _STATISTICS_FUNCTIONS_H_

#include "../../../DSP/Include/adi_sharcfx_math_types.h"


#ifdef   __cplusplus
extern "C"
{
#endif


/**
 * @brief Maximum value of absolute values of a floating-point vector.
 * @param[in]  pSrc       points to the input buffer
 * @param[in]  blockSize  length of the input vector
 * @param[out] pResult    maximum value returned here
 * @param[out] pIndex     index of maximum value returned here
 */
  void adi_sharcfx_absmax_f32(
  const float32_t * pSrc,
        uint32_t blockSize,
        float32_t * pResult,
        uint32_t * pIndex);


#ifdef   __cplusplus
}
#endif

#endif /* ifndef _STATISTICS_FUNCTIONS_H_ */
