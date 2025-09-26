/* ----------------------------------------------------------------------
 * Title:        adi_sharcfx_absmax_f32.c
 * Description:  Maximum value of absolute values of a floating-point vector
 * $Date:        14 May 2024
 * $Revision:    V.1.0
 *
 * Target Processor: SHARC-FX Processor
 * -------------------------------------------------------------------- */
/*
 * Copyright(c) 2025 Analog Devices, Inc. All Rights Reserved.
 * Copyright (C) 2010-2021 ARM Limited or its affiliates. All rights reserved.
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


#include <libs/DSP/Include/adi_sharcfx_math_types.h>
#include <vector.h>
#include <stdint.h>

/**
  @brief         Maximum value of absolute values of a floating-point vector.
  @param[in]     pSrc       points to the input vector
  @param[in]     blockSize  number of samples in input vector
  @param[out]    pResult    maximum value returned here
  @param[out]    pIndex     index of maximum value returned here
  @return        none
 */

void adi_sharcfx_absmax_f32(
  const float32_t * pSrc,
        uint32_t blockSize,
        float32_t * pResult,
        uint32_t * pIndex)
{
  float32_t maxVal, out;                         //Temporary variables to store the output value.
  uint32_t blkCnt, outIndex;                     //Loop counter

  outIndex = 0U;//Initialise index value to zero.
  out = fabsf(*pSrc++); // Load first input value that act as reference value for comparision
  blkCnt = (blockSize - 1U); //Initialize blkCnt with number of samples
  
  while (blkCnt > 0U)
  {
    maxVal = fabsf(*pSrc++);// Initialize maxVal to the next consecutive values one by one
    if (out < maxVal) //compare for the maximum value
    {
      out = maxVal; // Update the maximum value and it's index
      outIndex = blockSize - blkCnt;
    }
    blkCnt--; // Decrement loop counter
  }
  *pResult = out; // Store the maximum value and it's index into destination pointers
  *pIndex = outIndex;
}
/**
  @} end of AbsMax group
 */
