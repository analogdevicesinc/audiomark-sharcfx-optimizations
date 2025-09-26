/* ----------------------------------------------------------------------
 * Title:        adi_sharcfx_copy_f32.c
 * Description:  Copies the elements of a floating-point vector
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

/**
  @brief         Copies the elements of a floating-point vector.
  @param[in]     pSrc       points to input vector
  @param[out]    pDst       points to output vector
  @param[in]     blockSize  number of samples in each vector
  @return        none
 */

void adi_sharcfx_copy_f32(
  const float32_t * pSrc,
        float32_t * pDst,
        uint32_t blockSize)
{
  
  uint32_t blkCnt = blockSize;
  while (blkCnt > 0U)
  {
    *pDst++ = *pSrc++; // Copy and store result in destination buffer
    blkCnt--;  // Decrement loop counter
  }
}


