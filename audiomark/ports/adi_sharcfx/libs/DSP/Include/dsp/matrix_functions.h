/******************************************************************************
 * @file     matrix_functions.h
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

 
#ifndef _MATRIX_FUNCTIONS_H_
#define _MATRIX_FUNCTIONS_H_


#include "../../../DSP/Include/adi_sharcfx_math_types.h"

#ifdef   __cplusplus
extern "C"
{
#endif


  /**
   * @brief Instance structure for the floating-point matrix structure.
   */
  typedef struct
  {
    uint16_t numRows;     /**< number of rows of the matrix.     */
    uint16_t numCols;     /**< number of columns of the matrix.  */
    float32_t *pData;     /**< points to the data of the matrix. */
  } adi_sharcfx_matrix_instance_f32;


  /**
   * @brief Floating-point matrix and vector multiplication
   * @param[in]  pSrcMat  points to the input matrix structure
   * @param[in]  pVec     points to vector
   * @param[out] pDst     points to output vector
   */
void adi_sharcfx_mat_vec_mult_f32(
  const adi_sharcfx_matrix_instance_f32 *pSrcMat, 
  const float32_t *pVec, 
  float32_t *pDst);


#ifdef   __cplusplus
}
#endif

#endif /* ifndef _MATRIX_FUNCTIONS_H_ */
