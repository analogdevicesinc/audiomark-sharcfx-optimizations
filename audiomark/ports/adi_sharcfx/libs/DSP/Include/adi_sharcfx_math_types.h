/******************************************************************************
 * @file     adi_sharcfx_math_types.h
 * @brief    Public header file for  DSP Library
 * $Date:        22 April 2025
 * $Revision:    V.1.0
 *
 * Target Processor: SHARC-FX Processor
 ******************************************************************************/
/*
 * Copyright(c) 2025 Analog Devices, Inc. All Rights Reserved.
 * Copyright (c) 2010-2021 Arm Limited or its affiliates. All rights reserved.
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



#ifndef _ADI_SHARCFX_MATH_TYPES_H_

#define _ADI_SHARCFX_MATH_TYPES_H_

#ifdef   __cplusplus
extern "C"
{
#endif


#if defined(SHARC_FX)
/* sharcfx compiler specific defines */
#include "../../Core/Include/sharcfx_compiler.h"
#else
  #error Unknown compiler
#endif

#include <string.h>
#include <math.h>
#include <float.h>
#include <limits.h>
#include <complex.h>

#ifdef   __cplusplus
}
#endif


#ifdef   __cplusplus
extern "C"
{
#endif

 /**
   * @brief 8-bit fractional data type in 1.7 format.
   */
  typedef int8_t q7_t;

  /**
   * @brief 16-bit fractional data type in 1.15 format.
   */
  typedef int16_t q15_t;

  /**
   * @brief 32-bit fractional data type in 1.31 format.
   */
  typedef int32_t q31_t;

  /**
   * @brief 64-bit fractional data type in 1.63 format.
   */
  typedef int64_t q63_t;

  /**
   * @brief 32-bit floating-point type definition.
   */
  typedef float float32_t;

  /**
   * @brief 64-bit floating-point type definition.
   */
  typedef double float64_t;



#define F64_MAX   ((float64_t)DBL_MAX)
#define F32_MAX   ((float32_t)FLT_MAX)



#define F64_MIN   (-DBL_MAX)
#define F32_MIN   (-FLT_MAX)



#define F64_ABSMAX   ((float64_t)DBL_MAX)
#define F32_ABSMAX   ((float32_t)FLT_MAX)



#define F64_ABSMIN   ((float64_t)0.0)
#define F32_ABSMIN   ((float32_t)0.0)


#define Q31_MAX   ((q31_t)(0x7FFFFFFFL))
#define Q15_MAX   ((q15_t)(0x7FFF))
#define Q7_MAX    ((q7_t)(0x7F))
#define Q31_MIN   ((q31_t)(0x80000000L))
#define Q15_MIN   ((q15_t)(0x8000))
#define Q7_MIN    ((q7_t)(0x80))

#define Q31_ABSMAX   ((q31_t)(0x7FFFFFFFL))
#define Q15_ABSMAX   ((q15_t)(0x7FFF))
#define Q7_ABSMAX    ((q7_t)(0x7F))
#define Q31_ABSMIN   ((q31_t)0)
#define Q15_ABSMIN   ((q15_t)0)
#define Q7_ABSMIN    ((q7_t)0)

  /* Dimension C vector space */
  #define CMPLX_DIM 2

  /**
   * @brief Error status returned by some functions in the library.
   */

  typedef enum
  {
    ADI_SHARCFX_MATH_SUCCESS                 =  0,        /**< No error */
    ADI_SHARCFX_MATH_ARGUMENT_ERROR          = -1,        /**< One or more arguments are incorrect */
    ADI_SHARCFX_MATH_LENGTH_ERROR            = -2,        /**< Length of data buffer is incorrect */
    ADI_SHARCFX_MATH_SIZE_MISMATCH           = -3,        /**< Size of matrices is not compatible with the operation */
    ADI_SHARCFX_MATH_NANINF                  = -4,        /**< Not-a-number (NaN) or infinity is generated */
    ADI_SHARCFX_MATH_SINGULAR                = -5,        /**< Input matrix is singular and cannot be inverted */
    ADI_SHARCFX_MATH_TEST_FAILURE            = -6,        /**< Test Failed */
    ADI_SHARCFX_MATH_DECOMPOSITION_FAILURE   = -7         /**< Decomposition Failed */
  } adi_sharcfx_status;


#ifdef   __cplusplus
}
#endif

#endif /*ifndef _ADI_SHARCFX_MATH_TYPES_H_ */
