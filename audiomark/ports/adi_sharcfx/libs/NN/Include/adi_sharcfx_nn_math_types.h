/*
 * Copyright(c) 2025 Analog Devices, Inc. All Rights Reserved

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
/******************************************************************************
 * Title:     adi_sharcfx_nn_math_types.h
 * Description:   Compiler include and basic types
 * $Date:        22 April 2025
 * $Revision:    V.1.0
 *
 * Target Processor: SHARC-FX Processor
 ******************************************************************************/

#ifndef _ADI_SHARCFX_NN_MATH_TYPES_H_

#define _ADI_SHARCFX_NN_MATH_TYPES_H_

#ifdef __cplusplus
extern "C" {
#endif

#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdint.h>
#include <string.h>

/* Compiler specific diagnostic adjustment */
#if defined(SHARC_FX)
/* sharcfx compiler specific defines */
#include "../../Core/Include/sharcfx_compiler.h"
#else
#error Unknown compiler
#endif

#ifdef __cplusplus
}
#endif


#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Add necessary typedefs
 */

#define NN_Q31_MAX ((int32_t)(0x7FFFFFFFL))
#define NN_Q15_MAX ((int16_t)(0x7FFF))
#define NN_Q7_MAX ((int8_t)(0x7F))
#define NN_Q31_MIN ((int32_t)(0x80000000L))
#define NN_Q15_MIN ((int16_t)(0x8000))
#define NN_Q7_MIN ((int8_t)(0x80))

/**
 * @brief Error status returned by some functions in the library.
 */

typedef enum
{
    ADI_SHARCFX_NN_SUCCESS = 0,        /**< No error */
    ADI_SHARCFX_NN_ARG_ERROR = -1,     /**< One or more arguments are incorrect */
    ADI_SHARCFX_NN_NO_IMPL_ERROR = -2, /**<  No implementation available */
} adi_sharcfx_nn_status;



#ifdef __cplusplus
}
#endif

#endif /*ifndef _ADI_SHARCFX_NN_MATH_TYPES_H_ */
