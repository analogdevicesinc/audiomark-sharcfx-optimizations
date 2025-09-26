/* ----------------------------------------------------------------------
 * Title:        adi_sharcfx_const_structs.h
 * Description:  Constant structs that are initialized for user convenience.
 *               For example, some can be given as arguments to the adi_sharcfx_cfft_f32() function.
 * $Date:        22 April 2025
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



#ifndef _ADI_SHARCFX_CONST_STRUCTS_H
#define _ADI_SHARCFX_CONST_STRUCTS_H

#include "../../DSP/Include/adi_sharcfx_common_tables.h"
#include "../../DSP/Include/adi_sharcfx_math_types.h"
#include "../../DSP/Include/dsp/transform_functions.h"

#ifdef   __cplusplus
extern "C"
{
#endif
   extern const adi_sharcfx_cfft_instance_f64 adi_sharcfx_cfft_sR_f64_len16;
   extern const adi_sharcfx_cfft_instance_f64 adi_sharcfx_cfft_sR_f64_len32;
   extern const adi_sharcfx_cfft_instance_f64 adi_sharcfx_cfft_sR_f64_len64;
   extern const adi_sharcfx_cfft_instance_f64 adi_sharcfx_cfft_sR_f64_len128;
   extern const adi_sharcfx_cfft_instance_f64 adi_sharcfx_cfft_sR_f64_len256;
   extern const adi_sharcfx_cfft_instance_f64 adi_sharcfx_cfft_sR_f64_len512;
   extern const adi_sharcfx_cfft_instance_f64 adi_sharcfx_cfft_sR_f64_len1024;
   extern const adi_sharcfx_cfft_instance_f64 adi_sharcfx_cfft_sR_f64_len2048;
   extern const adi_sharcfx_cfft_instance_f64 adi_sharcfx_cfft_sR_f64_len4096;

   extern const adi_sharcfx_cfft_instance_f32 adi_sharcfx_cfft_sR_f32_len16;
   extern const adi_sharcfx_cfft_instance_f32 adi_sharcfx_cfft_sR_f32_len32;
   extern const adi_sharcfx_cfft_instance_f32 adi_sharcfx_cfft_sR_f32_len64;
   extern const adi_sharcfx_cfft_instance_f32 adi_sharcfx_cfft_sR_f32_len128;
   extern const adi_sharcfx_cfft_instance_f32 adi_sharcfx_cfft_sR_f32_len256;
   extern const adi_sharcfx_cfft_instance_f32 adi_sharcfx_cfft_sR_f32_len512;
   extern const adi_sharcfx_cfft_instance_f32 adi_sharcfx_cfft_sR_f32_len1024;
   extern const adi_sharcfx_cfft_instance_f32 adi_sharcfx_cfft_sR_f32_len2048;
   extern const adi_sharcfx_cfft_instance_f32 adi_sharcfx_cfft_sR_f32_len4096;

   extern const adi_sharcfx_cfft_instance_q31 adi_sharcfx_cfft_sR_q31_len16;
   extern const adi_sharcfx_cfft_instance_q31 adi_sharcfx_cfft_sR_q31_len32;
   extern const adi_sharcfx_cfft_instance_q31 adi_sharcfx_cfft_sR_q31_len64;
   extern const adi_sharcfx_cfft_instance_q31 adi_sharcfx_cfft_sR_q31_len128;
   extern const adi_sharcfx_cfft_instance_q31 adi_sharcfx_cfft_sR_q31_len256;
   extern const adi_sharcfx_cfft_instance_q31 adi_sharcfx_cfft_sR_q31_len512;
   extern const adi_sharcfx_cfft_instance_q31 adi_sharcfx_cfft_sR_q31_len1024;
   extern const adi_sharcfx_cfft_instance_q31 adi_sharcfx_cfft_sR_q31_len2048;
   extern const adi_sharcfx_cfft_instance_q31 adi_sharcfx_cfft_sR_q31_len4096;

   extern const adi_sharcfx_cfft_instance_q15 adi_sharcfx_cfft_sR_q15_len16;
   extern const adi_sharcfx_cfft_instance_q15 adi_sharcfx_cfft_sR_q15_len32;
   extern const adi_sharcfx_cfft_instance_q15 adi_sharcfx_cfft_sR_q15_len64;
   extern const adi_sharcfx_cfft_instance_q15 adi_sharcfx_cfft_sR_q15_len128;
   extern const adi_sharcfx_cfft_instance_q15 adi_sharcfx_cfft_sR_q15_len256;
   extern const adi_sharcfx_cfft_instance_q15 adi_sharcfx_cfft_sR_q15_len512;
   extern const adi_sharcfx_cfft_instance_q15 adi_sharcfx_cfft_sR_q15_len1024;
   extern const adi_sharcfx_cfft_instance_q15 adi_sharcfx_cfft_sR_q15_len2048;
   extern const adi_sharcfx_cfft_instance_q15 adi_sharcfx_cfft_sR_q15_len4096;

#ifdef   __cplusplus
}
#endif

#endif

