/* ----------------------------------------------------------------------
 * Title:        adi_sharcfx_rfft_fast_init_f32.c
 * Description:  Split Radix Decimation in Frequency CFFT Floating point processing function
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



#include <libs/DSP/Include/adi_sharcfx_common_tables.h>
#include <libs/DSP/Include/dsp/transform_functions.h>

/**
  @ingroup groupTransforms
 */

/**
  @addtogroup RealFFT
  @{
 */

#if !defined(ARM_DSP_CONFIG_TABLES) || defined(ARM_ALL_FFT_TABLES) || (defined(ARM_TABLE_TWIDDLECOEF_F32_16) && (defined(ARM_TABLE_BITREVIDX_FLT_16) || defined(ARM_TABLE_BITREVIDX_FXT_16)) && defined(ARM_TABLE_TWIDDLECOEF_RFFT_F32_32))

/**
  @private
  @brief         Initialization function for the 32pt floating-point real FFT.
  @param[in,out] S  points to an adi_sharcfx_rfft_fast_instance_f32 structure
  @return        execution status
                   - \ref ARM_MATH_SUCCESS        : Operation successful
                   - \ref ARM_MATH_ARGUMENT_ERROR : an error is detected
 */

static adi_sharcfx_status adi_sharcfx_rfft_32_fast_init_f32( adi_sharcfx_rfft_fast_instance_f32 * S ) {

  adi_sharcfx_status status;

  if( !S ) return ADI_SHARCFX_MATH_ARGUMENT_ERROR;

  status=adi_sharcfx_cfft_init_f32(&(S->Sint),16);
  if (status != ADI_SHARCFX_MATH_SUCCESS)
  {
    return(status);
  }

  S->fftLenRFFT = 32U;
  S->pTwiddleRFFT    = (float32_t *) twiddleCoef_rfft_32;

  return ADI_SHARCFX_MATH_SUCCESS;
}
#endif 

#if !defined(ARM_DSP_CONFIG_TABLES) || defined(ARM_ALL_FFT_TABLES) || (defined(ARM_TABLE_TWIDDLECOEF_F32_32) && (defined(ARM_TABLE_BITREVIDX_FLT_32) || defined(ARM_TABLE_BITREVIDX_FXT_32)) && defined(ARM_TABLE_TWIDDLECOEF_RFFT_F32_64))

/**
  @private
  @brief         Initialization function for the 64pt floating-point real FFT.
  @param[in,out] S  points to an adi_sharcfx_rfft_fast_instance_f32 structure
  @return        execution status
                   - \ref ARM_MATH_SUCCESS        : Operation successful
                   - \ref ARM_MATH_ARGUMENT_ERROR : an error is detected
 */

static adi_sharcfx_status adi_sharcfx_rfft_64_fast_init_f32( adi_sharcfx_rfft_fast_instance_f32 * S ) {

  adi_sharcfx_status status;

  if( !S ) return ADI_SHARCFX_MATH_ARGUMENT_ERROR;

  status=adi_sharcfx_cfft_init_f32(&(S->Sint),32);
  if (status != ADI_SHARCFX_MATH_SUCCESS)
  {
    return(status);
  }
  S->fftLenRFFT = 64U;

  S->pTwiddleRFFT    = (float32_t *) twiddleCoef_rfft_64;

  return ADI_SHARCFX_MATH_SUCCESS;
}
#endif 

#if !defined(ARM_DSP_CONFIG_TABLES) || defined(ARM_ALL_FFT_TABLES) || (defined(ARM_TABLE_TWIDDLECOEF_F32_64) && (defined(ARM_TABLE_BITREVIDX_FLT_64) || defined(ARM_TABLE_BITREVIDX_FXT_64)) && defined(ARM_TABLE_TWIDDLECOEF_RFFT_F32_128))

/**
  @private
  @brief         Initialization function for the 128pt floating-point real FFT.
  @param[in,out] S  points to an adi_sharcfx_rfft_fast_instance_f32 structure
  @return        execution status
                   - \ref ARM_MATH_SUCCESS        : Operation successful
                   - \ref ARM_MATH_ARGUMENT_ERROR : an error is detected
 */

static adi_sharcfx_status adi_sharcfx_rfft_128_fast_init_f32( adi_sharcfx_rfft_fast_instance_f32 * S ) {

  adi_sharcfx_status status;

  if( !S ) return ADI_SHARCFX_MATH_ARGUMENT_ERROR;

  status=adi_sharcfx_cfft_init_f32(&(S->Sint),64);
  if (status != ADI_SHARCFX_MATH_SUCCESS)
  {
    return(status);
  }
  S->fftLenRFFT = 128;

  S->pTwiddleRFFT    = (float32_t *) twiddleCoef_rfft_128;

  return ADI_SHARCFX_MATH_SUCCESS;
}
#endif 

#if !defined(ARM_DSP_CONFIG_TABLES) || defined(ARM_ALL_FFT_TABLES) || (defined(ARM_TABLE_TWIDDLECOEF_F32_128) && (defined(ARM_TABLE_BITREVIDX_FLT_128) || defined(ARM_TABLE_BITREVIDX_FXT_128)) && defined(ARM_TABLE_TWIDDLECOEF_RFFT_F32_256))

/**
  @private
  @brief         Initialization function for the 256pt floating-point real FFT.
  @param[in,out] S  points to an adi_sharcfx_rfft_fast_instance_f32 structure
  @return        execution status
                   - \ref ARM_MATH_SUCCESS        : Operation successful
                   - \ref ARM_MATH_ARGUMENT_ERROR : an error is detected
*/

static adi_sharcfx_status adi_sharcfx_rfft_256_fast_init_f32( adi_sharcfx_rfft_fast_instance_f32 * S ) {

  adi_sharcfx_status status;

  if( !S ) return ADI_SHARCFX_MATH_ARGUMENT_ERROR;

  status=adi_sharcfx_cfft_init_f32(&(S->Sint),128);
  if (status != ADI_SHARCFX_MATH_SUCCESS)
  {
    return(status);
  }
  S->fftLenRFFT = 256U;

  S->pTwiddleRFFT    = (float32_t *) twiddleCoef_rfft_256;

  return ADI_SHARCFX_MATH_SUCCESS;
}
#endif 

#if !defined(ARM_DSP_CONFIG_TABLES) || defined(ARM_ALL_FFT_TABLES) || (defined(ARM_TABLE_TWIDDLECOEF_F32_256) && (defined(ARM_TABLE_BITREVIDX_FLT_256) || defined(ARM_TABLE_BITREVIDX_FXT_256)) && defined(ARM_TABLE_TWIDDLECOEF_RFFT_F32_512))

/**
  @private
  @brief         Initialization function for the 512pt floating-point real FFT.
  @param[in,out] S  points to an adi_sharcfx_rfft_fast_instance_f32 structure
  @return        execution status
                   - \ref ARM_MATH_SUCCESS        : Operation successful
                   - \ref ARM_MATH_ARGUMENT_ERROR : an error is detected
 */

static adi_sharcfx_status adi_sharcfx_rfft_512_fast_init_f32( adi_sharcfx_rfft_fast_instance_f32 * S ) {

  adi_sharcfx_status status;

  if( !S ) return ADI_SHARCFX_MATH_ARGUMENT_ERROR;

  status=adi_sharcfx_cfft_init_f32(&(S->Sint),256);
  if (status != ADI_SHARCFX_MATH_SUCCESS)
  {
    return(status);
  }
  S->fftLenRFFT = 512U;

  S->pTwiddleRFFT    = (float32_t *) twiddleCoef_rfft_512;

  return ADI_SHARCFX_MATH_SUCCESS;
}
#endif 

#if !defined(ARM_DSP_CONFIG_TABLES) || defined(ARM_ALL_FFT_TABLES) || (defined(ARM_TABLE_TWIDDLECOEF_F32_512) && (defined(ARM_TABLE_BITREVIDX_FLT_512) || defined(ARM_TABLE_BITREVIDX_FXT_512)) && defined(ARM_TABLE_TWIDDLECOEF_RFFT_F32_1024))
/**
  @private
  @brief         Initialization function for the 1024pt floating-point real FFT.
  @param[in,out] S  points to an adi_sharcfx_rfft_fast_instance_f32 structure
  @return        execution status
                   - \ref ARM_MATH_SUCCESS        : Operation successful
                   - \ref ARM_MATH_ARGUMENT_ERROR : an error is detected
 */

static adi_sharcfx_status adi_sharcfx_rfft_1024_fast_init_f32( adi_sharcfx_rfft_fast_instance_f32 * S ) {

  adi_sharcfx_status status;

  if( !S ) return ADI_SHARCFX_MATH_ARGUMENT_ERROR;

  status=adi_sharcfx_cfft_init_f32(&(S->Sint),512);
  if (status != ADI_SHARCFX_MATH_SUCCESS)
  {
    return(status);
  }
  S->fftLenRFFT = 1024U;

  S->pTwiddleRFFT    = (float32_t *) twiddleCoef_rfft_1024;

  return ADI_SHARCFX_MATH_SUCCESS;
}
#endif

#if !defined(ARM_DSP_CONFIG_TABLES) || defined(ARM_ALL_FFT_TABLES) || (defined(ARM_TABLE_TWIDDLECOEF_F32_1024) && (defined(ARM_TABLE_BITREVIDX_FLT_1024) || defined(ARM_TABLE_BITREVIDX_FXT_1024)) && defined(ARM_TABLE_TWIDDLECOEF_RFFT_F32_2048))
/**
  @private
  @brief         Initialization function for the 2048pt floating-point real FFT.
  @param[in,out] S  points to an adi_sharcfx_rfft_fast_instance_f32 structure
  @return        execution status
                   - \ref ARM_MATH_SUCCESS        : Operation successful
                   - \ref ARM_MATH_ARGUMENT_ERROR : an error is detected
 */
static adi_sharcfx_status adi_sharcfx_rfft_2048_fast_init_f32( adi_sharcfx_rfft_fast_instance_f32 * S ) {

  adi_sharcfx_status status;

  if( !S ) return ADI_SHARCFX_MATH_ARGUMENT_ERROR;

  status=adi_sharcfx_cfft_init_f32(&(S->Sint),1024);
  if (status != ADI_SHARCFX_MATH_SUCCESS)
  {
    return(status);
  }
  S->fftLenRFFT = 2048U;

  S->pTwiddleRFFT    = (float32_t *) twiddleCoef_rfft_2048;

  return ADI_SHARCFX_MATH_SUCCESS;
}
#endif

#if !defined(ARM_DSP_CONFIG_TABLES) || defined(ARM_ALL_FFT_TABLES) || (defined(ARM_TABLE_TWIDDLECOEF_F32_2048) && (defined(ARM_TABLE_BITREVIDX_FLT_2048) || defined(ARM_TABLE_BITREVIDX_FXT_2048)) && defined(ARM_TABLE_TWIDDLECOEF_RFFT_F32_4096))
/**
  @private
* @brief         Initialization function for the 4096pt floating-point real FFT.
* @param[in,out] S  points to an adi_sharcfx_rfft_fast_instance_f32 structure
  @return        execution status
                   - \ref ARM_MATH_SUCCESS        : Operation successful
                   - \ref ARM_MATH_ARGUMENT_ERROR : an error is detected
 */

static adi_sharcfx_status adi_sharcfx_rfft_4096_fast_init_f32( adi_sharcfx_rfft_fast_instance_f32 * S ) {

  adi_sharcfx_status status;

  if( !S ) return ADI_SHARCFX_MATH_ARGUMENT_ERROR;

  status=adi_sharcfx_cfft_init_f32(&(S->Sint),2048);
  if (status != ADI_SHARCFX_MATH_SUCCESS)
  {
    return(status);
  }
  S->fftLenRFFT = 4096U;

  S->pTwiddleRFFT    = (float32_t *) twiddleCoef_rfft_4096;

  return ADI_SHARCFX_MATH_SUCCESS;
}
#endif 

/**
  @brief         Initialization function for the floating-point real FFT.
  @param[in,out] S       points to an adi_sharcfx_rfft_fast_instance_f32 structure
  @param[in]     fftLen  length of the Real Sequence
  @return        execution status
                   - \ref ARM_MATH_SUCCESS        : Operation successful
                   - \ref ARM_MATH_ARGUMENT_ERROR : <code>fftLen</code> is not a supported length

  @par           Description
                   The parameter <code>fftLen</code> specifies the length of RFFT/CIFFT process.
                   Supported FFT Lengths are 32, 64, 128, 256, 512, 1024, 2048, 4096.
  @par
                   This Function also initializes Twiddle factor table pointer and Bit reversal table pointer.
 */

adi_sharcfx_status adi_sharcfx_rfft_fast_init_f32(
  adi_sharcfx_rfft_fast_instance_f32 * S,
  uint16_t fftLen)
{
	/* Since both cfft and rfftf use same twiddle tables use cfft_init */
	adi_sharcfx_cfft_instance_f32 cfft_inst;
	adi_sharcfx_status status = adi_sharcfx_cfft_init_f32(&cfft_inst, fftLen);
	S->pTwiddleRFFT = cfft_inst.pTwiddle;
	S->fftLenRFFT = fftLen;
	return status;
}

/**
  @} end of RealFFT group
 */
