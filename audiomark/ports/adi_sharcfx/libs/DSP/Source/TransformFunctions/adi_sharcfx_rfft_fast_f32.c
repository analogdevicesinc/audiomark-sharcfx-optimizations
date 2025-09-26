/* ----------------------------------------------------------------------
 * Project:      ADI DSP Library
 * Title:        adi_sharcfx_rfft_fast_f32.c
 * Description:  RFFT & RIFFT Floating point process function
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


#include <libs/DSP/Include/dsp/transform_functions.h>

/**
  @brief         Processing function for the floating-point real FFT.
  @param[in]     S         points to an adi_sharcfx_rfft_fast_instance_f32 structure
  @param[in]     p         points to input buffer (Source buffer is modified by this function.)
  @param[in]     pOut      points to output buffer
  @param[in]     ifftFlag
                   - value = 0: RFFT
                   - value = 1: RIFFT
  @return        none
*/

void adi_sharcfx_rfft_fast_f32(
  const adi_sharcfx_rfft_fast_instance_f32 * S,
  float32_t * p,
  float32_t * pOut,
  uint8_t ifftFlag)
{
	if (ifftFlag)
	{
		complex_float *p_complex; p_complex=(complex_float *)p;
		rifftf(p_complex, pOut,
				(complex_float *)(S->pTwiddleRFFT), 1, (int)(S->fftLenRFFT));
	}
	else
	{
		rfftf((float32_t *)p, (complex_float *)pOut,
				(const complex_float *)(S->pTwiddleRFFT), 1, (int)(S->fftLenRFFT));
	}
	return;
}

/**
* @} end of RealFFT group
*/
