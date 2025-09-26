/* ----------------------------------------------------------------------
 * Title:        adi_sharcfx_cfft_radix2_f32.c
 * Description:  Radix-2 Decimation in Frequency CFFT & CIFFT Floating point processing function
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

void adi_sharcfx_radix2_butterfly_f32(
        float32_t * pSrc,
        uint32_t fftLen,
  const float32_t * pCoef,
        uint16_t twidCoefModifier);

void adi_sharcfx_radix2_butterfly_inverse_f32(
        float32_t * pSrc,
        uint32_t fftLen,
  const float32_t * pCoef,
        uint16_t twidCoefModifier,
        float32_t onebyfftLen);

extern void adi_sharcfx_bitreversal_f32(
        float32_t * pSrc,
        uint16_t fftSize,
        uint16_t bitRevFactor,
  const uint16_t * pBitRevTab);

/**
  @ingroup groupTransforms
 */

/**
  @addtogroup ComplexFFT
  @{
 */

/**
  @brief         Radix-2 CFFT/CIFFT.
  @deprecated    Do not use this function. It has been superseded by \ref adi_sharcfx_cfft_f32 and will be removed in the future
  @param[in]     S    points to an instance of the floating-point Radix-2 CFFT/CIFFT structure
  @param[in,out] pSrc points to the complex data buffer of size <code>2*fftLen</code>. Processing occurs in-place
  @return        none
 */

void adi_sharcfx_cfft_radix2_f32(
const adi_sharcfx_cfft_radix2_instance_f32 * S,
      float32_t * pSrc)
{

   if (S->ifftFlag == 1U)
   {
      /* Complex IFFT radix-2 */
      adi_sharcfx_radix2_butterfly_inverse_f32(pSrc, S->fftLen, S->pTwiddle,
      S->twidCoefModifier, S->onebyfftLen);
   }
   else
   {
      /* Complex FFT radix-2 */
      adi_sharcfx_radix2_butterfly_f32(pSrc, S->fftLen, S->pTwiddle,
      S->twidCoefModifier);
   }

   if (S->bitReverseFlag == 1U)
   {
      /* Bit Reversal */
      adi_sharcfx_bitreversal_f32(pSrc, S->fftLen, S->bitRevFactor, S->pBitRevTable);
   }

}


/**
  @} end of ComplexFFT group
 */



/* ----------------------------------------------------------------------
 ** Internal helper function used by the FFTs
 ** ------------------------------------------------------------------- */

/**
  brief  Core function for the floating-point CFFT butterfly process.
  param[in,out] pSrc             points to in-place buffer of floating-point data type
  param[in]     fftLen           length of the FFT
  param[in]     pCoef            points to twiddle coefficient buffer
  param[in]     twidCoefModifier twiddle coefficient modifier that supports different size FFTs with the same twiddle factor table
  return        none
 */

void adi_sharcfx_radix2_butterfly_f32(
        float32_t * pSrc,
        uint32_t fftLen,
  const float32_t * pCoef,
        uint16_t twidCoefModifier)
{

        uint32_t i, j, k, l;
        uint32_t n1, n2, ia;
        float32_t xt, yt, cosVal, sinVal;
        float32_t p0, p1, p2, p3;
        float32_t a0, a1;


   n2 = fftLen;

   // loop for stage
   for (k = fftLen; k > 1; k = k >> 1)
   {
      n1 = n2;
      n2 = n2 >> 1;
      ia = 0;

      // loop for groups
      j = 0;
      do
      {
         cosVal = pCoef[ia * 2];
         sinVal = pCoef[(ia * 2) + 1];
         ia += twidCoefModifier;

         // loop for butterfly
         i = j;
         do
         {
            l = i + n2;
            a0 = pSrc[2 * i] + pSrc[2 * l];
            xt = pSrc[2 * i] - pSrc[2 * l];

            yt = pSrc[2 * i + 1] - pSrc[2 * l + 1];
            a1 = pSrc[2 * l + 1] + pSrc[2 * i + 1];

            p0 = xt * cosVal;
            p1 = yt * sinVal;
            p2 = yt * cosVal;
            p3 = xt * sinVal;

            pSrc[2 * i] = a0;
            pSrc[2 * i + 1] = a1;

            pSrc[2 * l]     = p0 + p1;
            pSrc[2 * l + 1] = p2 - p3;

            i += n1;
         } while (i < fftLen);
         j++;
      } while (j < n2);
      twidCoefModifier <<= 1U;
   }

}


void adi_sharcfx_radix2_butterfly_inverse_f32(
        float32_t * pSrc,
        uint32_t fftLen,
  const float32_t * pCoef,
        uint16_t twidCoefModifier,
        float32_t onebyfftLen)
{

        uint32_t i, j, k, l;
        uint32_t n1, n2, ia;
        float32_t xt, yt, cosVal, sinVal;
        float32_t p0, p1, p2, p3;
        float32_t a0, a1;


   n2 = fftLen;

   // loop for stage
   for (k = fftLen; k > 2; k = k >> 1)
   {
      n1 = n2;
      n2 = n2 >> 1;
      ia = 0;

      // loop for groups
      j = 0;
      do
      {
         cosVal = pCoef[ia * 2];
         sinVal = pCoef[(ia * 2) + 1];
         ia = ia + twidCoefModifier;

         // loop for butterfly
         i = j;
         do
         {
            l = i + n2;
            a0 = pSrc[2 * i] + pSrc[2 * l];
            xt = pSrc[2 * i] - pSrc[2 * l];

            yt = pSrc[2 * i + 1] - pSrc[2 * l + 1];
            a1 = pSrc[2 * l + 1] + pSrc[2 * i + 1];

            p0 = xt * cosVal;
            p1 = yt * sinVal;
            p2 = yt * cosVal;
            p3 = xt * sinVal;

            pSrc[2 * i] = a0;
            pSrc[2 * i + 1] = a1;

            pSrc[2 * l]     = p0 - p1;
            pSrc[2 * l + 1] = p2 + p3;

            i += n1;
         } while ( i < fftLen );                    // butterfly loop end
         j++;
      } while ( j < n2 );                      // groups loop end

      twidCoefModifier = twidCoefModifier << 1U;
   }                             // stages loop end

   n1 = n2;
   n2 = n2 >> 1;

   // loop for butterfly
   for (i = 0; i < fftLen; i += n1)
   {
      l = i + n2;

      a0 = pSrc[2 * i] + pSrc[2 * l];
      xt = pSrc[2 * i] - pSrc[2 * l];

      a1 = pSrc[2 * l + 1] + pSrc[2 * i + 1];
      yt = pSrc[2 * i + 1] - pSrc[2 * l + 1];

      p0 = a0 * onebyfftLen;
      p2 = xt * onebyfftLen;
      p1 = a1 * onebyfftLen;
      p3 = yt * onebyfftLen;

      pSrc[2 * i] = p0;
      pSrc[2 * l] = p2;

      pSrc[2 * i + 1] = p1;
      pSrc[2 * l + 1] = p3;
   }                             // butterfly loop end

}
