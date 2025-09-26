
/*
 * Copyright(c) 2025 Analog Devices, Inc. All Rights Reserved.

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

/* ----------------------------------------------------------------------
 * Title:        adi_sharcfx_depthwise_conv_3x3_s8.c
 * Description:  Optimized s8 depthwise convolution function for channel
 *               multiplier of 1 and 3x3 kernel size.
 * $Date:        22 April 2025
 * $Revision:    V.1.0
 *
 * Target Processor:  SHARC-FX Processor
 *
 * -------------------------------------------------------------------- */


#include <sys/platform.h>
#include <stdio.h>

#include <stdio.h>
#include <math.h>
#include <float.h>
#include <complex.h>

#include <matrix.h>
#include <filter.h>
#include <vector.h>

#define NVEC 32

#include "math_fixedpoint_vec.h"
/* Cross-platform data type definitions. */
#include "libdsp_types.h"


#pragma once
#include <inttypes.h>

#ifdef __XTENSA__
#include <xtensa/sim.h>
#include <xtensa/tie/xt_pdxn.h>
#endif

#include <libs/NN/Include/adi_sharcfx_nnfunctions.h>
#include <libs/NN/Include/adi_sharcfx_nnsupportfunctions.h>


#include <stdint.h>
#include <xtensa/tie/xt_pdxn.h>


#include <alloca.h>
#include <stdint.h>
#include <stdio.h>

/* Cross-platform data type definitions. */
#include "libdsp_types.h"


#define TEMP_BUFFER_SIZE 96*96*20

#define MIN(X, Y) (((X) < (Y)) ? (X) : (Y))
#define MAX(X, Y) (((X) > (Y)) ? (X) : (Y))

#define ROUNDING_MODE                    1
//Temp scratch buffer used to introduce padding in depthwise convolution
int8_t pTemp[TEMP_BUFFER_SIZE];
/* Add your custom header content here */

//---------DATATYPES-------------
//input data type is int8, bias/mult are int32
typedef int8_t ADI_ACTIVATION_DATATYPE;
typedef int8_t ADI_WEIGHTS_DATATYPE;
typedef int32_t ADI_CONV_BUFF_DATATYPE;
typedef uint32_t ADI_QMULTIPLIER_DATATYPE;
typedef int16_t ADI_16BIT_CONV_BUFF_DATATYPE;

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define ABS(a) (((a) < 0) ? -(a) : (a))

// Unaligned 16-way by 8b load
#define LOADU_2MX8(dst,ptr,inc) \
	ptr ## a = PDX_LA_2MX8_PP (ptr); \
	PDX_LA16_2MX8_XP (dst, ptr ## a, ptr, inc)

// Unaligned 32-way by 8b load
#define LOADU_4MX8(dst,ptr,inc) \
	ptr ## a = PDX_LA_4MX8_PP (ptr); \
	PDX_LA_4MX8_XP (dst, ptr ## a, ptr, inc)


/* Cross-platform data type definitions. */
#include "libdsp_types.h"

#define USE_FULL_TF_ROUND_IN_REQUANT


#ifdef USE_TIE
// Use a TIE instruction to speed up this operation
// Do the requantization bias, multiply, shift, pack, ozp, sat8 and save
#define REQUANT(i) \
        /* Multiply */ \
        vScaled = vDot ## i * vScale ## i; \
        /* Shift for fixed point, st */ \
        vScaled = PDX_SLS_MX80 (vScaled, vShift ## i); \
        /* Saturate, round, remove low 32 bits, add the output zero point, saturate to 8 bits */ \
        vOut = ADI_PACK_ADD_SAT8_MX32 (vScaled, vOutZeroPoint); \
        /* Save the remaining bytes in this pass */ \
        PDX_SAV32_MX8_XP (vOut, pvOuta, pvOut, c8_left); \
        c8_left -= PDX_M
#elif defined (USE_FULL_TF_ROUND_IN_REQUANT)
// Do the requantization bias, multiply, shift, pack, ozp, sat8 and save
#define REQUANT(i) \
        /* Multiply */ \
        vScaled = vDot ## i * vScale ## i; \
        /* Handle left shifts */ \
        vScaled1 = PDX_SLS_MX80 (vScaled, vShift ## i); \
        /* Get left shift outputs and cut off low bits for right shifts */ \
        vOut_ge_zero = PDX_PACKQSRV_MX80 (vScaled1, 0); \
        /* Put back for right shifts */ \
		vScaled2 = vOut_ge_zero; \
		/* Move up to Q32 and do the right shifts */ \
		vScaled3 = PDX_SLS_MX80 (vScaled2, v32_m_shift ## i); \
        /* Saturate, round, remove low 32 bits */ \
        vOut_lt_zero = PDX_PACKQSRV_MX80 (vScaled3, 1); \
        /* Select between outputs for shift >= 0 and < 0 */ \
        vOut = PDX_MOV_MX32_T (vOut_ge_zero, vOut_lt_zero, en_sh_ge_zero ## i); \
        /* Add yet another bias */ \
        vOut += vOutZeroPoint; \
        /* Saturate to 8 bits */ \
        vOut = PDX_MIN_MX32 (PDX_MAX_MX32 (vOut, -128), 127); \
        /* Save the remaining bytes in this pass */ \
        PDX_SAV32_MX8_XP (vOut, pvOuta, pvOut, c8_left); \
        c8_left -= PDX_M
#else
// Do the requantization bias, multiply, shift, pack, ozp, sat8 and save
#define REQUANT(i) \
        /* Multiply */ \
        vScaled = vDot ## i * vScale ## i; \
        /* Shift for fixed point */ \
        vScaled = PDX_SLS_MX80 (vScaled, vShift ## i); \
        /* Saturate, round, remove low 32 bits */ \
        vOut = PDX_PACKQSRV_MX80 (vScaled, 0); \
        /* Add yet another bias */ \
        vOut += vOutZeroPoint; \
        /* Saturate to 8 bits */ \
        vOut = PDX_MIN_MX32 (PDX_MAX_MX32 (vOut, -128), 127); \
        /* Save the remaining bytes in this pass */ \
        PDX_SAV32_MX8_XP (vOut, pvOuta, pvOut, c8_left); \
        c8_left -= PDX_M
#endif


/*
 * Optimized s8 depthwise convolution function with constraint that
 * in_channel == out_channel and kernel_x == kernel_y == 3 with pads at most 1
 *
 *  Refer prototype header file for details.
 *
 */
adi_sharcfx_nn_status adi_sharcfx_depthwise_conv_3x3_s8(const nn_context *ctx,
                                              const nn_dw_conv_params *dw_conv_params,
                                              const nn_per_channel_quant_params *quant_params,
                                              const nn_dims *input_dims,
                                              const int8_t *pInputBuffer,  /* [in] input buffer of size [nInputHeight][nInputWidth][nChannels] */
                                              const nn_dims *filter_dims,
											  const int8_t *pWeightsBuffer,   /* [in] weight buffer of size [9][nChannels]*/
                                              const nn_dims *bias_dims,
                                              const int32_t *pBiasBuffer,  /* [in] bias buffer of size [nChannels]  */
                                              const nn_dims *output_dims,
                                              int8_t *pOutputBuffer ) /* [out] output buffer of size [nInputHeight][nInputWidth][nChannels] */
{
		    (void)ctx;
		    (void)bias_dims;

		    const int32_t nInputWidth = input_dims->w; /* [in] sets of of input buffer channels in one row */
		    const int32_t nInputHeight = input_dims->h;  /* [in] number of rows  */
		    const int32_t nInChannels = input_dims->c; /* [in] number of input channels */
		    const int32_t nOutChannels = output_dims->c;  /* [in] number of output channels, must equal nInChannels */
		    const int32_t nTotalPaddingWidth = 2*(dw_conv_params->padding.w); /* [in] unused - padding is always 1 pixel of 0 */
		    const int32_t nTotalPaddingHeight = 2*(dw_conv_params->padding.h);/* [in] unused - padding is always 1 pixel of 0 */
		    const int32_t nStrideWidth = dw_conv_params->stride.w;  /* [in] step from one set of input channels to the next */
		    const int32_t nStrideHeight = dw_conv_params->stride.h; /* [in] step from one row to the next */
		    const int32_t *pQuantizedShift = quant_params->shift;  /* [in] shift buffer of size [nChannels] */
		    const int32_t *pQuantizedMultiplier = quant_params->multiplier; /* [in] multiplier buffer of size [nChannels] */
		    const int32_t output_x = output_dims->w;
		    const int32_t output_y = output_dims->h;
		    const int32_t nOutZeroPoint = dw_conv_params->output_offset;  /* [in] post-scaling zero-point bias for all outputs */
		    const int32_t nInZeroPoint  = dw_conv_params->input_offset; /* [in] zero-point bias for all inputs */
		    const int32_t ACT_MIN = dw_conv_params->activation.min;
		    const int32_t ACT_MAX = dw_conv_params->activation.max;
		    const int32_t nKernelSizeWidth = filter_dims->w; /* [in] Kernel width, must be 3 */
		    const int32_t nKernelSizeHeight = filter_dims->h; /* [in] Kernel height, must be 3 */

			int32_t nDepthMult = 0; /* [in] if 1, use 1 input channel for all outputs.  Must be 0. */

			 /* Assertions for this particular routine */
			    ASSERT (nDepthMult != 1);
			    ASSERT (nOutChannels == nInChannels);
			    ASSERT (nKernelSizeWidth == 3);
			    ASSERT (nKernelSizeHeight == 3);
			    ASSERT (nTotalPaddingWidth == 2 || nTotalPaddingWidth == 0);
			    ASSERT (nTotalPaddingHeight == 2 || nTotalPaddingHeight == 0);
			    /* Generic good-practice assertions */
			    // No NULL pointers
			    ASSERT (pInputBuffer && pOutputBuffer && pWeightsBuffer && pBiasBuffer && pQuantizedMultiplier && pQuantizedShift);
			    // All positive sizes
			    ASSERT (nInputWidth > 0 && nInputHeight > 0 && nInChannels > 0);
			    // All positive strides
			    ASSERT (nStrideWidth > 0 && nStrideHeight > 0);

			    int c, h, ho, w;
			    int WS = nStrideWidth;
			    int HS = nStrideHeight;
			    int W = nInputWidth;
			    int H = nInputHeight;
			    int C = nInChannels;

			    xb_vec2Mx16 vInZeroPoint = nInZeroPoint;
			    xb_vecMx8 vOutZeroPoint = nOutZeroPoint;

			    // Work on a vector's-worth of channels at a time
			    for (c = 0; c < C; c += (2*PDX_M)) {
			        // Find how many channels are left in this pass
			        int c_left = MIN ((2*PDX_M), C - c);
			        int8_t *pIn  = (int8_t *) &pInputBuffer[c];
			        int8_t *pOut = &pOutputBuffer[c];
			        // Keep the 9 weights in registers
			        xb_vec2Mx8 *pW   = (xb_vec2Mx8 *) &pWeightsBuffer[c];
			        valign pWa;
			        xb_vec2Mx16 vW0, vW1, vW2, vW3, vW4, vW5, vW6, vW7, vW8;
			        LOADU_2MX8 (vW0, pW, C);
			        LOADU_2MX8 (vW1, pW, C);
			        LOADU_2MX8 (vW2, pW, C);
			        LOADU_2MX8 (vW3, pW, C);
			        LOADU_2MX8 (vW4, pW, C);
			        LOADU_2MX8 (vW5, pW, C);
			        LOADU_2MX8 (vW6, pW, C);
			        LOADU_2MX8 (vW7, pW, C);
			        LOADU_2MX8 (vW8, pW, C);

			        // Get the bias for requant
			        xb_vecMx32 *pBias  = (xb_vecMx32 *) &pBiasBuffer[c];
			        valign pBiasa = PDX_LA_MX32_PP (pBias);
			        xb_vecMx32 vBias0, vBias1;
			        PDX_LA_MX32_IP (vBias0, pBiasa, pBias);
			        PDX_LA_MX32_IP (vBias1, pBiasa, pBias);

			        // Get the scale for requant
			        xb_vecMx32 *pScale  = (xb_vecMx32 *) &pQuantizedMultiplier[c];
			        valign pScalea = PDX_LA_MX32_PP (pScale);
			        xb_vecMx32 vScale0, vScale1;
			        PDX_LA_MX32_IP (vScale0, pScalea, pScale);
			        PDX_LA_MX32_IP (vScale1, pScalea, pScale);

			        // Get the shifts for requant
			        xb_vecMx32 *pShift  = (xb_vecMx32 *) &pQuantizedShift[c];
			        valign pShifta = PDX_LA_MX32_PP (pShift);
			        xb_vecMx32 vShift0, vShift1;
			        PDX_LA_MX32_IP (vShift0, pShifta, pShift);
			        PDX_LA_MX32_IP (vShift1, pShifta, pShift);
			#ifdef USE_FULL_TF_ROUND_IN_REQUANT
			        vboolM en_sh_ge_zero0 = PDX_GE_MX32 (vShift0, 0);
			        vboolM en_sh_ge_zero1 = PDX_GE_MX32 (vShift1, 0);
			        xb_vecMx32 vOut_ge_zero, vOut_lt_zero;
			        xb_vecMx32 v32 = 32;
			        xb_vecMx32 v32_m_shift0 = v32 + vShift0;
			        xb_vecMx32 v32_m_shift1 = v32 + vShift1;
			        // Add one to convert Q31 to Q32, use 1 if shift < 0
			        xb_vecMx32 vOne = 1;
			        vShift0 = PDX_MOV_MX32_T (vShift0 + vOne, 1, en_sh_ge_zero0);
			        vShift1 = PDX_MOV_MX32_T (vShift1 + vOne, 1, en_sh_ge_zero0);
			        xb_vecMx80 vScaled1, vScaled2, vScaled3;
			#endif

			        // Enables for rows and columns that handle padded pixels
			        vbool2M en_row0, en_row1, en_row2, en_col0, en_col2;
			        xb_vec2Mx16 vh = 0;
			        xb_vec2Mx16 vHS = HS;
			        xb_vec2Mx16 vHm1 = H - 1;

			        en_row1 = PDX_EQ_2MX16 (vh, vh);   // Always true

			        xb_vecMx80 vScaled;
			        int c8_left;

			        // Step down the rows by the row stride
			        for (h = ho = 0; h < H; h += HS, ho++) {

			            en_row0 = PDX_GE_2MX16 (vh, vHS);  // Always clear on first pass
			            en_row2 = PDX_LT_2MX16 (vh, vHm1);  // True when h < H-HS, before the last row
			            vh += vHS;  // Increment for next pass

			            // Accumulator
			            xb_vec2Mx40 vDot;

			            xb_vec2Mx16 vw = 0;
			            xb_vec2Mx16 vWS = WS;
			            xb_vec2Mx16 vWm1 = W - 1;

			            // Line pointers.  Set to one before the start of the first, second, and third row
			            xb_vec2Mx8 * pvIn0 = (xb_vec2Mx8 *) &pIn[(h - 1) * W * C - C];
			            xb_vec2Mx8 * pvIn1 = (xb_vec2Mx8 *) &pIn[(h)     * W * C - C];
			            xb_vec2Mx8 * pvIn2 = (xb_vec2Mx8 *) &pIn[(h + 1) * W * C - C];
			            valign pvIn0a, pvIn1a, pvIn2a;

			            // Output pointer
			            xb_vecMx8 *pvOut = (xb_vecMx8 *) &pOut[ho*W/WS*C];
			            valign pvOuta = PDX_Z_ALIGN ();
			            // Output data
			            xb_vecMx32 vOut;

			            // Input data for each of the 9 pixels
			            xb_vec2Mx16 vIn0, vIn1, vIn2, vIn3, vIn4, vIn5, vIn6, vIn7, vIn8;
			            // Inc from the left column to the next right column
			            int in_inc = C * WS - C - C;

			            if (WS > 1) {
			                // With a width stride of more than one, the pixels cannot be
			                // assumed to be adjacent, and so all must be read in on each pass
			#pragma concurrent
			                for (w = 0; w < (W + WS - 1)/WS; w ++) {

			                    en_col0 = PDX_GE_2MX16 (vw, vWS);
			                    // True when w < W-WS, before the last column
			                    en_col2 = PDX_LT_2MX16 (vw, vWm1);
			                    // Increment to next column
			                    vw += vWS;

			                    // Do the 3x3 convolution

			                    // Init accumulator
			                    vDot = PDX_CVT40_MX32D (vBias1, vBias0);

			                    // Get the left column
			                    LOADU_2MX8 (vIn0, pvIn0, C);
			                    vIn0 += vInZeroPoint;
			                    LOADU_2MX8 (vIn3, pvIn1, C);
			                    vIn3 += vInZeroPoint;
			                    LOADU_2MX8 (vIn6, pvIn2, C);
			                    vIn6 += vInZeroPoint;

			                    // Add in left column
			                    PDX_MULAW_2MX16_T (vDot, vW0, vIn0, en_row0 & en_col0);
			                    PDX_MULAW_2MX16_T (vDot, vW3, vIn3,           en_col0);
			                    PDX_MULAW_2MX16_T (vDot, vW6, vIn6, en_row2 & en_col0);

			                    // Get the middle column
			                    LOADU_2MX8 (vIn1, pvIn0, C);
			                    vIn1 += vInZeroPoint;
			                    LOADU_2MX8 (vIn4, pvIn1, C);
			                    vIn4 += vInZeroPoint;
			                    LOADU_2MX8 (vIn7, pvIn2, C);
			                    vIn7 += vInZeroPoint;

			                    // Add in middle column
			                    PDX_MULAW_2MX16_T (vDot, vW1, vIn1, en_row0);
			                    PDX_MULAW_2MX16_T (vDot, vW4, vIn4, en_row1);
			                    PDX_MULAW_2MX16_T (vDot, vW7, vIn7, en_row2);

			                    // Get the right column
			                    LOADU_2MX8 (vIn2, pvIn0, in_inc);
			                    vIn2 += vInZeroPoint;
			                    LOADU_2MX8 (vIn5, pvIn1, in_inc);
			                    vIn5 += vInZeroPoint;
			                    LOADU_2MX8 (vIn8, pvIn2, in_inc);
			                    vIn8 += vInZeroPoint;

			                    // Add in right column
			                    PDX_MULAW_2MX16_T (vDot, vW2, vIn2, en_row0 & en_col2);
			                    PDX_MULAW_2MX16_T (vDot, vW5, vIn5,           en_col2);
			                    PDX_MULAW_2MX16_T (vDot, vW8, vIn8, en_row2 & en_col2);

			                    // Move from the 16-way 40b accumulator to two 8-way 32b vectors
			                    xb_vecMx32 vDot0, vDot1;
			                    PDX_CVT32D_2MX40 (vDot1, vDot0, vDot);

			                    // Do the requantization on the two halfs of the accumulator

			                    // Save only the channels left in the set
			                    c8_left = c_left;

			                    // Requantize each half
			                    REQUANT(0);
			                    REQUANT(1);
			                    PDX_SAPOS_MX8_FP (pvOuta, pvOut);  // Flush out last bytes
			                    pvOut = (xb_vecMx8 *) ((int8_t *) pvOut + C - c_left);  // Inc to next set of channels
			                }
			            } else {
			                // With a width stride of one, the code
			                // Can re-use two columns from the previous pass, which
			                // saves loads and zero-point adds

			                // Prefetch the first two columns
			                // Get the left column
			                LOADU_2MX8 (vIn0, pvIn0, C);
			                vIn0 += vInZeroPoint;
			                LOADU_2MX8 (vIn3, pvIn1, C);
			                vIn3 += vInZeroPoint;
			                LOADU_2MX8 (vIn6, pvIn2, C);
			                vIn6 += vInZeroPoint;
			                if (1) {
			                    // Clear if padding
			                    vIn0 = vIn3 = vIn6 = 0;
			                }

			                // Get the middle column
			                LOADU_2MX8 (vIn1, pvIn0, C);
			                vIn1 += vInZeroPoint;
			                LOADU_2MX8 (vIn4, pvIn1, C);
			                vIn4 += vInZeroPoint;
			                LOADU_2MX8 (vIn7, pvIn2, C);
			                vIn7 += vInZeroPoint;

			                in_inc = C;
			                int out_inc = C - c_left;
			#pragma concurrent
			                for (w = 0; w < W; w ++) {
			                    // True when w < W-WS, before the last column
			                    en_col2 = PDX_LT_2MX16 (vw, vWm1);
			                    // Increment to next column
			                    vw += vWS;

			                    // Do the 3x3 convolution

			                    // Init accumulator
			                    vDot = PDX_CVT40_MX32D (vBias1, vBias0);

			                    // Add in left column
			                    PDX_MULAW_2MX16_T (vDot, vW0, vIn0, en_row0);
			                    PDX_MULAW_2MX16_T (vDot, vW3, vIn3, en_row1);
			                    PDX_MULAW_2MX16_T (vDot, vW6, vIn6, en_row2);

			                    // Add in middle column
			                    PDX_MULAW_2MX16_T (vDot, vW1, vIn1, en_row0);
			                    PDX_MULAW_2MX16_T (vDot, vW4, vIn4, en_row1);
			                    PDX_MULAW_2MX16_T (vDot, vW7, vIn7, en_row2);

			                    // Get the right column
			                    LOADU_2MX8 (vIn2, pvIn0, in_inc);
			                    vIn2 += vInZeroPoint;
			                    LOADU_2MX8 (vIn5, pvIn1, in_inc);
			                    vIn5 += vInZeroPoint;
			                    LOADU_2MX8 (vIn8, pvIn2, in_inc);
			                    vIn8 += vInZeroPoint;

			                    // Shift over the columns already used
			                    vIn0 = vIn1; vIn1 = vIn2;
			                    vIn3 = vIn4; vIn4 = vIn5;
			                    vIn6 = vIn7; vIn7 = vIn8;

			                    // Add in right column
			                    PDX_MULAW_2MX16_T (vDot, vW2, vIn2, en_row0 & en_col2);
			                    PDX_MULAW_2MX16_T (vDot, vW5, vIn5,           en_col2);
			                    PDX_MULAW_2MX16_T (vDot, vW8, vIn8, en_row2 & en_col2);

			                    // Move from the 16-way 40b accumulator to two 8-way 32b vectors
			                    xb_vecMx32 vDot0, vDot1;
			                    PDX_CVT32D_2MX40 (vDot1, vDot0, vDot);

			                    // Do the requantization on the two halfs of the accumulator

			                    // Save only the channels left in the set
			                    c8_left = c_left;

			                    // Requantize each half
			                    REQUANT(0);
			                    REQUANT(1);
			                    PDX_SAPOS_MX8_FP (pvOuta, pvOut);  // Flush out last bytes
			                    pvOut = (xb_vecMx8 *) ((int8_t *) pvOut + out_inc);  // Inc to next set of channels
			                }
			            }
			        }
			    }
	return ADI_SHARCFX_NN_SUCCESS;
}
