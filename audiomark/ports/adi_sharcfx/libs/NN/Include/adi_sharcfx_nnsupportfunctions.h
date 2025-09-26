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
 * Title:        adi_sharcfx_nnsupportfunctions.h
 * Description:  Public header file of support functions for NN Library
 *
 * $Date:        22 April 2025
 * $Revision:    V.1.0
 *
 * Target Processor:  SHARC-FX Processor
 * -------------------------------------------------------------------- */
 
#ifndef _ADI_SHARCFX_NNSUPPORTFUNCTIONS_H_
#define _ADI_SHARCFX_NNSUPPORTFUNCTIONS_H_

#include <stdbool.h>

#include "../../NN/Include/adi_sharcfx_nn_math_types.h"
#include "../../NN/Include/adi_sharcfx_nn_types.h"

#ifdef __cplusplus
extern "C" {
#endif

//extern int global_count;
//extern int64_t total_cycles;

#define LEFT_SHIFT(_shift) (_shift > 0 ? _shift : 0)
#define RIGHT_SHIFT(_shift) (_shift > 0 ? 0 : -_shift)
#define MASK_IF_ZERO(x) (x) == 0 ? ~0 : 0
#define MASK_IF_NON_ZERO(x) (x) != 0 ? ~0 : 0
#define SELECT_USING_MASK(mask, a, b) ((mask) & (a)) ^ (~(mask) & (b))

#define MAX(A, B) ((A) > (B) ? (A) : (B))
#define MIN(A, B) ((A) < (B) ? (A) : (B))
#define CLAMP(x, h, l) MAX(MIN((x), (h)), (l))
#define REDUCE_MULTIPLIER(_mult) ((_mult < 0x7FFF0000) ? ((_mult + (1 << 15)) >> 16) : 0x7FFF)

// Number of channels processed in a block for DW Conv(MVE)
// Requirement: Greater than 0 & less than 128
// This can be fine tuned to match number of input channels for best performance.
// A layer with lower number of channels than CH_IN_BLOCK_MVE will result in higher
// scratch buffer usage and a layer with higher number of channels than CH_IN_BLOCK_MVE
// will result in lower scratch buffer usage.
#define CH_IN_BLOCK_MVE (124)

/**
 * @brief definition to pack four 8 bit values.
 */
#define PACK_S8x4_32x1(v0, v1, v2, v3)                                                                                 \
    ((((int32_t)(v0) << 0) & (int32_t)0x000000FF) | (((int32_t)(v1) << 8) & (int32_t)0x0000FF00) |                     \
     (((int32_t)(v2) << 16) & (int32_t)0x00FF0000) | (((int32_t)(v3) << 24) & (int32_t)0xFF000000))

/**
 * @brief definition to pack two 16 bit values.
 */
#define PACK_Q15x2_32x1(v0, v1) (((int32_t)v0 & (int32_t)0xFFFF) | ((int32_t)v1 << 16))

/**
 * @brief Union for SIMD access of q31/s16/s8 types
 */
union adi_sharcfx_nnword
{
    int32_t word;
    /**< q31 type */
    int16_t half_words[2];
    /**< s16 type */
    int8_t bytes[4];
    /**< s8 type */
};

/**
 * @brief Union for data type long long
 */
struct adi_sharcfx_nn_double
{
    uint32_t low;
    int32_t high;
};

union adi_sharcfx_nn_long_long
{
    int64_t long_long;
    struct adi_sharcfx_nn_double word;
};

/**
 * @defgroup groupSupport Private
 *
 * Internal Support functions. Not intended to be called direclty by a NN user.
 *
 */

/**
 * @defgroup supportConversion Data Conversion
 *
 * Perform data type conversion in-between neural network operations
 *
 */

/**
 * @brief Converts the elements from a s8 vector to a s16 vector with an added offset
 * @param[in]    src        pointer to the s8 input vector
 * @param[out]   dst        pointer to the s16 output vector
 * @param[in]    block_size length of the input vector
 * @param[in]    offset     s8 offset to be added to each input vector element.
 *
 * \par Description:
 *
 * The equation used for the conversion process is:
 *
 * <pre>
 *  dst[n] = (int16_t) src[n] + offset;   0 <= n < block_size.
 * </pre>
 *
 */
void adi_sharcfx_q7_to_q15_with_offset(const int8_t *src, int16_t *dst, uint32_t block_size, int16_t offset);




/**
 * @brief General Matrix-multiplication function with per-channel requantization.
 *        This function assumes:
 *        - LHS input matrix NOT transposed (nt)
 *        - RHS input matrix transposed (t)
 *
 *  @note This operation also performs the broadcast bias addition before the requantization
 *
 * @param[in]  lhs                Pointer to the LHS input matrix
 * @param[in]  rhs                Pointer to the RHS input matrix
 * @param[in]  bias               Pointer to the bias vector. The length of this vector is equal to the number of
 * output columns (or RHS input rows)
 * @param[out] dst                Pointer to the output matrix with "m" rows and "n" columns
 * @param[in]  dst_multipliers    Pointer to the multipliers vector needed for the per-channel requantization.
 *                                The length of this vector is equal to the number of output columns (or RHS input
 * rows)
 * @param[in]  dst_shifts         Pointer to the shifts vector needed for the per-channel requantization. The length
 * of this vector is equal to the number of output columns (or RHS input rows)
 * @param[in]  lhs_rows           Number of LHS input rows
 * @param[in]  rhs_rows           Number of RHS input rows
 * @param[in]  rhs_cols           Number of LHS/RHS input columns
 * @param[in]  lhs_offset         Offset to be applied to the LHS input value
 * @param[in]  dst_offset         Offset to be applied the output result
 * @param[in]  activation_min     Minimum value to clamp down the output. Range : int8
 * @param[in]  activation_max     Maximum value to clamp up the output. Range : int8
 * @param[in]  rhs_cols_offset    Offset between input columns. Used to handle non-unity strides
 *                                Expected value : x * rhs_cols, where x >= 1
 *
 * @return     The function returns <code>NN_SUCCESS</code>
 *
 */
adi_sharcfx_nn_status adi_sharcfx_nn_mat_mult_nt_t_s8(const int8_t *lhs,
                                            const int8_t *rhs,
                                            const int32_t *bias,
                                            int8_t *dst,
                                            const int32_t *dst_multipliers,
                                            const int32_t *dst_shifts,
                                            const int32_t lhs_rows,
                                            const int32_t rhs_rows,
                                            const int32_t rhs_cols,
                                            const int32_t lhs_offset,
                                            const int32_t dst_offset,
                                            const int32_t activation_min,
                                            const int32_t activation_max,
                                            const int32_t rhs_cols_offset);

/*
 * @brief s8 Vector by Matrix (transposed) multiplication
 *
 * @param[in]      lhs             Input left-hand side vector
 * @param[in]      rhs             Input right-hand side matrix (transposed)
 * @param[in]      bias            Input bias
 * @param[out]     dst             Output vector
 * @param[in]      lhs_offset      Offset to be added to the input values of the left-hand side vector.
 *                                 Range: -127 to 128
 * @param[in]      dst_offset      Offset to be added to the output values. Range: -127 to 128
 * @param[in]      dst_multiplier  Output multiplier
 * @param[in]      dst_shift       Output shift
 * @param[in]      rhs_cols        Number of columns in the right-hand side input matrix
 * @param[in]      rhs_rows        Number of rows in the right-hand side input matrix
 * @param[in]      activation_min  Minimum value to clamp the output to. Range: int8
 * @param[in]      activation_max  Maximum value to clamp the output to. Range: int8
 * @param[in]      address_offset  Memory position offset for dst. First output is stored at 'dst', the
 *                                 second at 'dst + address_offset' and so on. Default value is typically 1.
 *
 * @return         The function returns <code>NN_SUCCESS</code>
 *
 */
adi_sharcfx_nn_status adi_sharcfx_nn_vec_mat_mult_t_s8(const int8_t *lhs,
                                             const int8_t *rhs,
                                             const int32_t *bias,
                                             int8_t *dst,
                                             const int32_t lhs_offset,
                                             const int32_t dst_offset,
                                             const int32_t dst_multiplier,
                                             const int32_t dst_shift,
                                             const int32_t rhs_cols,
                                             const int32_t rhs_rows,
                                             const int32_t activation_min,
                                             const int32_t activation_max,
                                             const int32_t address_offset);



/**
 * @brief           memset optimized for MVE
 * @param[in, out]  dst         Destination pointer
 * @param[in]       val         Value to set
 * @param[in]       block_size  Number of bytes to copy.
 *
 */
__STATIC_FORCEINLINE void adi_sharcfx_memset_s8(int8_t *dst, const int8_t val, uint32_t block_size)
{

    memset(dst, val, block_size);

}


/**
 * @brief Matrix-multiplication function for convolution with per-channel requantization.
 * @param[in]       input_a     pointer to operand A
 * @param[in]       input_b     pointer to operand B, always consists of 2 vectors.
 * @param[in]       output_ch   number of rows of A
 * @param[in]       out_shift  pointer to per output channel requantization shift parameter.
 * @param[in]       out_mult   pointer to per output channel requantization multiplier parameter.
 * @param[in]       out_offset      output tensor offset.
 * @param[in]       activation_min   minimum value to clamp the output to. Range : int8
 * @param[in]       activation_max   maximum value to clamp the output to. Range : int8
 * @param[in]       num_col_a   number of columns of A
 * @param[in]       output_bias per output channel bias. Range : int32
 * @param[in,out]   out_0       pointer to output
 * @return     The function returns one of the two
 *              1. The incremented output pointer for a successful operation or
 *              2. NULL if implementation is not available.
 *
 * @details   This function does the matrix multiplication of weight matrix for all output channels
 *            with 2 columns from im2col and produces two elements/output_channel. The outputs are
 *            clamped in the range provided by activation min and max.
 *            Supported framework: TensorFlow Lite micro.
 */
int8_t *adi_sharcfx_nn_mat_mult_kernel_s8_s16(const int8_t *input_a,
                                      const int16_t *input_b,
                                      const uint16_t output_ch,
                                      const int32_t *out_shift,
                                      const int32_t *out_mult,
                                      const int32_t out_offset,
                                      const int16_t activation_min,
                                      const int16_t activation_max,
                                      const uint16_t num_col_a,
                                      const int32_t *const output_bias,
                                      int8_t *out_0);

/**
 * @brief Common softmax function for s8 input and s8 or s16 output
 * @param[in]  input          Pointer to the input tensor
 * @param[in]  num_rows       Number of rows in the input tensor
 * @param[in]  row_size       Number of elements in each input row
 * @param[in]  mult           Input quantization multiplier
 * @param[in]  shift          Input quantization shift within the range [0, 31]
 * @param[in]  diff_min       Minimum difference with max in row. Used to check if
 *                            the quantized exponential operation can be performed
 * @param[in]  int16_output   Indicating s8 output if 0 else s16 output
 * @param[out] output         Pointer to the output tensor
 *
 * @note Supported framework: TensorFlow Lite micro (bit-accurate)
 *
 */
void adi_sharcfx_nn_softmax_common_s8(const int8_t *input,
                              const int32_t num_rows,
                              const int32_t row_size,
                              const int32_t mult,
                              const int32_t shift,
                              const int32_t diff_min,
                              const bool int16_output,
                              void *output);

/**
 * @brief macro for adding rounding offset
 */
#ifndef ARM_NN_TRUNCATE
#define NN_ROUND(out_shift) ((0x1 << out_shift) >> 1)
#else
#define NN_ROUND(out_shift) 0
#endif

// Macros for shortening quantization functions' names and avoid long lines
#define MUL_SAT(a, b) adi_sharcfx_nn_doubling_high_mult((a), (b))
#define MUL_SAT_MVE(a, b) adi_sharcfx_doubling_high_mult_mve_32x4((a), (b))
#define MUL_POW2(a, b) adi_sharcfx_nn_mult_by_power_of_two((a), (b))

#define DIV_POW2(a, b) adi_sharcfx_nn_divide_by_power_of_two((a), (b))
#define DIV_POW2_MVE(a, b) adi_sharcfx_divide_by_power_of_two_mve((a), (b))

#define EXP_ON_NEG(x) adi_sharcfx_nn_exp_on_negative_values((x))
#define ONE_OVER1(x) adi_sharcfx_nn_one_over_one_plus_x_for_x_in_0_1((x))

/**
 * @brief           Saturating doubling high multiply. Result matches
 *                  NEON instruction VQRDMULH.
 * @param[in]       m1        Multiplicand. Range: {NN_Q31_MIN, NN_Q31_MAX}
 * @param[in]       m2        Multiplier. Range: {NN_Q31_MIN, NN_Q31_MAX}
 * @return          Result of multiplication.
 *
 */
__STATIC_FORCEINLINE int32_t adi_sharcfx_nn_doubling_high_mult(const int32_t m1, const int32_t m2)
{
    int32_t result = 0;
    // Rounding offset to add for a right shift of 31
    int64_t mult = 1 << 30;

    if ((m1 < 0) ^ (m2 < 0))
    {
        mult = 1 - mult;
    }
    // Gets resolved as a SMLAL instruction
    mult = mult + (int64_t)m1 * m2;

    // Utilize all of the upper 32 bits. This is the doubling step
    // as well.
    result = (int32_t)(mult / (1ll << 31));

    if ((m1 == m2) && (m1 == (int32_t)NN_Q31_MIN))
    {
        result = NN_Q31_MAX;
    }
    return result;
}

/**
 * @brief           Doubling high multiply without saturation. This is intended
 *                  for requantization where the scale is a positive integer
 *
 * @param[in]       m1        Multiplicand. Range: {NN_Q31_MIN, NN_Q31_MAX}
 * @param[in]       m2        Multiplier Range: {NN_Q31_MIN, NN_Q31_MAX}
 * @return          Result of multiplication.
 * @note            The result of this matches that of neon instruction
 *                  VQRDMULH for m1 in range {NN_Q31_MIN, NN_Q31_MAX} and m2 in
 *                  range {NN_Q31_MIN + 1, NN_Q31_MAX}. Saturation occurs when
 *                  m1 equals m2 equals NN_Q31_MIN and that is not handled by
 *                  this function.
 *
 */
__STATIC_FORCEINLINE int32_t adi_sharcfx_nn_doubling_high_mult_no_sat(const int32_t m1, const int32_t m2)
{
    int32_t result = 0;
    union adi_sharcfx_nn_long_long mult;

    // Rounding offset to add for a right shift of 31
    mult.word.low = 1 << 30;
    mult.word.high = 0;

    // Gets resolved as a SMLAL instruction
    mult.long_long = mult.long_long + (int64_t)m1 * m2;

    // Utilize all of the upper 32 bits. This is the doubling step
    // as well.
    result = (int32_t)(mult.long_long >> 31);

    return result;
}

/**
 * @brief           Rounding divide by power of two.
 * @param[in]       dividend - Dividend
 * @param[in]       exponent - Divisor = power(2, exponent)
 *                             Range: [0, 31]
 * @return          Rounded result of division. Midpoint is rounded away from zero.
 *
 */
__STATIC_FORCEINLINE int32_t adi_sharcfx_nn_divide_by_power_of_two(const int32_t dividend, const int32_t exponent)
{
    int32_t result = 0;
    const int32_t remainder_mask = (1 << exponent) - 1;
    int32_t remainder = remainder_mask & dividend;

    // Basic division
    result = dividend >> exponent;

    // Adjust 'result' for rounding (mid point away from zero)
    int32_t threshold = remainder_mask >> 1;
    if (result < 0)
    {
        threshold++;
    }
    if (remainder > threshold)
    {
        result++;
    }

    return result;
}

/**
 * @brief           Requantize a given value.
 * @param[in]       val         Value to be requantized
 * @param[in]       multiplier  multiplier. Range {NN_Q31_MIN + 1, Q32_MAX}
 * @param[in]       shift       left or right shift for 'val * multiplier'
 *
 * @return          Returns (val * multiplier)/(2 ^ shift)
 *
 */
__STATIC_FORCEINLINE int32_t adi_sharcfx_nn_requantize(const int32_t val, const int32_t multiplier, const int32_t shift)
{
#ifdef CMSIS_NN_USE_SINGLE_ROUNDING
    const int64_t total_shift = 31 - shift;
    const int64_t new_val = val * (int64_t)multiplier;

    int32_t result = new_val >> (total_shift - 1);
    result = (result + 1) >> 1;

    return result;
#else
    return adi_sharcfx_nn_divide_by_power_of_two(adi_sharcfx_nn_doubling_high_mult_no_sat(val * (1 << LEFT_SHIFT(shift)), multiplier),
                                         RIGHT_SHIFT(shift));
#endif
}

 

/**
 * @brief           memcpy optimized for MVE
 * @param[in, out]  dst         Destination pointer
 * @param[in]       src         Source pointer.
 * @param[in]       block_size  Number of bytes to copy.
 *
 */
__STATIC_FORCEINLINE void adi_sharcfx_memcpy_s8(int8_t *__RESTRICT dst, const int8_t *__RESTRICT src, uint32_t block_size)
{

    memcpy(dst, src, block_size);

}

 

// @note The following functions are used only for softmax layer, scaled bits = 5 assumed

__STATIC_FORCEINLINE int32_t adi_sharcfx_nn_exp_on_negative_values(int32_t val)
{
    int32_t mask = 0;
    int32_t shift = 24;

    const int32_t val_mod_minus_quarter = (val & ((1 << shift) - 1)) - (1 << shift);
    const int32_t remainder = val_mod_minus_quarter - val;
    const int32_t x = (val_mod_minus_quarter << 5) + (1 << 28);
    const int32_t x2 = MUL_SAT(x, x);

    int32_t result = 1895147668 +
        MUL_SAT(1895147668, x + DIV_POW2(MUL_SAT(DIV_POW2(MUL_SAT(x2, x2), 2) + MUL_SAT(x2, x), 715827883) + x2, 1));

#define SELECT_IF_NON_ZERO(x)                                                                                          \
    {                                                                                                                  \
        mask = MASK_IF_NON_ZERO(remainder & (1 << shift++));                                                           \
        result = SELECT_USING_MASK(mask, MUL_SAT(result, x), result);                                                  \
    }

    SELECT_IF_NON_ZERO(1672461947)
    SELECT_IF_NON_ZERO(1302514674)
    SELECT_IF_NON_ZERO(790015084)
    SELECT_IF_NON_ZERO(290630308)
    SELECT_IF_NON_ZERO(39332535)
    SELECT_IF_NON_ZERO(720401)
    SELECT_IF_NON_ZERO(242)

#undef SELECT_IF_NON_ZERO

    mask = MASK_IF_ZERO(val);
    return SELECT_USING_MASK(mask, NN_Q31_MAX, result);
}

__STATIC_FORCEINLINE int32_t adi_sharcfx_nn_mult_by_power_of_two(const int32_t val, const int32_t exp)
{
    const int32_t thresh = ((1 << (31 - exp)) - 1);
    int32_t result = val << exp;
    result = SELECT_USING_MASK(MASK_IF_NON_ZERO(val > thresh), NN_Q31_MAX, result);
    result = SELECT_USING_MASK(MASK_IF_NON_ZERO(val < -thresh), NN_Q31_MIN, result);
    return result;
}

__STATIC_FORCEINLINE int32_t adi_sharcfx_nn_one_over_one_plus_x_for_x_in_0_1(int32_t val)
{
    const int64_t sum = (int64_t)val + (int64_t)NN_Q31_MAX;
    const int32_t half_denominator = (int32_t)((sum + (sum >= 0 ? 1 : -1)) / 2L);
    int32_t x = 1515870810 + MUL_SAT(half_denominator, -1010580540);

    const int32_t shift = (1 << 29);
    x += MUL_POW2(MUL_SAT(x, shift - MUL_SAT(half_denominator, x)), 2);
    x += MUL_POW2(MUL_SAT(x, shift - MUL_SAT(half_denominator, x)), 2);
    x += MUL_POW2(MUL_SAT(x, shift - MUL_SAT(half_denominator, x)), 2);

    return MUL_POW2(x, 1);
}

 __STATIC_FORCEINLINE uint8_t __CLZ(uint32_t data)
    {
      if (data == 0U) { return 32U; }

      uint32_t count = 0U;
      uint32_t mask = 0x80000000U;

      while ((data & mask) == 0U)
      {
        count += 1U;
        mask = mask >> 1U;
      }
      return count;
    }

#ifdef __cplusplus
}
#endif

#endif
