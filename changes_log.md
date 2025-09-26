### Overview
This file documents all significant changes made to the original source code, in accordance with the requirements of the Apache 2.0 License. The modifications primarily include performance and architecture-specific optimizations for the SHARC-FX processor by Analog Devices.

### Summary of changes
| Date             | Release Version | Component / File(s) Modified           | Description of Change                                                                 |
|------------------|-----------------|----------------------------------------|----------------------------------------------------------------------------------------|
| 2025-04-22 | V.1.0           | audiomark/main.c                 | Added SHARC-FX compatible cycle counting and profiling to compute AudioMark score for sharc-fx processor                          |
| 2025-04-22 | V.1.0           | audiomark/lib/speexdsp/libspeexdsp/filterbank_opt.c              | Added conditional logic to use SHARC-FX optimized functions from filterbank_opt_sharcfx.c         |
| 2025-04-22 | V.1.0           | audiomark/lib/speexdsp/libspeexdsp/mdf_opt.c                     | Added conditional logic to use SHARC-FX optimized functions from mdf_opt_sharcfx.c                |
| 2025-04-22 | V.1.0           | audiomark/lib/speexdsp/libspeexdsp/preprocess_opt.c              | Added conditional logic to use SHARC-FX optimized functions from preprocess_opt_sharcfx.c         |
| 2025-04-22 | V.1.0           | audiomark/lib/speexdsp/libspeexdsp/filterbank_opt_sharcfx.c      | Created new file with SHARC-FX optimized filterbank functions using override macros               |
| 2025-04-22 | V.1.0           | audiomark/lib/speexdsp/libspeexdsp/mdf_opt_sharcfx.c            | Created new file with SHARC-FX optimized MDF functions using override macros                    |
| 2025-04-22 | V.1.0           | audiomark/lib/speexdsp/libspeexdsp/preprocess_opt_sharcfx.c    | Created new file with SHARC-FX optimized preprocess functions using override macros               |
| 2025-04-22 | V.1.0           | audiomark/hardware/adi/adsp-21835-kit/system/ | Added board-specific system folder for ADSP-21835 kit containing startup and initialization files |
| 2025-04-22 | V.1.0           | audiomark/ports/adi_sharcfx/libs/NN/Source/, DSP/Source/            | Added SHARC-FX optimized NN and DSP folders inspired by CMSIS-NN and CMSIS-DSP; code rewritten for SHARC-FX architecture with appropriate copyright attribution |
| 2025-04-22 | V.1.0            | audiomark/ports/adi_sharcfx/libs/NN/Include/, DSP/Include/, DSP/Include/dsp/ | Modified header files to rename functions with `adi_sharcfx_` prefix for SHARC-FX compatibility and namespace clarity |
| 2025-04-22 | V.1.0          |  audiomark/ports/adi_sharcfx/th_api.c, adi_sharcfx/th_api.h                    | Modified function names to use SHARC-FX specific implementations |
| 2025-04-22 | V.1.0           | audiomark/ports/adi_sharcfx/libs/DSP/Source/BasicMathFunctions/adi_sharcfx_abs_f32.c | Modified to use ADI SHARC-FX library function `svecabsf` to compute an absolute value of 32-bit float in a vector. |
| 2025-04-22 | V.1.0           | audiomark/ports/adi_sharcfx/libs/DSP/Source/BasicMathFunctions/adi_sharcfx_add_f32.c | Modified to use ADI SHARC-FX library function `vecvaddf` to compute addition of two 32-bit float vectors. |
| 2025-04-22 | V.1.0           | audiomark/ports/adi_sharcfx/libs/DSP/Source/BasicMathFunctions/adi_sharcfx_dot_prod_f32.c | Modified to use ADI SHARC-FX library function `vecdotf` to compute vector dot product. |
| 2025-04-22 | V.1.0           | audiomark/ports/adi_sharcfx/libs/DSP/Source/BasicMathFunctions/adi_sharcfx_mult_f32.c | Modified to use ADI SHARC-FX library function `vecvmltf` to compute multiplication of two 32-bit float vectors. |
| 2025-04-22 | V.1.0           | audiomark/ports/adi_sharcfx/libs/DSP/Source/BasicMathFunctions/adi_sharcfx_offset_f32.c | Modified to use ADI SHARC-FX library function `vecsaddf` to add offset value to the vector. |
| 2025-04-22 | V.1.0           | audiomark/ports/adi_sharcfx/libs/DSP/Source/BasicMathFunctions/adi_sharcfx_scale_f32.c | Modified to use ADI SHARC-FX library function `vecsmltf` to multiply vector with a scalar value. |
| 2025-04-22 | V.1.0           | audiomark/ports/adi_sharcfx/libs/DSP/Source/BasicMathFunctions/adi_sharcfx_sub_f32.c | Modified to use ADI SHARC-FX library function `vecvsubf` for vector subtraction. |
| 2025-04-22 | V.1.0           | audiomark/ports/adi_sharcfx/libs/DSP/Source/CommonTables/adi_sharcfx_common_Tables.c | Modified to add SHARC-FX FFT library compatible twiddle tables and FFT tables for different FFT sizes. |
| 2025-04-22 | V.1.0           | audiomark/ports/adi_sharcfx/libs/DSP/Source/CommonTables/adi_sharcfx_const_structs.c | Added the definition of twiddle coefficients for SHARC-FX. |
| 2025-04-22 | V.1.0           | audiomark/ports/adi_sharcfx/libs/DSP/Source/ComplexMathFunctions/adi_sharcfx_cmplx_conj_f32.c | Modified to use ADI SHARC-FX library function `cvecconjf` to calculate floating point complex conjugate. |
| 2025-04-22 | V.1.0           | audiomark/ports/adi_sharcfx/libs/DSP/Source/ComplexMathFunctions/adi_sharcfx_cmplx_dot_prod_f32.c | Wrote intrinsics optimized code to compute complex dot product of two vectors. |
| 2025-04-22 | V.1.0           | audiomark/ports/adi_sharcfx/libs/DSP/Source/ComplexMathFunctions/adi_sharcfx_cmplx_mag_f32.c | Modified to use ADI SHARC-FX library function `cvecabsf` to calculate absolute magnitude of complex vector. |
| 2025-04-22 | V.1.0           | audiomark/ports/adi_sharcfx/libs/DSP/Source/ComplexMathFunctions/adi_sharcfx_cmplx_mag_squared_f32.c | Modified to use ADI SHARC-FX library function `cvecabsf` to get floating point complex magnitude squared. |
| 2025-04-22 | V.1.0           | audiomark/ports/adi_sharcfx/libs/DSP/Source/ComplexMathFunctions/adi_sharcfx_cmplx_mult_cmplx_f32.c | Modified to use ADI SHARC-FX library function `cvecvmltf` for complex-to-complex multiplication. |
| 2025-04-22 | V.1.0           | audiomark/ports/adi_sharcfx/libs/DSP/Source/FastMathFunctions/adi_sharcfx_vexp_f32.c | Modified to use ADI SHARC-FX library function `vecexpf` for calculating exponential of vector. |
| 2025-04-22 | V.1.0           | audiomark/ports/adi_sharcfx/libs/DSP/Source/FastMathFunctions/adi_sharcfx_vlog_f32.c | Modified to use ADI SHARC-FX function `veclogf ` for calculating logarithm of vector. |
| 2025-04-22 | V.1.0           | audiomark/ports/adi_sharcfx/libs/DSP/Source/MatrixFunctions/adi_sharcfx_mat_vec_mult_f32.c | Wrote intrinsics optimized code to compute floating point matrix and vector multiplication. |
| 2025-04-22 | V.1.0           | audiomark/ports/adi_sharcfx/libs/DSP/Source/StatisticsFunctions/adi_sharcfx_absmax_f32.c | Updated C code optimization for calculating maximum absolute value. |
| 2025-04-22 | V.1.0           | audiomark/ports/adi_sharcfx/libs/DSP/Source/SupportFunctions/adi_sharcfx_copy_f32.c | Updated C code optimization for copying elements from one vector to another. |
| 2025-04-22 | V.1.0           | audiomark/ports/adi_sharcfx/libs/DSP/Source/SupportFunctions/adi_sharcfx_fill_f32.c | Wrote intrinsics optimized code for filling constant value into a floating point vector. |
| 2025-04-22 | V.1.0           | audiomark/ports/adi_sharcfx/libs/DSP/Source/SupportFunctions/adi_sharcfx_float_to_q15.c | Wrote intrinsics optimized code for converting elements of float point to Q15 vector. |
| 2025-04-22 | V.1.0           | audiomark/ports/adi_sharcfx/libs/DSP/Source/SupportFunctions/adi_sharcfx_q15_to_float.c | Wrote intrinsics optimized code for converting elements of Q15 to float point vector. |
| 2025-04-22 | V.1.0           | audiomark/ports/adi_sharcfx/libs/DSP/Source/TransformFunctions/adi_sharcfx__cfft_f32.c | Modified to use ADI SHARC-FX library function `cfft` for optimized FFT code. |
| 2025-04-22 | V.1.0           | audiomark/ports/adi_sharcfx/libs/DSP/Source/TransformFunctions/adi_sharcfx_rfft_fast_f32.c | Modified to use ADI SHARC-FX library function `rfftf` for optimized real FFT code. |
| 2025-04-22 | V.1.0           | audiomark/ports/adi_sharcfx/libs/NN/Source/ConvolutionFunctions/adi_sharcfx_convolve_s8.c | Wrote intrinsic optimized code to compute 8-bit integer convolution. |
| 2025-04-22 | V.1.0           | audiomark/ports/adi_sharcfx/libs/NN/Source/ConvolutionFunctions/adi_sharcfx_depthwise_conv_3x3_s8.c | Wrote intrinsic optimized code for depthwise 3x3 convolution. |
| 2025-04-22 | V.1.0           | audiomark/ports/adi_sharcfx/libs/NN/Source/ConvolutionFunctions/adi_sharcfx_nn_mat_mult_kernel_s8_s16.c | Wrote optimized intrinsic code for 8-bit and 16-bit integer matrix multiplication. |
| 2025-04-22 | V.1.0           | audiomark/ports/adi_sharcfx/libs/NN/Source/NNSupportFunctions/adi_sharcfx_nn_mat_mult_nt_t_s8.c | Wrote intrinsic optimized code for matrix multiplication with RHS matrix transposed. |
| 2025-04-22 | V.1.0           | audiomark/ports/adi_sharcfx/libs/NN/Source/NNSupportFunctions/adi_sharcfx_nn_vec_mat_mult_t_s8.c | Wrote intrinsic optimized code for matrix and vector multiplication. |
| 2025-04-22 | V.1.0           | audiomark/ports/adi_sharcfx/libs/NN/Source/NNSupportFunctions/adi_sharcfx_q7_to_q15_with_offset.c | Wrote intrinsic optimized code for Q7 to Q15 vector conversion with offset. |
| 2025-04-22 | V.1.0           | audiomark/ports/adi_sharcfx/libs/NN/Source/PoolingFunctions/adi_sharcfx_avgpool_s8.c | Wrote intrinsic optimized code for average pooling function. |
| 2025-04-22 | V.1.0           | audiomark/ports/adi_sharcfx/libs/NN/Source/PoolingFunctions/adi_sharcfx_softmax_s8.c | Wrote intrinsic optimized code for softmax function. |

### Notes
- Optimizations are tailored for SHARC-FX architecture and may not be portable to other platforms.
- Original authorship and license headers have been preserved in all modified files.