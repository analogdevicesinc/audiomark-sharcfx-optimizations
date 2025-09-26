/**
 * Copyright(c) 2025 Analog Devices, Inc. All Rights Reserved.
 * Copyright (C) 2022 EEMBC
 * Copyright (C) 2022 Arm Limited
 *
 * All EEMBC Benchmark Software are products of EEMBC and are provided under the
 * terms of the EEMBC Benchmark License Agreements. The EEMBC Benchmark Software
 * are proprietary intellectual properties of EEMBC and its Members and is
 * protected under all applicable laws, including all applicable copyright laws.
 *
 * If you received this EEMBC Benchmark Software without having a currently
 * effective EEMBC Benchmark License Agreement, you must discontinue use.
 */

#ifndef __TH_TYPES_H
#define __TH_TYPES_H

#include "libs/DSP/Include/dsp/matrix_functions.h"
#include "libs/DSP/Include/dsp/statistics_functions.h"
#include "libs/DSP/Include/dsp/support_functions.h"
#include "libs/DSP/Include/dsp/transform_functions.h"

#define TH_FLOAT32_TYPE                 float
#define TH_MATRIX_INSTANCE_FLOAT32_TYPE adi_sharcfx_matrix_instance_f32
#define TH_RFFT_INSTANCE_FLOAT32_TYPE   adi_sharcfx_rfft_fast_instance_f32
#define TH_CFFT_INSTANCE_FLOAT32_TYPE   adi_sharcfx_cfft_instance_f32

#endif /* __TH_TYPES_H */
