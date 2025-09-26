/* Copyright(c) 2025 Analog Devices, Inc. All Rights Reserved. 
   Copyright (C) 2003 Epic Games (written by Jean-Marc Valin)
   Copyright (C) 2004-2006 Epic Games

   File: preprocess.c
   Preprocessor with denoising based on the algorithm by Ephraim and Malah

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions are
   met:

   1. Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.

   2. Redistributions in binary form must reproduce the above copyright
   notice, this list of conditions and the following disclaimer in the
   documentation and/or other materials provided with the distribution.

   3. The name of the author may not be used to endorse or promote products
   derived from this software without specific prior written permission.

   THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR
   IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
   OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
   DISCLAIMED. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT,
   INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
   (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
   SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
   HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
   STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
   ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
   POSSIBILITY OF SUCH DAMAGE.
*/



#ifndef PREPROCESS_OPT_SHARCFX_H
#define PREPROCESS_OPT_SHARCFX_H

#include "filterbank.h"
#include "speex_echo.h"
#pragma once
/** Speex pre-processor state. */
struct SpeexPreprocessState__ {
   /* Basic info */
   int    frame_size;        /**< Number of samples processed each time */
   int    ps_size;           /**< Number of points in the power spectrum */
   int    sampling_rate;     /**< Sampling rate of the input/output */
   int    nbands;
   FilterBank *bank;

   /* Parameters */
   int    denoise_enabled;
   int    vad_enabled;
   int    dereverb_enabled;
   spx_word16_t  reverb_decay;
   spx_word16_t  reverb_level;
   spx_word16_t speech_prob_start;
   spx_word16_t speech_prob_continue;
   int    noise_suppress;
   int    echo_suppress;
   int    echo_suppress_active;
   SpeexEchoState *echo_state;

   spx_word16_t	speech_prob;  /**< Probability last frame was speech */

   /* DSP-related arrays */
   spx_word16_t *frame;      /**< Processing frame (2*ps_size) */
   spx_word16_t *ft;         /**< Processing frame in freq domain (2*ps_size) */
   spx_word32_t *ps;         /**< Current power spectrum */
   spx_word16_t *gain2;      /**< Adjusted gains */
   spx_word16_t *gain_floor; /**< Minimum gain allowed */
   spx_word16_t *window;     /**< Analysis/Synthesis window */
   spx_word32_t *noise;      /**< Noise estimate */
   spx_word32_t *reverb_estimate; /**< Estimate of reverb energy */
   spx_word32_t *old_ps;     /**< Power spectrum for last frame */
   spx_word16_t *gain;       /**< Ephraim Malah gain */
   spx_word16_t *prior;      /**< A-priori SNR */
   spx_word16_t *post;       /**< A-posteriori SNR */

   spx_word32_t *S;          /**< Smoothed power spectrum */
   spx_word32_t *Smin;       /**< See Cohen paper */
   spx_word32_t *Stmp;       /**< See Cohen paper */
   int *update_prob;         /**< Probability of speech presence for noise update */

   spx_word16_t *zeta;       /**< Smoothed a priori SNR */
   spx_word32_t *echo_noise;
   spx_word32_t *residual_echo;

   /* Misc */
   spx_word16_t *inbuf;      /**< Input buffer (overlapped analysis) */
   spx_word16_t *outbuf;     /**< Output buffer (for overlap and add) */

   /* AGC stuff, only for floating point for now */
#ifndef FIXED_POINT
   int    agc_enabled;
   float  agc_level;
   float  loudness_accum;
   float *loudness_weight;   /**< Perceptual loudness curve */
   float  loudness;          /**< Loudness estimate */
   float  agc_gain;          /**< Current AGC gain */
   float  max_gain;          /**< Maximum gain allowed */
   float  max_increase_step; /**< Maximum increase in gain from one frame to another */
   float  max_decrease_step; /**< Maximum decrease in gain from one frame to another */
   float  prev_loudness;     /**< Loudness of previous frame */
   float  init_max;          /**< Current gain limit during initialisation */
#endif
   int    nb_adapt;          /**< Number of frames used for adaptation so far */
   int    was_speech;
   int    min_count;         /**< Number of frames processed so far */
   void  *fft_lookup;        /**< Lookup table for the FFT */
#ifdef FIXED_POINT
   int    frame_shift;
#endif
};

/** State of the preprocessor (one per channel). Should never be accessed directly. */
typedef struct SpeexPreprocessState__ SpeexPreprocessState2;

#endif
