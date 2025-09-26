
/* Copyright(c) 2025 Analog Devices, Inc. All Rights Reserved. */
/* Copyright (C) 2006 Jean-Marc Valin */
/**
   @file filterbank.c
   @brief Converting between psd and filterbank
 */
/*
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


#include <math.h>
#include <stdint.h>
#include "arch.h"
#include "filterbank.h"


#ifdef OVERRIDE_FB_COMPUTE_BANK32
void filterbank_compute_bank32(FilterBank * bank, spx_word32_t * ps, spx_word32_t * mel)
{
    //Convert linear power spectrum in MEL perceptual scale.
    for (int i = 0; i < bank->nb_banks; i++)
        mel[i] = 0;
    for (int i = 0; i < bank->len; i++) 
    {
        int id1 = bank->bank_left[i];
        int id2 = bank->bank_right[i];
        mel[id1] += bank->filter_left[i]*ps[i];
        mel[id2] += bank->filter_right[i]*ps[i];
    }
}
#endif

#ifdef OVERRIDE_FB_COMPUTE_PSD16
void filterbank_compute_psd16(FilterBank * bank, spx_word16_t * mel, spx_word16_t * ps)
{
    //Compute the linear power spectral density from MEL perceptual scale.
    for (int i = 0; i < bank->len; i++) 
    {
        int id1 = bank->bank_left[i];
        int id2 = bank->bank_right[i];
        spx_word32_t tmp = mel[id1] * bank->filter_left[i] + mel[id2] *bank->filter_right[i];
        ps[i] = EXTRACT16(tmp);
    }
}

#endif
