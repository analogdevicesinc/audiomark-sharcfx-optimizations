/**
 * Copyright (C) 2024 SPEC Embedded Group
 * Copyright (C) 2022 EEMBC
 *
 * All EEMBC Benchmark Software are products of EEMBC and are provided under the
 * terms of the EEMBC Benchmark License Agreements. The EEMBC Benchmark Software
 * are proprietary intellectual properties of EEMBC and its Members and is
 * protected under all applicable laws, including all applicable copyright laws.
 *
 * If you received this EEMBC Benchmark Software without having a currently
 * effective EEMBC Benchmark License Agreement, you must discontinue use.
 */

#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include "ee_types.h"
#include "ee_audiomark.h"

#define NBUFFERS 93
#define NINFERS  73
#define NSAMPLES 256
#define NCLASSES 12

/* Noise to signal ratio */
#define NSRM35DB 0.017783f

//#define DEBUG_EXACT_BITS
#define MAX(a,b) (((a)>(b))?(a):(b))

extern const int16_t p_input[NBUFFERS][NSAMPLES];
extern const int8_t  p_expected[NINFERS][NCLASSES];

// Used deep inside audiomark core
//char *spxGlobalHeapPtr;
//char *spxGlobalHeapEnd;

int32_t ee_kws_f32(int32_t command,
                   void  **pp_instance,
                   void   *p_data,
                   void   *p_params);

static int16_t        aec_output[256];     // 5
static int16_t        audio_fifo[13 * 64]; // 6
static int8_t         mfcc_fifo[490];      // 7
static int8_t         classes[12];         // 8
static xdais_buffer_t xdais[4];

int
main(int argc, char *argv[])
{
    int           err           = 0;
    int           new_inference = 0;
    const int8_t *p_check       = NULL;
    int           idx_check     = 0;
    uint32_t      memreq        = 0;
    uint32_t     *p_req         = &memreq;
    void         *memory        = NULL;
    void         *inst          = NULL;
    uint32_t     A              = 0;
    uint32_t     B              = 0;
    float        ratio          = 0.0f;
    int          i, j;

    int inferences = 0;

    ee_kws_f32(NODE_MEMREQ, (void **)&p_req, NULL, NULL);

    printf("KWS F32 MEMREQ = %d bytes\n", memreq);
    memory = malloc(memreq);
    if (!memory)
    {
        printf("malloc() fail\n");
        return -1;
    }
    inst = (void *)memory;
    SETUP_XDAIS(xdais[0], aec_output, 512);
    SETUP_XDAIS(xdais[1], audio_fifo, 13 * 64 * 2);
    SETUP_XDAIS(xdais[2], mfcc_fifo, 490);
    SETUP_XDAIS(xdais[3], classes, 12);

    ee_kws_f32(NODE_RESET, (void **)&inst, NULL, NULL);

    for (i = 0; i < NBUFFERS; ++i)
    {
        memcpy(aec_output, p_input[i], 512 /* 256 samples @ 2bytes@ */);
        ee_kws_f32(NODE_RUN, (void **)&inst, xdais, &new_inference);

        /* printf("inferences=%d, i=%d, idx_check=%d\n", inferences, i, idx_check); */

        /* check both classes are noises */
        A = B = -127;
        p_check = p_expected[idx_check];
        for (j = 0; j < NCLASSES; ++j)
            {  A = MAX(A, classes[j]); /* Look for max value in the calculated result */
               B = MAX(B, p_check[j]); /* Look for max value in the expected result */
            }
        if ( (A < 0)  && (B < 0)) {
          if (new_inference) {
            ++inferences;
            ++idx_check;
          }
          continue; /* Both are less than 0, considered as noise and skip */
        }
        A = 0; /* sum of abs(signals) */
        B = 0; /* sum of abs(errors) */

        if (new_inference)
        {
            ++inferences;
            p_check = p_expected[idx_check];
            ++idx_check;

            for (int j = 0; j < NCLASSES; ++j)
            {
            A += abs(128 + ((int32_t) classes[j])); /* Shift to eliminate noises */
            B += abs(((int32_t) classes[j]) - ((int32_t)p_check[j]));

#ifdef DEBUG_EXACT_BITS
                if (classes[j] != p_check[j])
                {
                    err = 1;
                    printf("buffer[%d]class[%d]: Got %d, expected %d - FAIL\n",
                           i,
                           j,
                           classes[j],
                           p_check[j]);
                }
#endif
            }
            ratio = (float)B / (float)A; /* Noise to signal ratio */
            if (ratio > NSRM35DB)
            {
                err = true;
                printf("KWS FAIL: Inference #%d exceeded -35 dB SNR\n", i);
            }

        }
    }

    if (inferences == 0)
    {
        err = 1;
        printf("KWS did not perform any inferences\n");
    }

    if (inferences != NINFERS)
    {
        err = 1;
        printf("KWS expected %d inferences but got %d\n", NINFERS, inferences);
    }

    if (err)
    {
        printf("KWS test failed\n");
        return -1;
    }

    printf("KWS test passed\n");
    return 0;
}
