/**
 * Copyright(c) 2025 Analog Devices, Inc. All Rights Reserved.
 * Copyright (C) 2022 EEMBC
 * Copyright (C) 2022 Arm Limited
 * 
 *
 * All EEMBC Benchmark Software are products of EEMBC and are provided under the
 * terms of the EEMBC Benchmark License Agreements. The EEMBC Benchmark Software
 * are proprietary intellectual properties of EEMBC and its Members and is
 * protected under all applicable laws, including all applicable copyright laws.
 *
 * If you received this EEMBC Benchmark Software without having a currently
 * effective EEMBC Benchmark License Agreement, you must discontinue use.
 */



#include <stdio.h>
#include <stdint.h>
#include <stdbool.h>

#include "ee_audiomark.h"

// There are several POSIX assumptions in this implementation.
#if defined __linux__ || __APPLE__
#include <time.h>
#elif defined _WIN32
#include <sys\timeb.h>
#elif defined __arm__
#include <RTE_Components.h>
#if defined __PERF_COUNTER__
#include "perf_counter.h"
#endif
#elif defined SHARC_FX
#include "adi_pwr_SC8xx_config.h" //For getting the core clock frequency //#define CORE_CLOCK 1000000000
#include <xtensa/xtruntime.h>
#define TIMER_INTERVAL CORE_CLOCK/10
#include <inttypes.h>
//#include <libs/NN/Include/adi_sharcfx_nnsupportfunctions.h>
#else
#error "Operating system not recognized"
#endif
#include <assert.h>
#include <processor_include.h>

#if defined(SHARC_FX) && defined(HWPROF)
#include <xtensa/xt_profiling.h>

void xt_profile_setup(void)
{
	unsigned int xt_prof_sel = XTPERF_CNT_CYCLES;
	unsigned int xt_prof_msk = XTPERF_MASK_CYCLES;
	int xt_prof_level = -2;
	unsigned int xt_prof_period = 4096;//1024;
	int xt_prof_ret = 0;

	xt_prof_ret = xt_profile_config_counter(xt_prof_sel, xt_prof_msk, xt_prof_level, xt_prof_period);
	xt_profile_enable();

	return;
}

void xt_profile_end(void)
{
	xt_profile_config_clear();
	return;
}
#endif//SHARC-FX HWPROF

uint64_t
th_microseconds(void)
{
    uint64_t usec = 0;
#if defined __linux__ || __APPLE__
    const long      NSEC_PER_SEC      = 1000 * 1000 * 1000;
    const long      TIMER_RES_DIVIDER = 1000;
    struct timespec t;
    clock_gettime(CLOCK_REALTIME, &t);
    usec = t.tv_sec * (NSEC_PER_SEC / TIMER_RES_DIVIDER)
           + t.tv_nsec / TIMER_RES_DIVIDER;
#elif defined _WIN32
    struct timeb t;
    ftime(&t);
    usec = ((uint64_t)t.time) * 1000 * 1000 + ((uint64_t)t.millitm) * 1000;
#elif defined __arm__ && defined __PERF_COUNTER__
    usec = (uint64_t)get_system_us();
#elif defined SHARC_FX
#warning "ADI SHARC-FX uses its cycle counter"
#else
#error "Operating system not recognized"
#endif
    return usec;
}


uint32_t count;
void timer_handler(void *arg) {
//increasing the count by 1 after TIMER_INTERVAL cycles are completed, and re-setting the compare.
  xthal_set_ccompare(0, xthal_get_ccompare(0) + TIMER_INTERVAL);
  count += 1;
}

bool
time_audiomark_run(uint32_t iterations, uint64_t *dt)
{
    uint64_t t0  = 0;
    uint64_t t1  = 0;
    bool     err = false;

#if defined SHARC_FX
	  count = 0;
	  xtos_set_interrupt_handler(XCHAL_TIMER0_INTERRUPT, timer_handler, NULL, NULL); // calling the interrupt handler, which tracks the cycle count
	  xtos_interrupt_enable(XCHAL_TIMER0_INTERRUPT); //enable the timer interrupt
	  xthal_set_ccompare(0, xthal_get_ccount() + TIMER_INTERVAL); //set the compare value, which calls the interrupt when cycle count = timer_interval value

#else
	  t0 = th_microseconds();
#endif
	/*
	 * Loop through 'iterations' calling the 'ee_audiomark_run()' function.
	 * If 'ee_audiomark_run()' returns -1, set 'err' flag to true and exit.
	 */
    for (uint32_t i = 0; i < iterations; ++i)
    {
        if (ee_audiomark_run())
        {
            err = true;
            break;
        }
    }

#if defined SHARC_FX
    xtos_interrupt_disable(XCHAL_TIMER0_INTERRUPT); //disable the timer interrupt value after the runtime is calculated
    *dt = (count) * (TIMER_INTERVAL/(CORE_CLOCK / 1e6)); //so number of cycles = count*timer_interval, divided by core clock then multiplied by 1e6 to get microseconds
#else
    t1  = th_microseconds();
    *dt = t1 - t0;
#endif
    return err;
}

int
main(void)
{
#if defined(SHARC_FX) && defined(HWPROF)
	xt_profile_setup();
#endif//SHARC_FX HWPROF


    bool     err        = false;
    uint32_t iterations = 1;
    uint64_t dt         = 0;

    printf("Initializing\n");

#if defined SHARC_FX
    printf("Core Clock Frequency %u Hz\n", CORE_CLOCK);
#endif

    if (ee_audiomark_initialize())
    {
        printf("Failed to initialize\n");
        return -1;
    }

    printf("Computing run speed\n");

    /*
     * Keep doubling 'iterations' until the elapsed time by 'time_audiomark_run()' function exceeds 1 sec.
     * Update the 'dt' with the elased time in each iteration.
     */
    do
    {
        iterations *= 2;
        //dt = 0;
        err = time_audiomark_run(iterations, &dt);
        if (err)
        {
            break;
        }
    } while (dt < 1e6);

    if (err)
    {
        printf("Failed to compute iteration speed\n");
        goto exit;
    }

    // Must run for 10 sec. or at least 10 iterations
    float scale = 11e6 / dt;
    iterations  = (uint32_t)((float)iterations * scale);
    iterations  = iterations < 10 ? 10 : iterations;


    printf("Measuring\n");

    err = time_audiomark_run(iterations, &dt);
    if (err)
    {
        printf("Failed main performance run\n");
        goto exit;
    }

    /**
     * The input stream is 24e3 samples at 16 kHz, which means to exactly
     * match the throughput of the stream the score would be one iteration
     * per 1.5 seconds. The score is how many times faster than the ADC
     * the pipeline runs. x 1000 to make it a bigger number.
     */
    float sec   = (float)dt / 1.0e6f;
    float score = (float)iterations / sec * 1000.f * (1 / 1.5f);

    printf("Total runtime    : %.3f seconds\n", sec);
    printf("Total iterations : %d iterations\n", iterations);
    printf("Score            : %f AudioMarks\n", score);
#if defined SHARC_FX
    printf("Score            : %f AudioMarks/MHz\n", score/(CORE_CLOCK / 1000000));
#endif
exit:
    ee_audiomark_release();

#if defined(SHARC_FX) && defined(HWPROF)
    xt_profile_end();
#endif// SHARC_FX HWPROF

    return err ? -1 : 0;
}
