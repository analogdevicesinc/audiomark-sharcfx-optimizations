/*********************************************************************************
Copyright(c) 2023 Analog Devices, Inc. All Rights Reserved.

This software is proprietary. By using this software you agree
to the terms of the associated Analog Devices License Agreement.
*********************************************************************************/

/*!
* @file      adi_pwr_SC8xx_config.c
*
* @brief     power Service configuration file
*
* @details
*            power Service configuration file
*/

#include <sys/platform.h>
#include <stdint.h>
#include <stdlib.h>
#include <services/pwr/adi_pwr.h>
#include "adi_pwr_SC8xx_config.h"
#include <time.h>
//#include "CGUPowerServiceConfig.h" //some part of code from  : C:\Analog Devices\EV-SC8xx_EZ-KIT-Rel1.0.0\EV-SC8xx_EZ-KIT\Examples\services\pwr\CGUPowerServiceConfig.h


#ifdef _MISRA_RULES
#pragma diag(push)
#pragma diag(suppress:misra_rule_14_7:"Allow functions to have multiple exits for better readability and optimized code")
#endif

//#define ADI_PWR_DEBUG_CLOCK_PRINTF

/**
 * @brief    Initializes clocks, including CGU and CDU modules.
 *
 * @return   Status
 *           - 0: Successful in all the initializations.
 *           - 1: Error.

 */
uint32_t adi_pwr_cfg0_init()
{
    uint32_t status = 0u; /*Return zero if there are no errors*/

    /* Structure pointer for CGU0 and CGU1 parameters*/
    ADI_PWR_CGU_PARAM_LIST pADI_CGU_Param_List;

    /* Structure pointer for CDU parameters*/
    ADI_PWR_CDU_PARAM_LIST pADI_CDU_Param_List;

    /* CDU Configuration*/
    pADI_CDU_Param_List.cdu_settings[0].cfg_SEL                     =       (ADI_PWR_CDU_CLKIN)CFG0_BIT_CDU0_CFG0_SEL_VALUE;
    pADI_CDU_Param_List.cdu_settings[0].cfg_EN                      =       true;

    pADI_CDU_Param_List.cdu_settings[2].cfg_SEL                     =       (ADI_PWR_CDU_CLKIN)CFG0_BIT_CDU0_CFG2_SEL_VALUE;
    pADI_CDU_Param_List.cdu_settings[2].cfg_EN                      =       true;

    pADI_CDU_Param_List.cdu_settings[3].cfg_SEL                     =       (ADI_PWR_CDU_CLKIN)CFG0_BIT_CDU0_CFG3_SEL_VALUE;
    pADI_CDU_Param_List.cdu_settings[3].cfg_EN                      =       true;

    pADI_CDU_Param_List.cdu_settings[4].cfg_SEL                     =       (ADI_PWR_CDU_CLKIN)CFG0_BIT_CDU0_CFG4_SEL_VALUE;
    pADI_CDU_Param_List.cdu_settings[4].cfg_EN                      =       true;

    pADI_CDU_Param_List.cdu_settings[5].cfg_SEL                     =       (ADI_PWR_CDU_CLKIN)CFG0_BIT_CDU0_CFG5_SEL_VALUE;
    pADI_CDU_Param_List.cdu_settings[5].cfg_EN                      =       true;

    pADI_CDU_Param_List.cdu_settings[6].cfg_SEL                     =       (ADI_PWR_CDU_CLKIN)CFG0_BIT_CDU0_CFG6_SEL_VALUE;
    pADI_CDU_Param_List.cdu_settings[6].cfg_EN                      =       true;

    pADI_CDU_Param_List.cdu_settings[7].cfg_SEL                     =       (ADI_PWR_CDU_CLKIN)CFG0_BIT_CDU0_CFG7_SEL_VALUE;
    pADI_CDU_Param_List.cdu_settings[7].cfg_EN                      =       true;

    pADI_CDU_Param_List.cdu_settings[8].cfg_SEL                     =       (ADI_PWR_CDU_CLKIN)CFG0_BIT_CDU0_CFG8_SEL_VALUE;
    pADI_CDU_Param_List.cdu_settings[8].cfg_EN                      =       true;

    pADI_CDU_Param_List.cdu_settings[9].cfg_SEL                     =       (ADI_PWR_CDU_CLKIN)CFG0_BIT_CDU0_CFG9_SEL_VALUE;
    pADI_CDU_Param_List.cdu_settings[9].cfg_EN                      =       true;

    pADI_CDU_Param_List.cdu_settings[10].cfg_SEL                    =       (ADI_PWR_CDU_CLKIN)CFG0_BIT_CDU0_CFG10_SEL_VALUE;
    pADI_CDU_Param_List.cdu_settings[10].cfg_EN                     =       true;


    pADI_CDU_Param_List.cdu_settings[12].cfg_SEL                    =       (ADI_PWR_CDU_CLKIN)CFG0_BIT_CDU0_CFG12_SEL_VALUE;
    pADI_CDU_Param_List.cdu_settings[12].cfg_EN                     =       true;

    pADI_CDU_Param_List.cdu_settings[13].cfg_SEL                    =       (ADI_PWR_CDU_CLKIN)CFG0_BIT_CDU0_CFG13_SEL_VALUE;
    pADI_CDU_Param_List.cdu_settings[13].cfg_EN                     =       true;


    /* CGU0 Configuration*/
    pADI_CGU_Param_List.cgu0_settings.clocksettings.ctl_MSEL        =       (uint32_t)CFG0_BIT_CGU0_CTL_MSEL;
    pADI_CGU_Param_List.cgu0_settings.clocksettings.ctl_DF          =       (uint32_t)CFG0_BIT_CGU0_CTL_DF;
    pADI_CGU_Param_List.cgu0_settings.clocksettings.div_CSEL        =       (uint32_t)CFG0_BIT_CGU0_DIV_CSEL;
    pADI_CGU_Param_List.cgu0_settings.clocksettings.div_SYSSEL      =       (uint32_t)CFG0_BIT_CGU0_DIV_SYSSEL;
    pADI_CGU_Param_List.cgu0_settings.clocksettings.div_S0SEL       =       (uint32_t)CFG0_BIT_CGU0_DIV_S0SEL;
    pADI_CGU_Param_List.cgu0_settings.clocksettings.div_S1SEL       =       (uint32_t)CFG0_BIT_CGU0_DIV_S1SEL;
    pADI_CGU_Param_List.cgu0_settings.clocksettings.divex_S1SELEX   =       (uint32_t)CFG0_BIT_CGU0_DIV_S1SELEX;
    pADI_CGU_Param_List.cgu0_settings.clocksettings.div_DSEL        =       (uint32_t)CFG0_BIT_CGU0_DIV_DSEL;
    pADI_CGU_Param_List.cgu0_settings.clocksettings.div_OSEL        =       (uint32_t)CFG0_BIT_CGU0_DIV_OSEL;
    pADI_CGU_Param_List.cgu0_settings.clkin                         =       (uint32_t)CFG0_BIT_CGU0_CLKIN;
    pADI_CGU_Param_List.cgu0_settings.enable_SCLK1ExDiv             =       true;

    /* CGU1 Configuration*/
    pADI_CGU_Param_List.cgu1_settings.clocksettings.ctl_MSEL        =       (uint32_t)CFG0_BIT_CGU1_CTL_MSEL;
    pADI_CGU_Param_List.cgu1_settings.clocksettings.ctl_DF          =       (uint32_t)CFG0_BIT_CGU1_CTL_DF;
    pADI_CGU_Param_List.cgu1_settings.clocksettings.div_CSEL        =       (uint32_t)CFG0_BIT_CGU1_DIV_CSEL;
    pADI_CGU_Param_List.cgu1_settings.clocksettings.div_SYSSEL      =       (uint32_t)CFG0_BIT_CGU1_DIV_SYSSEL;
    pADI_CGU_Param_List.cgu1_settings.clocksettings.div_S0SEL       =       (uint32_t)CFG0_BIT_CGU1_DIV_S0SEL;
    pADI_CGU_Param_List.cgu1_settings.clocksettings.div_S1SEL       =       (uint32_t)CFG0_BIT_CGU1_DIV_S1SEL;
    pADI_CGU_Param_List.cgu1_settings.clocksettings.divex_S0SELEX   =       (uint32_t)CFG0_BIT_CGU1_DIV_S0SELEX;
    pADI_CGU_Param_List.cgu1_settings.clocksettings.divex_S1SELEX   =       (uint32_t)CFG0_BIT_CGU1_DIV_S1SELEX;
    pADI_CGU_Param_List.cgu1_settings.clocksettings.div_DSEL        =       (uint32_t)CFG0_BIT_CGU1_DIV_DSEL;
    pADI_CGU_Param_List.cgu1_settings.clocksettings.div_OSEL        =       (uint32_t)CFG0_BIT_CGU1_DIV_OSEL;
    pADI_CGU_Param_List.cgu1_settings.clkin                         =       (uint32_t)CFG0_BIT_CGU1_CLKIN;
    pADI_CGU_Param_List.cgu1_settings.enable_SCLK0ExDiv             =       true;
    pADI_CGU_Param_List.cgu1_settings.enable_SCLK1ExDiv             =       true;

    pADI_CGU_Param_List.cgu1_settings.cgu1_clkinsel                 =       (ADI_PWR_CDU_CLK_SELECT)ADI_PWR_CDU_CLK_SELECT_CLKIN0;


    /* Initialize all the clocks */
    if (adi_pwr_ClockInit(&pADI_CGU_Param_List, &pADI_CDU_Param_List) != ADI_PWR_SUCCESS)
    {
       /* Return non-zero */
       status = 1u;
    }

    return status;
   }


/* Power service return code*/
ADI_PWR_RESULT result;


typedef struct
{
  uint32_t  sysclk; /**< SYSCLK Frequency */
  uint32_t  sclk0;  /**< SCLK0 Frequency */
  uint32_t  sclk1;  /**< SCLK1 Frequency */
  uint32_t  dclk;   /**< DCLK Frequency */
  uint32_t  oclk;   /**< OCLK Frequency */
  uint32_t  cclk;   /**< CCLK Frequency */
  uint32_t  pllclk; /**< PLLCLK Frequency */
} ADI_PWR_CLK_STRUCT_VAL;


#define CGUParams 2

ADI_PWR_CLK_STRUCT_VAL fclkstructval[CGUParams];
ADI_PWR_CLK_STRUCT     fclkstruct0 =
{
  &fclkstructval[0].sysclk,
  &fclkstructval[0].sclk0,
  &fclkstructval[0].sclk1,
  &fclkstructval[0].dclk,
  &fclkstructval[0].oclk,
  &fclkstructval[0].cclk,
  &fclkstructval[0].pllclk
};

ADI_PWR_CLK_STRUCT     fclkstruct1 =
{
  &fclkstructval[1].sysclk,
  &fclkstructval[1].sclk0,
  &fclkstructval[1].sclk1,
  &fclkstructval[1].dclk,
  &fclkstructval[1].oclk,
  &fclkstructval[1].cclk,
  &fclkstructval[1].pllclk
};

void adi_initPwrClock(void)
{

	  /* Initialize the power service with configuration provided*/
	  if(adi_pwr_cfg0_init() != 0)
	  {
	     printf("Failed to configure PLL \n");
	  }

	  /* Read back the CGU0 frequencies and print to console */
	  result = adi_pwr_GetFreq(ADI_PWR_CGU0, fclkstruct0);
	  if(result != ADI_PWR_SUCCESS)
	  {
	     printf("Failed to read frequencies \n");
	  }
	  /* Read back the CGU1 frequencies and print to console */
	  result = adi_pwr_GetFreq(ADI_PWR_CGU1,fclkstruct1);
	  if(result != ADI_PWR_SUCCESS)
	  {
	     printf("Failed to read frequencies \n");
	  }

	  printf("All done\n");
}

uint32_t adi_clk_frequency(void)
{
	  /* Initialize the power service with configuration provided*/
	  if(adi_pwr_cfg0_init() != 0)
	  {
	     printf("Failed to configure PLL \n");
	  }

	  /* Read back the CGU0 frequencies and print to console */
	  result = adi_pwr_GetFreq(ADI_PWR_CGU0, fclkstruct0);
	  uint32_t core_clock = *fclkstruct0.fcclk;
	  if(result != ADI_PWR_SUCCESS)
	  {
	     printf("Failed to read frequencies \n");
	  }
	  return core_clock;
}

