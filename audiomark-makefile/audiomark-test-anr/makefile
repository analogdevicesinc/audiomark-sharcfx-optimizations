################################################################################
# Project makefile
################################################################################

# Build tool settings
RM := rm
SHARCFX_CC := "C:\analog\cces\3.0.2\ccfx"
SHARCFX_CCX := "C:\analog\cces\3.0.2\ccfx++"
# Project settings
PROJECT = audiomark
PROC = ADSP-21835
SI_REVISION = 0.0

# Set optimizer flags
SHARCFX_OPTIMIZE ?= -O3

# retrieve the command-line release flag (default no release)
ifeq ($(RELEASE), yes)
	BUILD_RELEASE =
else
	BUILD_RELEASE = -g
endif

################################################################################
# SHARCFX section
################################################################################

# All of the SHARCFX sources participating in the build are defined here
SRC_PREFIX = ..\..
SHARCFX_SRC_DIRS = \
	audiomark/tests \
	audiomark/tests/data \
	audiomark/hardware/adi/adsp-21835-ezkit/system \
	audiomark/hardware/adi/adsp-21835-ezkit/system/startup_lsp/lsp \


# Include directories
SHARCFX_INCLUDE_DIRS = \
	-I"../../audiomark/lib/speexdsp/include" \
	-I"../../audiomark/lib/speexdsp/include/speex" \
	-I"../../audiomark/lib/speexdsp/libspeexdsp" \
	-I"../../audiomark/src" \
	-I"../../audiomark/src/ee_data" \
	-I"../../audiomark/hardware/adi/adsp-21835-ezkit/system" \
	-I"../../audiomark/hardware/adi/adsp-21835-ezkit/system/startup_lsp/lsp" \
	-I"../../audiomark" \
	-I"../../audiomark/ports/adi_sharcfx" \
	-I"../../audiomark/ports/adi_sharcfx/libs/Core/Include" \
	-I"../../audiomark/ports/adi_sharcfx/libs/NN/Include" \
	-I"../../audiomark/ports/adi_sharcfx/libs/DSP/Include" \
	-I"../../audiomark/ports/adi_sharcfx/libs/DSP/Include/dsp" \
	-I"C:\analog\cces\3.0.0\Xtensa\SHARC-FX\include\sharc" \
	-I"C:\analog\cces\3.0.0\Xtensa\XtensaSrc\libdsp\include" \
	-I"C:\analog\cces\3.0.0\Xtensa\XtensaSrc\libdsp\include_private" \

EXCLUDE_FILE = \
	$(SRC_PREFIX)/audiomark/hardware/adi/adsp-21835-ezkit/system/adi_pwr_SC8xx_config.c \
	$(SRC_PREFIX)/audiomark/tests/test_abf_f32.c \
	$(SRC_PREFIX)/audiomark/tests/test_aec_f32.c \
	$(SRC_PREFIX)/audiomark/tests/test_kws.c \
	$(SRC_PREFIX)/audiomark/tests/test_mfcc_f32.c \
	$(SRC_PREFIX)/audiomark/tests/data/abf_f32_expected.c \
	$(SRC_PREFIX)/audiomark/tests/data/abf_f32_input_ch1.c \
	$(SRC_PREFIX)/audiomark/tests/data/abf_f32_input_ch2.c \
	$(SRC_PREFIX)/audiomark/tests/data/aec_f32_expected.c \
	$(SRC_PREFIX)/audiomark/tests/data/aec_f32_input_echo.c \
	$(SRC_PREFIX)/audiomark/tests/data/aec_f32_input_source.c \
	$(SRC_PREFIX)/audiomark/tests/data/kws_expected.c \
	$(SRC_PREFIX)/audiomark/tests/data/kws_input.c \
	$(SRC_PREFIX)/audiomark/tests/data/mfcc_f32_all.c

# SHARCFX_C_SRC := $(foreach srcdir, $(SHARCFX_SRC_DIRS), $(wildcard $(SRC_PREFIX)/$(srcdir)/*.c))
SHARCFX_C_SRC := $(filter-out $(EXCLUDE_FILE),$(foreach srcdir, $(SHARCFX_SRC_DIRS), $(wildcard $(SRC_PREFIX)/$(srcdir)/*.c)))
SHARCFX_C_SRC += $(SRC_PREFIX)/audiomark/tests/test_anr_f32.c
SHARCFX_C_OBJ := $(subst $(SRC_PREFIX)/,,${SHARCFX_C_SRC:%.c=%.o})

# retrieve dependencies
SHARCFX_DEPS = $(foreach depdir, $(SHARCFX_SRC_DIRS), $(wildcard $(depdir)/*.d))

# define the compiler flags
SHARCFX_CFLAGS = $(SHARCFX_OPTIMIZE) $(BUILD_RELEASE) $(SHARCFX_INCLUDE_DIRS)
SHARCFX_CFLAGS += -proc $(PROC) -si-revision $(SI_REVISION) -c -LNO:simd -ffunction-sections -fdata-sections -mlongcalls -g
SHARCFX_CFLAGS += -DSHARC_FX -DNDEBUG -DFLOATING_POINT -DD__RESTRICT= -DOS_SUPPORT_CUSTOM -DCORE1 -DEXPORT=

SHARCFX_OBJS = $(SHARCFX_ASM_OBJ) $(SHARCFX_C_OBJ)
SHARCFX_LIB = ../audiomark-cces-library/audiomark-cces-library.a
SHARCFX_EXE = $(PROJECT)-test-anr.dxe

# compile SHARCFX 'C' files
audiomark/%.o: $(SRC_PREFIX)/audiomark/%.c
	$(SHARCFX_CC) $(SHARCFX_CFLAGS) -I"$(SRC_PREFIX)/$(dir $@)" -MMD -MP -MF "$(basename $@).d" -o "$@" -c "$<"

.DEFAULT_GOAL = all
all: builddirs $(SHARCFX_EXE)


$(SHARCFX_EXE): $(SHARCFX_OBJS) $(SHARCFX_LIB) 
	$(SHARCFX_CCX) -proc ADSP-21835 -si-revision 0.0 -Wl,--gc-sections -Wl,--defsym=__prefctl_default=0x88 -mlsp="../../audiomark/hardware/adi/adsp-21835-ezkit/system/startup_lsp/lsp" -lm -o "$@" $^

# Other Targets
clean:
	$(RM) -f *.xml $(PROJECT).* ../*.xml *.map
	$(RM) -rf $(SHARCFX_EXE) $(SHARCFX_SRC_DIRS)
	$(RM) -rf $(PROJECT)

builddirs:
	@mkdir -p $(SHARCFX_SRC_DIRS)

help:
	@echo 'usage:'
	@echo '    make [all|clean] [SHARCFX_OPTIMIZE=<-Og,-O1,etc.>] [RELEASE=<yes/no>]'
	@echo ''
	@echo 'DEBUG:'
	@echo '    make all'
	@echo ''
	@echo 'RELEASE:'
	@echo '    make all SHARCFX_OPTIMIZE=-O1 RELEASE=yes'


.PHONY: all clean help builddirs
.SECONDARY:

# pull in and check dependencies
-include $(SHARCFX_C_OBJ:.o=.d)