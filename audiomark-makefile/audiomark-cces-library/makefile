################################################################################
# Project makefile
################################################################################

# Build tool settings
RM := rm
SHARCFX_CC := "C:\analog\cces\3.0.2\ccfx"
SHARCFX_CCX := "C:\analog\cces\3.0.2\ccfx++"
AR = "C:\analog\cces\3.0.2\Xtensa\XtensaTools\bin\xt-ar"
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
	audiomark/lib/speexdsp/libspeexdsp \
	audiomark/ports/adi_sharcfx \
	audiomark/ports/adi_sharcfx/libs/NN/Source/ConvolutionFunctions \
	audiomark/ports/adi_sharcfx/libs/NN/Source/FullyConnectedFunctions \
	audiomark/ports/adi_sharcfx/libs/NN/Source/NNSupportFunctions \
	audiomark/ports/adi_sharcfx/libs/NN/Source/PoolingFunctions \
	audiomark/ports/adi_sharcfx/libs/NN/Source/SoftmaxFunctions \
	audiomark/ports/adi_sharcfx/libs/DSP/Source/BasicMathFunctions \
	audiomark/ports/adi_sharcfx/libs/DSP/Source/CommonTables \
	audiomark/ports/adi_sharcfx/libs/DSP/Source/ComplexMathFunctions \
	audiomark/ports/adi_sharcfx/libs/DSP/Source/FastMathFunctions \
	audiomark/ports/adi_sharcfx/libs/DSP/Source/MatrixFunctions \
	audiomark/ports/adi_sharcfx/libs/DSP/Source/StatisticsFunctions \
	audiomark/ports/adi_sharcfx/libs/DSP/Source/SupportFunctions \
	audiomark/ports/adi_sharcfx/libs/DSP/Source/TransformFunctions \
	audiomark/src \
	audiomark/src/ee_data/wav \


# Include directories
SHARCFX_INCLUDE_DIRS = \
	-I"../../audiomark/lib/speexdsp/include" \
	-I"../../audiomark/lib/speexdsp/include/speex" \
	-I"../../audiomark/lib/speexdsp/libspeexdsp" \
	-I"../../audiomark/src" \
	-I"../../audiomark/src/ee_data" \
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
	$(SRC_PREFIX)/audiomark/main.c \
	$(SRC_PREFIX)/audiomark/lib/speexdsp/libspeexdsp/filterbank_opt_generic.c \
	$(SRC_PREFIX)/audiomark/lib/speexdsp/libspeexdsp/filterbank_opt_helium.c \
	$(SRC_PREFIX)/audiomark/lib/speexdsp/libspeexdsp/filterbank_opt.c \
	$(SRC_PREFIX)/audiomark/lib/speexdsp/libspeexdsp/filterbank_opt_sharcfx.c \
	$(SRC_PREFIX)/audiomark/lib/speexdsp/libspeexdsp/mdf_opt_generic.c \
	$(SRC_PREFIX)/audiomark/lib/speexdsp/libspeexdsp/mdf_opt_helium.c \
	$(SRC_PREFIX)/audiomark/lib/speexdsp/libspeexdsp/preprocess_opt_generic.c \
	$(SRC_PREFIX)/audiomark/lib/speexdsp/libspeexdsp/preprocess_opt_helium.c \
	$(SRC_PREFIX)/audiomark/lib/speexdsp/libspeexdsp/testdenoise.c \
	$(SRC_PREFIX)/audiomark/lib/speexdsp/libspeexdsp/testecho.c \
	$(SRC_PREFIX)/audiomark/lib/speexdsp/libspeexdsp/testjitter.c \
	$(SRC_PREFIX)/audiomark/lib/speexdsp/libspeexdsp/testresample.c \
	$(SRC_PREFIX)/audiomark/lib/speexdsp/libspeexdsp/testresample2.c

# SHARCFX_C_SRC := $(foreach srcdir, $(SHARCFX_SRC_DIRS), $(wildcard $(SRC_PREFIX)/$(srcdir)/*.c))
SHARCFX_C_SRC := $(filter-out $(EXCLUDE_FILE),$(foreach srcdir, $(SHARCFX_SRC_DIRS), $(wildcard $(SRC_PREFIX)/$(srcdir)/*.c)))
SHARCFX_C_OBJ := $(subst $(SRC_PREFIX)/,,${SHARCFX_C_SRC:%.c=%.o})

# retrieve dependencies
SHARCFX_DEPS = $(foreach depdir, $(SHARCFX_SRC_DIRS), $(wildcard $(depdir)/*.d))

# define the compiler flags
SHARCFX_CFLAGS = $(SHARCFX_OPTIMIZE) $(BUILD_RELEASE) $(SHARCFX_INCLUDE_DIRS)
SHARCFX_CFLAGS += -proc $(PROC) -si-revision $(SI_REVISION) -c -LNO:simd -ffunction-sections -fdata-sections -mlongcalls -g
SHARCFX_CFLAGS += -DSHARC_FX -DNDEBUG -DFLOATING_POINT -D__RESTRICT= -DOS_SUPPORT_CUSTOM -DCORE1 -DEXPORT=
SHARCFX_CFLAGS += -DOVERRIDE_SPEEX_ALLOC -DOVERRIDE_MDF_SMOOTHED_ADD -DOVERRIDE_MDF_VEC_SCALE -DOVERRIDE_MDF_VEC_MULT -DOVERRIDE_MDF_VEC_ADD -DOVERRIDE_MDF_VEC_SUB -DOVERRIDE_MDF_POWER_SPECTRUM -DOVERRIDE_MDF_INNER_PROD
SHARCFX_CFLAGS += -DOVERRIDE_MDF_PREEMPH_FLT -DOVERRIDE_MDF_SMOOTH_FE_NRG -DOVERRIDE_MDF_POWER_SPECTRUM_ACCUM -DOVERRIDE_MDF_SPECTRAL_MUL_ACCUM -DOVERRIDE_MDF_VEC_CLEAR -DOVERRIDE_MDF_VEC_SUB_INT16 -D OVERRIDE_MDF_DEEMPH
SHARCFX_CFLAGS += -DOVERRIDE_MDF_WEIGHT_SPECT_MUL_CONJ -DOVERRIDE_MDF_FILTERED_SPEC_AD_XCORR -DOVERRIDE_MDF_ADJUST_PROP -DOVERRIDE_MDF_CONVERG_LEARN_RATE_CALC -DOVERRIDE_MDF_NORM_LEARN_RATE_CALC -DOVERRIDE_MDF_DC_NOTCH
SHARCFX_CFLAGS += -DOVERRIDE_ANR_VEC_COPY -DOVERRIDE_ANR_VEC_MUL -DOVERRIDE_ANR_OLA -DOVERRIDE_ANR_UPDATE_NOISE_ESTIMATE -DOVERRIDE_ANR_APOSTERIORI_SNR -DOVERRIDE_ANR_POWER_SPECTRUM -DOVERRIDE_ANR_COMPUTE_GAIN_FLOOR
SHARCFX_CFLAGS += -DOVERRIDE_ANR_VEC_CONV_FROM_INT16 -DOVERRIDE_ANR_UPDATE_ZETA -DOVERRIDE_ANR_HYPERGEOM_GAIN -DOVERRIDE_ANR_UPDATE_GAINS_CRITICAL_BANDS -DOVERRIDE_ANR_UPDATE_GAINS_LINEAR -DOVERRIDE_ANR_UPDATE_NOISE_PROB
SHARCFX_CFLAGS += -DOVERRIDE_ANR_APPLY_SPEC_GAIN -DOVERRIDE_FB_COMPUTE_BANK32 -DOVERRIDE_FB_COMPUTE_PSD16

SHARCFX_OBJS = $(SHARCFX_ASM_OBJ) $(SHARCFX_C_OBJ)

# compile SHARCFX 'C' files
audiomark/%.o: $(SRC_PREFIX)/audiomark/%.c
	$(SHARCFX_CC) $(SHARCFX_CFLAGS) -I"$(SRC_PREFIX)/$(dir $@)" -MMD -MP -MF "$(basename $@).d" -o "$@" -c "$<"

SHARCFX_LIB = $(PROJECT)-cces-library.a

.DEFAULT_GOAL = all
all: builddirs $(SHARCFX_LIB)

# Create static library from object files
$(SHARCFX_LIB): $(SHARCFX_OBJS)
	$(AR) -rcs -o $@ $^

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