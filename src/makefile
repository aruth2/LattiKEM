INC_DIR     = Core
SRC_DIR     = Core
OBJ_DIR     = Executables
BG_DIR      = EnergyModels/Bandgap
SRS         = ${SRC_DIR}/crystal.c ${SRC_DIR}/crystalnetwork.c ${SRC_DIR}/supp.c ${SRC_DIR}/matrix.c ${SRC_DIR}/perovskite.c ${SRC_DIR}/pkmc.c ${SRC_DIR}/latticedynamics.c ${SRC_DIR}/lattikem.c ${SRC_DIR}/settings.c ${SRC_DIR}/sourcecontrol.c ${BG_DIR}/mixedhalide.c ${BG_DIR}/bandgap.c ${EF_DIR}/polarization.c ${EF_DIR}/electrostatics.c ${EF_DIR}/chargecationhalidevacancies.c ${EF_DIR}/chargedvacanciespolarFAI.c ${STACD_DIR}/stacd.c ${STACD_DIR}/devicestack.c 
DEPS        = ${INC_DIR}/crystal.h ${INC_DIR}/crystalnetwork.h ${INC_DIR}/supp.h ${INC_DIR}/matrix.h ${INC_DIR}/perovskite.h ${INC_DIR}/pkmc.h ${INC_DIR}/latticedynamics.h ${INC_DIR}/lattikem.h ${INC_DIR}/settings.h ${SRC_DIR}/sourcecontrol.h ${BG_DIR}/bandgap.h ${EF_DIR}/polarization.h ${EF_DIR}/electrostatics.h ${STACD_DIR}/stacd.h
OBJS        = ${OBJ_DIR}/crystal.o ${OBJ_DIR}/crystalnetwork.o ${OBJ_DIR}/supp.o ${OBJ_DIR}/matrix.o ${OBJ_DIR}/perovskite.o ${OBJ_DIR}/pkmc.o ${OBJ_DIR}/latticedynamics.o ${OBJ_DIR}/lattikem.o ${OBJ_DIR}/settings.o ${OBJ_DIR}/sourcecontrol.o
MH_OBJS     = ${OBJS} ${BG_DIR}/mixedhalide.o ${BG_DIR}/bandgap.o
MH_GPU_OBJS = ${MH_OBJS} ${BG_DIR}/bandgap_GPU.o
TH_OBJS     = ${OBJS} ${BG_DIR}/triplehalide.o ${BG_DIR}/bandgap.o
TH_GPU_OBJS = ${TH_OBJS} ${BG_DIR}/bandgap_GPU.o
MCH_OBJS    = ${OBJS} ${BG_DIR}/mixedcationhalide.o ${BG_DIR}/bandgap.o
all_OBJS    = ${MH_OBJS} ${MH_GPU_OBJS} ${MCH_OBJS}

$(info Compiling Lattikem target $(MAKECMDGOALS)!)

ifeq ($(MAKECMDGOALS),$(filter $(MAKECMDGOALS),mixedhalide_GPU triplehalide_GPU))
cc=nvcc
else
cc=gcc
endif

ifeq ($(cc),nvcc)
$(info Using Cuda Compiler nvcc)
CUDAFLAGS = --default-stream per-thread -diag-suppress 2464
CFLAGS = -x cu ${CUDAFLAGS}
else
$(info Using GNU Compiler gcc)

ifeq ($(MAKECMDGOALS),mixedhalide_profile)
CFLAGS = -O3 -Wno-incompatible-pointer-types -Wno-format-zero-length -Wno-format-security -Wno-int-to-pointer-cast -Wno-pointer-to-int-cast -no-pie -pg
else
CFLAGS = -O3 -Wno-incompatible-pointer-types -Wno-format-zero-length -Wno-format-security -Wno-int-to-pointer-cast -Wno-pointer-to-int-cast
endif
endif


clean::
	${RM} ${all_OBJS} 
	
${OBJ_DIR}/%.o: ${SRC_DIR}/%.c
	${cc} ${CFLAGS} -I. -c $< -o $@

${SRC_DIR}/%.o: ${SRC_DIR}/%.c
	${cc} ${CFLAGS} -I. -c $< -o $@

${BG_DIR}/%.o: ${BG_DIR}/%.c
	${cc} ${CFLAGS} -I${SRC_DIR} -c $< -o $@				

mixedhalide: ${MH_OBJS}
	-${cc} ${CFLAGS} -o ${OBJ_DIR}/mixedhalide ${MH_OBJS} ${ACCFLAGS} -lm -lpthread -lnuma
	${RM} ${MH_OBJS}
	
triplehalide: ${TH_OBJS}
	-${cc} ${CFLAGS} -o ${OBJ_DIR}/triplehalide ${TH_OBJS} ${ACCFLAGS} -lm -lpthread -lnuma
	${RM} ${TH_OBJS}
	
mixedhalide_profile: ${MH_OBJS}
	-${cc} ${CFLAGS} -o ${OBJ_DIR}/mixedhalide_profile ${MH_OBJS} ${ACCFLAGS} -lm -lpthread -lnuma
	${RM} ${MH_OBJS}	

${BG_DIR}/bandgap_GPU.o: ${BG_DIR}/bandgap_GPU.cu
	${cc} ${CFLAGS} -I${SRC_DIR} -c $< -o $@

mixedhalide_GPU: ${MH_GPU_OBJS}
	-${cc} ${CUDAFLAGS} -o ${OBJ_DIR}/mixedhalide_GPU ${MH_GPU_OBJS} -lm -lpthread -lnuma
	${RM} ${MH_GPU_OBJS}
	
triplehalide_GPU: ${TH_GPU_OBJS}
	-${cc} ${CUDAFLAGS} -o ${OBJ_DIR}/triplehalide_GPU ${TH_GPU_OBJS} -lm -lpthread -lnuma
	${RM} ${TH_GPU_OBJS}

mixedcationhalide: ${MCH_OBJS}
	-${cc} -o ${OBJ_DIR}/mixedcationhalide ${MCH_OBJS} -lm -lpthread
	${RM} ${MCH_OBJS}

.DEFAULT_GOAL := all
 
test:
	Executables/mixedhalide ../Examples/mixedhalide/4x4x4
 
profile: 
	gprof Executables/mixedhalide_profile gmon.out > runstats

GPUprofile:
	nsys profile -t nvtx,openacc --stats=true --force-overwrite true -o GPUstats Executables/mixedhalide_GPU
