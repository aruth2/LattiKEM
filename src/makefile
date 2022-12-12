#This is the makefile for spectrum-slicing interface. 

INC_DIR     = Core
SRC_DIR     = Core
OBJ_DIR     = Executables
BG_DIR      = EnergyModels/Bandgap
EF_DIR      = EnergyModels/ElectricField
STACD_DIR   = EnergyModels/DeviceSimulator
SRS         = ${SRC_DIR}/crystal.c ${SRC_DIR}/crystalnetwork.c ${SRC_DIR}/supp.c ${SRC_DIR}/matrix.c ${SRC_DIR}/perovskite.c ${SRC_DIR}/pkmc.c ${SRC_DIR}/latticedynamics.c ${SRC_DIR}/lattikem.c ${SRC_DIR}/settings.c ${SRC_DIR}/sourcecontrol.c ${BG_DIR}/mixedhalide.c ${BG_DIR}/bandgap.c ${EF_DIR}/polarization.c ${EF_DIR}/electrostatics.c ${EF_DIR}/chargecationhalidevacancies.c ${EF_DIR}/chargedvacanciespolarFAI.c ${STACD_DIR}/stacd.c ${STACD_DIR}/devicestack.c 
DEPS        = ${INC_DIR}/crystal.h ${INC_DIR}/crystalnetwork.h ${INC_DIR}/supp.h ${INC_DIR}/matrix.h ${INC_DIR}/perovskite.h ${INC_DIR}/pkmc.h ${INC_DIR}/latticedynamics.h ${INC_DIR}/lattikem.h ${INC_DIR}/settings.h ${SRC_DIR}/sourcecontrol.h ${BG_DIR}/bandgap.h ${EF_DIR}/polarization.h ${EF_DIR}/electrostatics.h ${STACD_DIR}/stacd.h
OBJS        = ${OBJ_DIR}/crystal.o ${OBJ_DIR}/crystalnetwork.o ${OBJ_DIR}/supp.o ${OBJ_DIR}/matrix.o ${OBJ_DIR}/perovskite.o ${OBJ_DIR}/pkmc.o ${OBJ_DIR}/latticedynamics.o ${OBJ_DIR}/lattikem.o ${OBJ_DIR}/settings.o ${OBJ_DIR}/sourcecontrol.o 
MH_OBJS     = ${OBJS} ${BG_DIR}/mixedhalide.o ${BG_DIR}/bandgap.o
MCH_OBJS    = ${OBJS} ${BG_DIR}/mixedcationhalide.o ${BG_DIR}/bandgap.o
CCHV_OBJS   = ${OBJS} ${EF_DIR}/chargedcationhalidevacancies.o ${EF_DIR}/electrostatics.o ${EF_DIR}/polarization.o
CVPF_OBJS   = ${OBJS} ${EF_DIR}/chargedvacanciespolarFAI.o ${EF_DIR}/electrostatics.o ${EF_DIR}/polarization.o
STACD_OBJS  = ${OBJS} ${STACD_DIR}/stacd.o ${STACD_DIR}/devicestack.o

#suppresses a warning from CHKERRQ
CFLAGS     = -O3 -Wno-incompatible-pointer-types -Wno-format-zero-length -Wno-format-security -Wno-int-to-pointer-cast -Wno-pointer-to-int-cast

#For Profiling
#CFLAGS     = -Wno-incompatible-pointer-types -Wno-format-zero-length -Wno-format-security -Wno-int-to-pointer-cast -Wno-pointer-to-int-cast -no-pie -pg
CC=gcc


clean::
	${RM} ${MH_OBJS}
	
${OBJ_DIR}/%.o: ${SRC_DIR}/%.c
	${CC} ${CFLAGS} -I. -c $< -o $@
	
${BG_DIR}/%.o: ${BG_DIR}/%.c
	${CC} ${CFLAGS} -I${SRC_DIR} -c $< -o $@	

${EF_DIR}/%.o: ${EF_DIR}/%.c
	${CC} ${CFLAGS} -I${SRC_DIR} -c $< -o $@
	
${STACD_DIR}/%.o: ${STACD_DIR}/%.c
	${CC} ${CFLAGS} -I${SRC_DIR} -c $< -o $@		

mixedhalide: ${MH_OBJS}
	-${CC} -o ${OBJ_DIR}/mixedhalide ${MH_OBJS} -lm -lpthread
	${RM} ${MH_OBJS}

mixedcationhalide: ${MCH_OBJS}
	-${CC} -o ${OBJ_DIR}/mixedcationhalide ${MCH_OBJS} -lm -lpthread
	${RM} ${MCH_OBJS}

#Currently pretty broken. Need to fix since upgrades in core modules.
chargedcationhalidevacancies: ${CCHV_OBJS}
	-${CC} -o ${OBJ_DIR}/chargedcationhalidevacancies ${CCHV_OBJS} -lm -lpthread
	${RM} ${CCHV_OBJS}

#Currently pretty broken. Need to fix since upgrades in core modules.	
chargedvacanciespolarFAI: ${CVPF_OBJS}
	-${CC} -o ${OBJ_DIR}/chargedvacanciespolarFAI ${CVPF_OBJS} -lm -lpthread
	${RM} ${CVPF_OBJS}	

stacd: ${STACD_OBJS}
	-${CC} -o ${OBJ_DIR}/devicestack ${STACD_OBJS} -lm -lpthread
	${RM} ${STACD_OBJS}	

.DEFAULT_GOAL := all
 
profile: 
	gprof ./mixedhalide gmon.out > runstats

