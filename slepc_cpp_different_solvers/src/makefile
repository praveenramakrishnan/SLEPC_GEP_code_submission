CFLAGS     =
FFLAGS     =
CPPFLAGS   =

main: main.o
			-${CLINKER} -o main main.o ${SLEPC_SYS_LIB}
			${RM} main.o
           
include ${SLEPC_DIR}/lib/slepc/conf/slepc_common
