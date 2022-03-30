CFLAGS     =
FFLAGS     = -cpp -dM -I${SLEPC_DIR}/include -I${PETSC_DIR}/include
CPPFLAGS   =

main: main.o
			-${CLINKER} -o main *.o ${SLEPC_SYS_LIB}
			${RM} main.o
           
include ${SLEPC_DIR}/lib/slepc/conf/slepc_common
