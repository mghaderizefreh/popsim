# blup library

SET(BLUP_src ${SRCBLUP}/BSRibsCalc.f90
    	     ${SRCBLUP}/relationMatrix.f90
			 ${SRCBLUP}/makegenotype.f90
	     ${SRCBLUP}/blup_module.f90
)
ADD_LIBRARY(${BLUPLIB} STATIC ${BLUP_src})

TARGET_LINK_LIBRARIES(${BLUPLIB} ${MATHLIB} ${SORTLIB}
				 ${BLAS_LIBRARIES}
                                 ${LAPACK_LIBRARIES}
                                 ${CMAKE_THREAD_LIBS_INIT})
