# user library

SET(USER_src 
		${SRCUSER}/user_type.f90
		${SRCUSER}/covariate.f90
		${SRCUSER}/setUserParameters.f90
		${SRCUSER}/simulatePhenotype.f90
)
ADD_LIBRARY(${USERLIB} STATIC ${USER_src})

TARGET_LINK_LIBRARIES(${USERLIB} ${MATHLIB} ${BLAS_LIBRARIES}
	${LAPACK_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
