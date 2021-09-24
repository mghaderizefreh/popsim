# quick sort library

SET(SORT_src ${SRCSORT}/quickSort.f90
)

ADD_LIBRARY(${SORTLIB} STATIC ${SORT_src} ${UTIL_src})
