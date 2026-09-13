

from openabel.abel.base cimport abel_plan



cdef int plan_fat_trapezoidalDesingConst(abel_plan* plan) except -1 nogil
cdef int execute_fat_trapezoidalDesingConst(abel_plan* plan, double* dataIn, double* dataOut, int leftBoundary, int rightBoundary) except -1 nogil
cdef int destroy_fat_trapezoidalDesingConst(abel_plan* plan) except -1 nogil

cdef int plan_fat_trapezoidalEndCorr(abel_plan* plan, int order = ?) except -1 nogil
cdef int execute_fat_trapezoidalEndCorr(abel_plan* plan, double* dataIn, double* dataOut, int leftBoundary, int rightBoundary) except -1 nogil
cdef int destroy_fat_trapezoidalEndCorr(abel_plan* plan) except -1 nogil
