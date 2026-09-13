

from openabel.abel.base cimport abel_plan



cdef int plan_fat_trapezoidal_desing_const(abel_plan* plan) except -1 nogil
cdef int execute_fat_trapezoidal_desing_const(abel_plan* plan, double* data_in, double* data_out, int left_boundary, int right_boundary) except -1 nogil
cdef int destroy_fat_trapezoidal_desing_const(abel_plan* plan) except -1 nogil

cdef int plan_fat_trapezoidal_end_corr(abel_plan* plan, int order = ?) except -1 nogil
cdef int execute_fat_trapezoidal_end_corr(abel_plan* plan, double* data_in, double* data_out, int left_boundary, int right_boundary) except -1 nogil
cdef int destroy_fat_trapezoidal_end_corr(abel_plan* plan) except -1 nogil
