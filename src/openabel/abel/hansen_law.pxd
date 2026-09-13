

from openabel.abel.base cimport abel_plan



cdef int plan_fat_hansen_law_org_lin(abel_plan* plan) except -1 nogil
cdef int execute_fat_hansen_law_linear(abel_plan* plan, double* data_in, double* data_out) except -1 nogil
cdef int destroy_fat_hansen_law_linear(abel_plan* plan) except -1 nogil
