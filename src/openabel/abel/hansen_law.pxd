

from openabel.abel.base cimport abel_plan



cdef int plan_fat_hansenLawOrgLin(abel_plan* plan) except -1 nogil
cdef int execute_fat_hansenLawLinear(abel_plan* plan, double* dataIn, double* dataOut) except -1 nogil
cdef int destroy_fat_hansenLawLinear(abel_plan* plan) except -1 nogil
