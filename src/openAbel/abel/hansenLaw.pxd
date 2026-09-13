

from openAbel.abel.base cimport abel_plan



cdef int plan_fat_hansenLawOrgLin(abel_plan* plan) nogil
cdef int execute_fat_hansenLawLinear(abel_plan* plan, double* dataIn, double* dataOut) nogil
cdef int destroy_fat_hansenLawLinear(abel_plan* plan) nogil
