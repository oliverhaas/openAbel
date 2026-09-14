cimport openabel.abel.base as base


cdef class Abel:

    cdef:
        base.abel_plan* plan
        int n_outside
