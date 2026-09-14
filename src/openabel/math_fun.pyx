



cdef unsigned int uint_max(unsigned int aa, unsigned int bb) nogil:

    if aa > bb:
        return aa

    return bb


cdef unsigned int uint_min(unsigned int aa, unsigned int bb) nogil:

    if aa < bb:
        return aa

    return bb
