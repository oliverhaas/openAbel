

from libc.string cimport memset
cdef extern from "stdlib.h":
    void* aligned_alloc(size_t alignment, size_t size) nogil


# Inline malloc with null check, kind of a simple "hacky" solution not to have to do null check every time manually.
# Good enough for me and saves a lot of lines.
# Just "cimport [...] nullCheckMalloc as malloc" to replace normal malloc
# https://stackoverflow.com/questions/26831981/should-i-check-if-malloc-was-successful/26844703
# I decided to force alignment for up to AVX512 here, since it's usually worth it and not much lost if not.
# Might change this in the future. So for very specific cases alignment should be chosen manually anyway.

cdef:
    size_t stdAlgn = max(sizeof(void*), 64)     # Ensures that nothing will completely break in the future (I hope).


cdef inline void* nullCheckMalloc(size_t MemSize, size_t alignment = stdAlgn) nogil:

    cdef:
        void* AllocMem = aligned_alloc(alignment, MemSize)
    
    if NULL == AllocMem and MemSize > 0:
        with gil:
            print('Malloc returned NULL pointer, probably not enough memory or wrong user input alignment.\
                   Will exit because there is no good way to handle this and continue.')
            exit(-1)
    
    return AllocMem


cdef inline void* nullCheckCalloc(size_t nn, size_t size, size_t alignment = stdAlgn) nogil:

    cdef:
        void* AllocMem = aligned_alloc(alignment, nn*size)
    
    if NULL == AllocMem and nn*size > 0:
        with gil:
            print('Malloc returned NULL pointer, probably not enough memory or wrong user input alignment.\
                   Will exit because there is no good way to handle this and continue.')
            exit(-1)
    
    memset(AllocMem, 0, nn*size)
    
    return AllocMem









