from libc.string cimport memset
from libc.stdint cimport SIZE_MAX
cdef extern from "stdlib.h":
    void* aligned_alloc(size_t alignment, size_t size) nogil


# Inline malloc with null check, kind of a simple "hacky" solution not to have to do null check every time manually.
# Good enough for me and saves a lot of lines.
# Just "cimport [...] null_check_malloc as malloc" to replace normal malloc
# https://stackoverflow.com/questions/26831981/should-i-check-if-malloc-was-successful/26844703
# I decided to force alignment for up to AVX512 here, since it's usually worth it and not much lost if not.
# Might change this in the future. So for very specific cases alignment should be chosen manually anyway.
#
# C11 requires the size passed to aligned_alloc to be a multiple of the alignment (macOS enforces this, glibc does not),
# so the size is rounded up, and a zero-byte request allocates one alignment block so that the returned pointer is
# always valid and can be passed to free(). A size that would wrap during the rounding raises MemoryError as well. On
# failure a MemoryError is raised (the functions are "except NULL"), so the callers' cleanup paths run instead of
# writing through a NULL pointer.


cdef inline void* null_check_malloc(size_t mem_size, size_t alignment = 64) except NULL nogil:

    cdef:
        size_t alloc_size
        void* alloc_mem

    if mem_size > SIZE_MAX - (alignment - 1):
        with gil:
            raise MemoryError('Requested allocation size overflows size_t.')

    alloc_size = ((mem_size + alignment - 1) // alignment) * alignment

    if alloc_size == 0:
        alloc_size = alignment

    alloc_mem = aligned_alloc(alignment, alloc_size)

    if NULL == alloc_mem:
        with gil:
            raise MemoryError('aligned_alloc returned NULL, probably not enough memory or an invalid alignment.')

    return alloc_mem


cdef inline void* null_check_calloc(size_t nn, size_t size, size_t alignment = 64) except NULL nogil:

    cdef:
        void* alloc_mem

    if nn != 0 and size > SIZE_MAX // nn:
        with gil:
            raise MemoryError('Requested allocation size overflows size_t.')

    alloc_mem = null_check_malloc(nn*size, alignment)
    memset(alloc_mem, 0, nn*size)

    return alloc_mem
