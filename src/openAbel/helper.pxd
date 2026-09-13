from libc.string cimport memset
cdef extern from "stdlib.h":
    void* aligned_alloc(size_t alignment, size_t size) nogil


# Inline malloc with null check, kind of a simple "hacky" solution not to have to do null check every time manually.
# Good enough for me and saves a lot of lines.
# Just "cimport [...] nullCheckMalloc as malloc" to replace normal malloc
# https://stackoverflow.com/questions/26831981/should-i-check-if-malloc-was-successful/26844703
# I decided to force alignment for up to AVX512 here, since it's usually worth it and not much lost if not.
# Might change this in the future. So for very specific cases alignment should be chosen manually anyway.
#
# C11 requires the size passed to aligned_alloc to be a multiple of the alignment (macOS enforces this, glibc does
# not), so the size is rounded up, and a zero-byte request allocates one alignment block so that the returned pointer
# is always valid and can be passed to free(). On failure a MemoryError is raised (the functions are "except NULL"),
# so the callers' cleanup paths run instead of writing through a NULL pointer.


cdef inline void* nullCheckMalloc(size_t MemSize, size_t alignment = 64) except NULL nogil:

    cdef:
        size_t allocSize = ((MemSize + alignment - 1) // alignment) * alignment
        void* AllocMem

    if allocSize == 0:
        allocSize = alignment

    AllocMem = aligned_alloc(alignment, allocSize)

    if NULL == AllocMem:
        with gil:
            raise MemoryError('aligned_alloc returned NULL, probably not enough memory or an invalid alignment.')

    return AllocMem


cdef inline void* nullCheckCalloc(size_t nn, size_t size, size_t alignment = 64) except NULL nogil:

    cdef:
        void* AllocMem = nullCheckMalloc(nn*size, alignment)

    memset(AllocMem, 0, nn*size)

    return AllocMem
