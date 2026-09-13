

ctypedef struct abel_plan:
    int nData
    int forwardBackward
    double shift
    double stepSize
    int method
    double* grid
    void* methodData

cdef abel_plan* plan_fat(int nData, int forwardBackward, double shift, double stepSize, int method = ?, int order = ?, double eps = ?) except NULL nogil
cdef int execute_fat(abel_plan* plan, double* dataIn, double* dataOut, int leftBoundary = ?, int rightBoundary = ?) except -1 nogil
cdef int destroy_fat(abel_plan* plan) except -1 nogil
