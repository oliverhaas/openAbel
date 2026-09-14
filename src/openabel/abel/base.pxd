

ctypedef struct abel_plan:
    int n_data
    int forward_backward
    double shift
    double step_size
    int method
    double* grid
    void* method_data

cdef abel_plan* plan_fat(int n_data, int forward_backward, double shift, double step_size, int method = ?, int order = ?, double eps = ?) except NULL nogil
cdef int execute_fat(abel_plan* plan, double* data_in, double* data_out, int left_boundary = ?, int right_boundary = ?) except -1 nogil
cdef int destroy_fat(abel_plan* plan) except -1 nogil
