

from libc.stdlib cimport free
from openabel.helper cimport null_check_malloc as malloc

from openabel.abel.hansen_law cimport plan_fat_hansen_law_org_lin, execute_fat_hansen_law_linear, destroy_fat_hansen_law_linear
from openabel.abel.trap cimport plan_fat_trapezoidal_desing_const, execute_fat_trapezoidal_desing_const, \
                                destroy_fat_trapezoidal_desing_const, \
                                plan_fat_trapezoidal_end_corr, execute_fat_trapezoidal_end_corr, \
                                destroy_fat_trapezoidal_end_corr
from openabel.abel.fmm cimport plan_fat_fmm_trap_end_corr, execute_fat_fmm_trap_end_corr, destroy_fat_fmm_trap_end_corr

cimport openabel.constants as const


ctypedef struct abel_plan:
    int n_data, forward_backward, method
    double shift, step_size
    double* grid
    void* method_data



########################################################################################################################
### Fast Abel transforms                                                                                             ###
########################################################################################################################


# Create plan for Abel transform
cdef abel_plan* plan_fat(int n_data, int forward_backward, double shift, double step_size, 
                         int method = 3, int order = 2, double eps = 1.e3*const.machine_epsilon) except NULL nogil:

    cdef:
        abel_plan* pl
        int ii

    pl = <abel_plan*> malloc(sizeof(abel_plan))
    pl.method_data = NULL
    pl.n_data = n_data
    pl.forward_backward = forward_backward
    pl.shift = shift
    pl.step_size = step_size
    pl.method = method
    pl.grid = <double*> malloc(n_data*sizeof(double))
    for ii in range(n_data):
        pl.grid[ii] = (ii+shift)*step_size

    with gil:
        try:
            if pl.method == 0:
                plan_fat_trapezoidal_desing_const(pl)
            elif pl.method == 1:
                plan_fat_hansen_law_org_lin(pl)
            elif pl.method == 2:
                plan_fat_trapezoidal_end_corr(pl, order = order)
            elif pl.method == 3:
                plan_fat_fmm_trap_end_corr(pl, order = order, eps = eps)
            else:
                raise NotImplementedError('Method not implemented for given parameters.')
        except:
            free(pl.grid)
            free(pl)
            raise

    return pl


# Execute given plan for Abel transform
cdef int execute_fat(abel_plan* pl, double* data_in, double* data_out, int left_boundary = 0, 
                     int right_boundary = 0) except -1 nogil:

    if NULL == pl:
        with gil:
            raise TypeError('Input plan is NULL.')

    if pl.method == 0:
        execute_fat_trapezoidal_desing_const(pl, data_in, data_out, left_boundary, right_boundary)
    elif pl.method == 1:
        execute_fat_hansen_law_linear(pl, data_in, data_out)
    elif pl.method == 2:
        execute_fat_trapezoidal_end_corr(pl, data_in, data_out, left_boundary, right_boundary)
    elif pl.method == 3:
        execute_fat_fmm_trap_end_corr(pl, data_in, data_out, left_boundary, right_boundary)
    else:
        with gil:
            raise NotImplementedError('Method not implemented for given parameters.')


# Destroy given plan for Abel transform
cdef int destroy_fat(abel_plan* pl) except -1 nogil:

    if NULL == pl:
        return 1

    if pl.method == 0:
        destroy_fat_trapezoidal_desing_const(pl)
    elif pl.method == 1:
        destroy_fat_hansen_law_linear(pl)
    elif pl.method == 2:
        destroy_fat_trapezoidal_end_corr(pl)
    elif pl.method == 3:
        destroy_fat_fmm_trap_end_corr(pl)
    else:
        with gil:
            raise NotImplementedError('Method not implemented for given parameters.')

    free(pl.grid)
    free(pl)    

    return 0
