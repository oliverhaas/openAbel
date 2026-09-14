from libc.stdlib cimport free
from openabel.helper cimport null_check_malloc as malloc
from libc.string cimport memset

import openabel.abel.coeffs as coeffs

from libc.math cimport sqrt, log, fmin, fmax
cimport openabel.constants as co
from openabel.abel.base cimport abel_plan





############################################################################################################################################
### Fast Abel transforms based on trapezoidal rules                                                                                      ###
############################################################################################################################################


############################################################################################################################################
### Trapezoidal rule with constant desingularization                                                                                     ###

ctypedef struct method_data_desing_const:
    double* desing
    double* coeffs_filter
    int order_filter

# Plan desingularized quadrature trapezoidal
cdef int plan_fat_trapezoidal_desing_const(abel_plan* pl) except -1 nogil:

    cdef:
        method_data_desing_const* md
        int ii
        double temp0, temp1

    # Input check
    if NULL == pl:
        with gil:
            raise ValueError('Illegal input argument.')   

    # Small data set
    if pl.n_data < 3:
        with gil:
            raise ValueError('Not enough data points for given parameters.')

    # Main method struct
    md = <method_data_desing_const*> malloc(sizeof(method_data_desing_const))
    md.desing = NULL
    md.coeffs_filter = NULL
    pl.method_data = <void*> md

    # Desingularize array
    md.desing = <double*> malloc((pl.n_data-1)*sizeof(double))
    if pl.forward_backward == -1:
        for ii in range(pl.n_data-1):
            temp0 = sqrt(pl.grid[pl.n_data-1]**2-pl.grid[ii]**2)
            md.desing[ii] = temp0/pl.step_size
    elif pl.forward_backward == 2 or pl.forward_backward == 1:
        for ii in range(1,pl.n_data-1):
            temp0 = sqrt(pl.grid[pl.n_data-1]**2-pl.grid[ii]**2)
            temp1 = log((pl.grid[pl.n_data-1]+temp0)/pl.grid[ii])
            md.desing[ii] = temp1/pl.step_size
        temp0 = sqrt(pl.grid[pl.n_data-1]**2-pl.grid[0]**2)
        if pl.shift == 0.:
            md.desing[0] = 0.
        else:
            temp1 = log((pl.grid[pl.n_data-1]+temp0)/pl.grid[0])
            md.desing[0] = temp1/pl.step_size
    elif pl.forward_backward == -2:
        for ii in range(pl.n_data-1):
            md.desing[ii] = sqrt(pl.grid[pl.n_data-1]**2-pl.grid[ii]**2)/pl.grid[pl.n_data-1]/pl.step_size
    else:
        destroy_fat_trapezoidal_desing_const(pl)
        with gil:
            raise NotImplementedError('Method not implemented for given parameters.')

    # Input modification filter, as for first order end corrections
    with gil:
        try:
            plan_filter_coeffs(pl, 1, &md.coeffs_filter, &md.order_filter)
        except:
            destroy_fat_trapezoidal_desing_const(pl)
            raise

    return 0


# Execute desingularized quadrature trapezoidal
cdef int execute_fat_trapezoidal_desing_const(abel_plan* pl, double* data_in, double* data_out, int left_boundary, 
                                            int right_boundary) except -1 nogil:

    cdef:
        int ii, jj
        method_data_desing_const* md
        double* data_in_temp1

    md = <method_data_desing_const*> pl.method_data
    data_in_temp1 = prepare_input_end_corr(pl, data_in, 1, md.order_filter, md.coeffs_filter, left_boundary,
                                           right_boundary)

    # Main trapezoidal rule
    memset(data_out, 0, pl.n_data*sizeof(double))
    if pl.forward_backward == -1:
        for ii in range(pl.n_data-1):
            for jj in range(ii+1, pl.n_data-1):
                data_out[ii] += (data_in_temp1[jj]-data_in_temp1[ii]) * pl.grid[jj]/sqrt(pl.grid[jj]**2-pl.grid[ii]**2)
            jj = pl.n_data-1
            data_out[ii] += 0.5*(data_in_temp1[jj]-data_in_temp1[ii]) * pl.grid[jj]/sqrt(pl.grid[jj]**2-pl.grid[ii]**2)
            data_out[ii] += data_in_temp1[ii]*md.desing[ii]
    elif pl.forward_backward == 1 or pl.forward_backward == 2:
        for ii in range(pl.n_data-1):
            for jj in range(ii+1, pl.n_data-1):
                data_out[ii] += (data_in_temp1[jj]-data_in_temp1[ii]) / sqrt(pl.grid[jj]**2-pl.grid[ii]**2)
            jj = pl.n_data-1
            data_out[ii] += 0.5*(data_in_temp1[jj]-data_in_temp1[ii]) / sqrt(pl.grid[jj]**2-pl.grid[ii]**2)
            data_out[ii] += data_in_temp1[ii]*md.desing[ii]
    elif pl.forward_backward == -2:
        for ii in range(pl.n_data-1):
            for jj in range(ii+1, pl.n_data-1):
                data_out[ii] += (data_in_temp1[jj]-data_in_temp1[ii]) * (pl.grid[ii]/pl.grid[jj])**2 / \
                               sqrt(pl.grid[jj]**2-pl.grid[ii]**2)
            jj = pl.n_data-1
            data_out[ii] += 0.5*(data_in_temp1[pl.n_data-1]-data_in_temp1[ii]) * (pl.grid[ii]/pl.grid[pl.n_data-1])**2 / \
                           sqrt(pl.grid[pl.n_data-1]**2-pl.grid[ii]**2)
            data_out[ii] += data_in_temp1[ii]*md.desing[ii]
    else:
        free(data_in_temp1)
        with gil:
            raise NotImplementedError('Method not implemented for given parameters.')

    free(data_in_temp1)

    return 0


cdef int destroy_fat_trapezoidal_desing_const(abel_plan* pl) except -1 nogil:

    cdef:
        method_data_desing_const* md = <method_data_desing_const*> pl.method_data

    free(md.desing)
    free(md.coeffs_filter)
    free(md)

    return 0


########################################################################################################################
### Trapezoidal rule with end corrections                                                                            ###


ctypedef struct method_data_end_corr:
    double* coeffs_sing
    double* coeffs_nonsing
    double* coeffs_filter
    int order
    int order_filter


# Plan end-corrected trapezoidal
cdef int plan_fat_trapezoidal_end_corr(abel_plan* pl, int order = 2) except -1 nogil:

    cdef:
        method_data_end_corr* md

    # Input check
    if NULL == pl or order <= 0:
        with gil:
            raise ValueError('Illegal input argument.')

    # Main method struct
    md = <method_data_end_corr*> malloc(sizeof(method_data_end_corr))
    md.coeffs_sing = NULL
    md.coeffs_nonsing = NULL
    md.coeffs_filter = NULL
    md.order = order
    pl.method_data = <void*> md

    # End correction coefficients and input modification filter
    with gil:
        try:
            plan_end_corr_coeffs(pl, order, &md.coeffs_sing, &md.coeffs_nonsing)
            plan_filter_coeffs(pl, order, &md.coeffs_filter, &md.order_filter)
        except:
            destroy_fat_trapezoidal_end_corr(pl)
            raise

    return 0


# Execute end-corrected trapezoidal
cdef int execute_fat_trapezoidal_end_corr(abel_plan* pl, double* data_in, double* data_out, int left_boundary, 
                                        int right_boundary) except -1 nogil:

    cdef:
        int ii, jj, nn
        method_data_end_corr* md
        double* data_in_temp1
        int order_m1_half, order_m1_half_inner

    md = <method_data_end_corr*> pl.method_data
    order_m1_half = <int> ((md.order-1)/2)
    order_m1_half_inner = <int> (md.order/2)
    data_in_temp1 = prepare_input_end_corr(pl, data_in, md.order, md.order_filter, md.coeffs_filter, left_boundary,
                                           right_boundary)

    # Main trapezoidal rule
    memset(data_out, 0, pl.n_data*sizeof(double))
    if pl.forward_backward == -1:
        for ii in range(pl.n_data-1):
            for jj in range(ii+1, pl.n_data):
                data_out[ii] += data_in_temp1[order_m1_half+jj]*pl.grid[jj]/sqrt(pl.grid[jj]**2-pl.grid[ii]**2)
    elif pl.forward_backward == 1 or pl.forward_backward == 2:
        for ii in range(pl.n_data-1):
            for jj in range(ii+1, pl.n_data):
                data_out[ii] += data_in_temp1[order_m1_half+jj]/sqrt(pl.grid[jj]**2-pl.grid[ii]**2)
    elif pl.forward_backward == -2:
        for ii in range(pl.n_data-1):
            for jj in range(ii+1, pl.n_data):
                data_out[ii] += data_in_temp1[order_m1_half+jj]*(pl.grid[ii]/pl.grid[jj])**2/sqrt(pl.grid[jj]**2-pl.grid[ii]**2)
    else:
        free(data_in_temp1)
        with gil:
            raise NotImplementedError('Method not implemented for given parameters.')

    # End correction right fairly smooth end
    for ii in range(pl.n_data-1):
        for nn in range(md.order):
            jj = pl.n_data-1+nn-order_m1_half_inner+order_m1_half
            data_out[ii] += md.coeffs_nonsing[md.order*ii+nn]*data_in_temp1[jj]

    # End correction left singular end
    for ii in range(pl.n_data-1):
        for nn in range(md.order):
            data_out[ii] += md.coeffs_sing[md.order*ii+nn]*data_in_temp1[ii+nn]

    free(data_in_temp1)

    return 0


cdef int destroy_fat_trapezoidal_end_corr(abel_plan* pl) except -1 nogil:

    cdef:
        method_data_end_corr* md

    md = <method_data_end_corr*> pl.method_data
    free(md.coeffs_sing)
    free(md.coeffs_nonsing)
    free(md.coeffs_filter)
    free(md)

    return 0


########################################################################################################################
### Preparation shared by the trapezoidal rules and the FMM                                                          ###


# End correction weights of every output point for the singular end (the lower integration limit) and the smooth end
# (the truncation radius), order weights each, loaded from the coefficient tables and scaled to the grid of pl. The
# caller frees both arrays, also after an exception.
cdef int plan_end_corr_coeffs(abel_plan* pl, int order, double** coeffs_sing, double** coeffs_nonsing) except -1 nogil:

    cdef:
        double[:,::1] coeffs_nonsing_sqrt_small_mv
        double[:,::1] coeffs_nonsing_sqrt_large_mv
        double[:,::1] coeffs_sing_small_mv
        double[:,::1] coeffs_sing_large_mv
        double* cs
        double* cn
        int ii, jj
        double n_inv_sca, y_inv_sca
        int n_cross, n_large, n_inv_sca_int, y_cross, y_large, y_inv_sca_int
        int order_m1_half_inner

    # Small data set
    if pl.n_data < order+2:
        with gil:
            raise ValueError('Not enough data points for given parameters.')

    order_m1_half_inner = <int> (order/2)
    coeffs_sing[0] = <double*> malloc(order*(pl.n_data-1)*sizeof(double))
    coeffs_nonsing[0] = <double*> malloc(order*(pl.n_data-1)*sizeof(double))
    cs = coeffs_sing[0]
    cn = coeffs_nonsing[0]
    if pl.forward_backward == -1:    # Forward transform
        with gil:
            coeffs_sing_large_mv = coeffs.get_coeffs('coeffs_inv_sqrt_diff_sq_lin_sing_large', order)
            coeffs_nonsing_sqrt_small_mv = coeffs.get_coeffs('coeffs_inv_sqrt_nonsing_small', order)
            coeffs_nonsing_sqrt_large_mv = coeffs.get_coeffs('coeffs_inv_sqrt_nonsing_large', order)
            if pl.shift == 0.:
                coeffs_sing_small_mv = coeffs.get_coeffs('coeffs_inv_sqrt_diff_sq_lin_sing_small', order)
            elif pl.shift == 0.5:
                coeffs_sing_small_mv = coeffs.get_coeffs('coeffs_inv_sqrt_diff_sq_lin_sing_small_half_shift', order)
            else:
                raise NotImplementedError('Method not implemented for given parameters.')
        y_cross = coeffs_sing_small_mv.shape[0]
        y_large = coeffs_sing_large_mv.shape[0]
        for ii in range(min(y_cross,pl.n_data-1)):
            for jj in range(order):
                cs[order*ii+jj] = coeffs_sing_small_mv[ii,jj]
        for ii in range(y_cross, pl.n_data-1):
            y_inv_sca = pl.step_size/pl.grid[ii]*(y_cross-1)*(y_large-1)
            y_inv_sca_int = <int> fmax(fmin(y_inv_sca,y_large-3),1)
            for jj in range(order):
                cs[order*ii+jj] = interp_cubic(y_inv_sca-y_inv_sca_int, order, 
                                             &coeffs_sing_large_mv[y_inv_sca_int-1,jj]) * \
                                 sqrt(pl.grid[ii]/2./pl.step_size)
        n_cross = coeffs_nonsing_sqrt_small_mv.shape[0]            
        for ii in range(max(pl.n_data-1-n_cross,0),pl.n_data-1):
            for jj in range(order):
                cn[order*ii+jj] = coeffs_nonsing_sqrt_small_mv[pl.n_data-2-ii,jj] * \
                                 (pl.grid[pl.n_data-1]+(jj-order_m1_half_inner)*pl.step_size) / \
                                 sqrt((pl.grid[pl.n_data-1] + \
                                          (jj-order_m1_half_inner)*pl.step_size+pl.grid[ii]) * \
                                         (pl.grid[pl.n_data-1]-pl.grid[ii]))
        n_large = coeffs_nonsing_sqrt_large_mv.shape[0]      
        for ii in range(max(pl.n_data-1-n_cross,0)):
            n_inv_sca = pl.step_size/(pl.grid[pl.n_data-1]-pl.grid[ii])*n_cross*(n_large-1)
            n_inv_sca_int = <int> fmax(fmin(n_inv_sca,n_large-3),1)
            for jj in range(order):
                cn[order*ii+jj] = interp_cubic(n_inv_sca-n_inv_sca_int, order,
                                             &coeffs_nonsing_sqrt_large_mv[n_inv_sca_int-1,jj]) * \
                                 (pl.grid[pl.n_data-1]+(jj-order_m1_half_inner)*pl.step_size) / \
                                 sqrt((pl.grid[pl.n_data-1] + \
                                          (jj-order_m1_half_inner)*pl.step_size+pl.grid[ii]) * \
                                         (pl.grid[pl.n_data-1]-pl.grid[ii]))
        for ii in range(pl.n_data-1):
            cn[order*ii+order_m1_half_inner] -= 0.5*pl.grid[pl.n_data-1] / \
                                            sqrt(pl.grid[pl.n_data-1]**2-pl.grid[ii]**2)

    elif pl.forward_backward == 1 or pl.forward_backward == 2:    # Backward transform
        with gil:
            coeffs_sing_large_mv = coeffs.get_coeffs('coeffs_inv_sqrt_diff_sq_sing_large', order)
            coeffs_nonsing_sqrt_small_mv = coeffs.get_coeffs('coeffs_inv_sqrt_nonsing_small', order)
            coeffs_nonsing_sqrt_large_mv = coeffs.get_coeffs('coeffs_inv_sqrt_nonsing_large', order)
            if pl.shift == 0.:
                coeffs_sing_small_mv = coeffs.get_coeffs('coeffs_inv_sqrt_diff_sq_sing_small', order)
            elif pl.shift == 0.5:
                coeffs_sing_small_mv = coeffs.get_coeffs('coeffs_inv_sqrt_diff_sq_sing_small_half_shift', order)
            else:
                raise NotImplementedError('Method not implemented for given parameters.')
        y_cross = coeffs_sing_small_mv.shape[0]
        y_large = coeffs_sing_large_mv.shape[0]
        for ii in range(min(y_cross,pl.n_data-1)):
            for jj in range(order):
                cs[order*ii+jj] = coeffs_sing_small_mv[ii,jj]/pl.step_size
        for ii in range(y_cross, pl.n_data-1):
            y_inv_sca = pl.step_size/pl.grid[ii]*(y_cross-1)*(y_large-1)
            y_inv_sca_int = <int> fmax(fmin(y_inv_sca,y_large-3),1)
            for jj in range(order):
                cs[order*ii+jj] = interp_cubic(y_inv_sca-y_inv_sca_int, order,
                                             &coeffs_sing_large_mv[y_inv_sca_int-1,jj]) / \
                                 sqrt(pl.grid[ii]*2.*pl.step_size)
        n_cross = coeffs_nonsing_sqrt_small_mv.shape[0]            
        for ii in range(max(pl.n_data-1-n_cross,0),pl.n_data-1):
            for jj in range(order):
                cn[order*ii+jj] = coeffs_nonsing_sqrt_small_mv[pl.n_data-2-ii,jj] / \
                                 sqrt((pl.grid[pl.n_data-1] + \
                                          (jj-order_m1_half_inner)*pl.step_size+pl.grid[ii]) * \
                                         (pl.grid[pl.n_data-1]-pl.grid[ii]))
        n_large = coeffs_nonsing_sqrt_large_mv.shape[0]      
        for ii in range(max(pl.n_data-1-n_cross,0)):
            n_inv_sca = pl.step_size/(pl.grid[pl.n_data-1]-pl.grid[ii])*n_cross*(n_large-1)
            n_inv_sca_int = <int> fmax(fmin(n_inv_sca,n_large-3),1)
            for jj in range(order):
                cn[order*ii+jj] = interp_cubic(n_inv_sca-n_inv_sca_int, order,
                                             &coeffs_nonsing_sqrt_large_mv[n_inv_sca_int-1,jj]) / \
                                 sqrt((pl.grid[pl.n_data-1] + \
                                          (jj-order_m1_half_inner)*pl.step_size+pl.grid[ii]) * \
                                         (pl.grid[pl.n_data-1]-pl.grid[ii]))
        for ii in range(pl.n_data-1):
            cn[order*ii+order_m1_half_inner] -= 0.5/sqrt(pl.grid[pl.n_data-1]**2-pl.grid[ii]**2)

    elif pl.forward_backward == -2:    # Modified forward transform for 1/r^2 singular functions
        with gil:
            coeffs_sing_large_mv = coeffs.get_coeffs('coeffs_inv_sqrt_diff_sq_y2_over_r2_sing_large', order)
            coeffs_nonsing_sqrt_small_mv = coeffs.get_coeffs('coeffs_inv_sqrt_nonsing_small', order)
            coeffs_nonsing_sqrt_large_mv = coeffs.get_coeffs('coeffs_inv_sqrt_nonsing_large', order)
            if pl.shift == 0.:
                coeffs_sing_small_mv = coeffs.get_coeffs('coeffs_inv_sqrt_diff_sq_y2_over_r2_sing_small', order)
            elif pl.shift == 0.5:
                coeffs_sing_small_mv = coeffs.get_coeffs('coeffs_inv_sqrt_diff_sq_y2_over_r2_sing_small_half_shift', order)
            else:
                raise NotImplementedError('Method not implemented for given parameters.')
        y_cross = coeffs_sing_small_mv.shape[0]
        y_large = coeffs_sing_large_mv.shape[0]
        for ii in range(min(y_cross,pl.n_data-1)):
            for jj in range(order):
                cs[order*ii+jj] = coeffs_sing_small_mv[ii,jj]/pl.step_size
        for ii in range(y_cross, pl.n_data-1):
            y_inv_sca = pl.step_size/pl.grid[ii]*(y_cross-1)*(y_large-1)
            y_inv_sca_int = <int> fmax(fmin(y_inv_sca,y_large-3),1)
            for jj in range(order):
                cs[order*ii+jj] = interp_cubic(y_inv_sca-y_inv_sca_int, order, 
                                             &coeffs_sing_large_mv[y_inv_sca_int-1,jj]) / \
                                 sqrt(pl.grid[ii]*2.*pl.step_size)
        n_cross = coeffs_nonsing_sqrt_small_mv.shape[0]            
        for ii in range(max(pl.n_data-1-n_cross,0),pl.n_data-1):
            for jj in range(order):
                cn[order*ii+jj] = coeffs_nonsing_sqrt_small_mv[pl.n_data-2-ii,jj] * \
                                 (pl.grid[ii]/(pl.grid[pl.n_data-1]+(jj-order_m1_half_inner)*pl.step_size))**2 / \
                                 sqrt( (pl.grid[pl.n_data-1]+(jj-order_m1_half_inner)*pl.step_size+pl.grid[ii]) *
                                               (pl.grid[pl.n_data-1]-pl.grid[ii]) )
        n_large = coeffs_nonsing_sqrt_large_mv.shape[0]      
        for ii in range(max(pl.n_data-1-n_cross,0)):
            n_inv_sca = pl.step_size/(pl.grid[pl.n_data-1]-pl.grid[ii])*n_cross*(n_large-1)
            n_inv_sca_int = <int> fmax(fmin(n_inv_sca,n_large-3),1)
            for jj in range(order):
                cn[order*ii+jj] = interp_cubic(n_inv_sca-n_inv_sca_int, order, &coeffs_nonsing_sqrt_large_mv[n_inv_sca_int-1,jj]) * \
                                 (pl.grid[ii]/(pl.grid[pl.n_data-1]+(jj-order_m1_half_inner)*pl.step_size))**2 / \
                                 sqrt( (pl.grid[pl.n_data-1]+(jj-order_m1_half_inner)*pl.step_size+pl.grid[ii]) *
                                               (pl.grid[pl.n_data-1]-pl.grid[ii]) )
        for ii in range(pl.n_data-1):
            cn[order*ii+order_m1_half_inner] -= 0.5*(pl.grid[ii]/pl.grid[pl.n_data-1])**2/sqrt(pl.grid[pl.n_data-1]**2-pl.grid[ii]**2)

    else:
        with gil:
            raise NotImplementedError('Method not implemented for given parameters.')

    return 0


# Filter applied to the input before the quadrature: the scaling of the forward transforms, or the smoothed
# numerical derivative for the backward transform with numerical derivative (order_filter points). The caller frees
# coeffs_filter, also after an exception.
cdef int plan_filter_coeffs(abel_plan* pl, int order, double** coeffs_filter, int* order_filter) except -1 nogil:

    cdef:
        int ii
        double[::1] coeffs_filter_mv

    if pl.forward_backward == 1:
        order_filter[0] = order+1 + (order % 2)
        coeffs_filter[0] = <double*> malloc(order_filter[0]*sizeof(double))
        with gil:
            coeffs_filter_mv = coeffs.get_coeffs('coeffs_deriv_smooth', order_filter[0]-1)
        for ii in range(order_filter[0]):
            coeffs_filter[0][ii] = coeffs_filter_mv[ii]*(-co.piinv)
    elif pl.forward_backward == 2 or pl.forward_backward == -1 or pl.forward_backward == -2:
        order_filter[0] = 1
        coeffs_filter[0] = <double*> malloc(1*sizeof(double))
        if pl.forward_backward == 2:            
            coeffs_filter[0][0] = -co.piinv*pl.step_size
        else:
            coeffs_filter[0][0] = 2.*pl.step_size
    else:
        with gil:
            raise NotImplementedError('Method not implemented for given parameters.')

    return 0


# Copy the input into a buffer extended by the reach of the end correction stencil and the filter on each side
# (extrapolated for boundary value 0, mirrored for 1 and 2, taken from the input for 3) and apply the filter. Returns
# the filtered buffer of n_data + 2*((order-1)/2) samples; the caller frees it.
cdef double* prepare_input_end_corr(abel_plan* pl, double* data_in, int order, int order_filter, double* coeffs_filter,
                                    int left_boundary, int right_boundary) except NULL nogil:

    cdef:
        int ii, nn
        (double*) data_in_temp0 = NULL, data_in_temp1 = NULL
        int order_m1_half, order_filter_m1_half, n_left_ext, n_right_ext

    order_m1_half = <int> ((order-1)/2)
    order_filter_m1_half = <int> ((order_filter-1)/2)

    # Allocate temporary data arrays
    data_in_temp0 = <double*> malloc((pl.n_data+2*(order_m1_half+order_filter_m1_half))*sizeof(double))
    data_in_temp1 = <double*> malloc((pl.n_data+2*order_m1_half)*sizeof(double))
    
    # Left boundary handling
    if left_boundary == 0 or left_boundary == 1 or left_boundary == 2:
        n_left_ext = order_m1_half + order_filter_m1_half
    elif left_boundary == 3:
        n_left_ext = 0
    else:
        free(data_in_temp0)
        free(data_in_temp1)
        with gil:
            raise NotImplementedError('Method not implemented for given parameters.')
    # Right boundary handling
    if right_boundary == 0:
        n_right_ext = order_m1_half + order_filter_m1_half
    elif right_boundary == 3:
        n_right_ext = 0
    else:
        free(data_in_temp0)
        free(data_in_temp1)
        with gil:
            raise NotImplementedError('Method not implemented for given parameters.')           
    # Copy and extend data if necessary
    nn = max(order, order_filter-1)
    for ii in range(pl.n_data+2*(order_m1_half+order_filter_m1_half)-n_left_ext-n_right_ext):
        data_in_temp0[n_left_ext+ii] = data_in[ii]
    if left_boundary == 0:
        for ii in range(n_left_ext):
            data_in_temp0[ii] = polint(&data_in_temp0[n_left_ext], nn, ii-n_left_ext)
    elif left_boundary == 1:
        if pl.shift == 0.:
            for ii in range(n_left_ext):
                data_in_temp0[n_left_ext-1-ii] = -data_in_temp0[n_left_ext+1+ii]
        elif pl.shift == 0.5:
            for ii in range(n_left_ext):
                data_in_temp0[n_left_ext-1-ii] = -data_in_temp0[n_left_ext+ii]
        else:
            free(data_in_temp0)
            free(data_in_temp1)
            with gil:
                raise NotImplementedError('Method not implemented for given parameters.')
    elif left_boundary == 2:
        if pl.shift == 0.:
            for ii in range(n_left_ext):
                data_in_temp0[n_left_ext-1-ii] = data_in_temp0[n_left_ext+1+ii]
        elif pl.shift == 0.5:
            for ii in range(n_left_ext):
                data_in_temp0[n_left_ext-1-ii] = data_in_temp0[n_left_ext+ii]
        else:
            free(data_in_temp0)
            free(data_in_temp1)
            with gil:
                raise NotImplementedError('Method not implemented for given parameters.')
    elif left_boundary == 3:
        pass
    else:
        free(data_in_temp0)
        free(data_in_temp1)
        with gil:
            raise NotImplementedError('Method not implemented for given parameters.')
    if right_boundary == 0:
        for ii in range(n_right_ext):
            data_in_temp0[pl.n_data+order_m1_half+order_filter_m1_half+ii] = polint(&data_in_temp0[pl.n_data+order_m1_half+order_filter_m1_half-nn], nn, ii+nn)
    elif right_boundary == 3:
        pass
    else:
        free(data_in_temp0)
        free(data_in_temp1)
        with gil:
            raise NotImplementedError('Method not implemented for given parameters.')

    # Do scaling or numerical derivative
    convolve(data_in_temp0, pl.n_data+2*order_m1_half, data_in_temp1, order_filter, coeffs_filter)
    free(data_in_temp0)

    return data_in_temp1


########################################################################################################################
# Cubic interpolation
cdef inline double interp_cubic(double x, int incx, double* p) nogil:

    return p[1*incx] + 0.5 * x*(p[2*incx] - p[0*incx] + 
                                x*(2.0*p[0*incx] - 5.0*p[1*incx] + 4.0*p[2*incx] - p[3*incx] + 
                                   x*(3.0*(p[1*incx] - p[2*incx]) + p[3*incx] - p[0*incx])))


# Polynomial inter-/extrapolation on equidistant grid
cdef double polint(double* data, int n_data, double xx) nogil:

    cdef:
        int ii, mm, ns
        double* cc = NULL
        double* dd = NULL
        double den, res

    cc = <double*> malloc(n_data*sizeof(double))
    dd = <double*> malloc(n_data*sizeof(double))
    ns = <int> (xx+0.5)
    ns = min(max(ns,0),n_data-1)
    for ii in range(n_data):
        cc[ii] = data[ii]
        dd[ii] = data[ii]

    res = data[ns]
    ns -= 1
    for mm in range(1,n_data):
        for ii in range(n_data-mm):
            den = (dd[ii]-cc[ii+1])/mm
            dd[ii] = (ii+mm-xx)*den
            cc[ii]= (ii-xx)*den
        if 2*(ns+1) < n_data-mm:
            res += cc[ns+1]
        else:
            res += dd[ns]
            ns -= 1
    free(cc)
    free(dd)

    return res


# Apply filter; possibly just numerical derivative
cdef int convolve(double* data_in, int n_data, double* data_out, int order, double* coeffs) nogil:

    cdef:
        int ii, jj

    memset(data_out, 0, n_data*sizeof(double))
    for ii in range(n_data):
        for jj in range(order):
            data_out[ii] += coeffs[jj]*data_in[ii+jj]

    return 0


########################################################################################################################
