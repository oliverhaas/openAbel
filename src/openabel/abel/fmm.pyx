

import numpy as np

from libc.stdlib cimport free
from openabel.helper cimport null_check_malloc as malloc, null_check_calloc as calloc
from libc.string cimport memset

cimport scipy.linalg.cython_blas as blas

import openabel.abel.coeffs as cffs

cimport openabel.math_fun as mf
cimport openabel.constants as co
from openabel.abel.base cimport abel_plan


########################################################################################################################
### Fast multipole method trapezoidal with end corrections                                                           ###


ctypedef struct method_data_fmm:
    (int*) kl, kl_cum
    (double*) cheb_roots, mtmp, mtmm, ltp, mtlk, direct, direct0, coeffs_sing, coeffs_nonsing, coeffs_filter
    int pp, pp1, ss, nlevs, k_total, order, order_filter

cdef double _kern_forward(double rr, double yy) nogil:
    return rr/mf.sqrt(rr**2-yy**2)
    
cdef double _kern_backward(double rr, double yy) nogil:
    return 1/mf.sqrt(rr**2-yy**2)

cdef double _kern_modified(double rr, double yy) nogil:
    return (yy/rr)**2/mf.sqrt(rr**2-yy**2)

# Plan FMM
cdef int plan_fat_fmm_trap_end_corr(abel_plan* pl, int order = 2, double eps = co.machine_epsilon) except -1 nogil:

    cdef:
        int ii, jj, ll, kk, mm
        method_data_fmm* md
        double ti, tauj, temp, n_inv_sca, y_inv_sca
        double[:,::1] cffs_s_sm_mv, cffs_s_la_mv, cffs_ns_sqrt_sm_mv, cffs_ns_sqrt_la_mv
        double[::1] cffs_f_mv
        int order_m1_half, order_m1_half_inner, y_cross, y_large, y_inv_sca_int, n_cross, n_large, n_inv_sca_int
        double (*kern)(double, double) nogil

    # Input check
    if NULL == pl or order <= 0 or eps < co.machine_epsilon:
        with gil:
            raise ValueError('Illegal input argument.')   

    # Main method struct
    md = <method_data_fmm*> malloc(sizeof(method_data_fmm))
    # Initialize to NULL so I can destroy properly when exception is raised
    md.cheb_roots = md.kl = md.kl_cum = md.mtmp = md.mtmm = md.mtlk = NULL 
    md.direct = md.direct0 = md.ltp = md.coeffs_sing = md.coeffs_nonsing = md.coeffs_filter = NULL
    pl.method_data = <void*> md

    # Small data set
    if pl.n_data < order+2:
        destroy_fat_fmm_trap_end_corr(pl)
        with gil:
            raise ValueError('Not enough data points for given parameters.')

    # Load and prepare end correction coefficients
    md.order = order
    order_m1_half = <int> ((md.order-1)/2)
    order_m1_half_inner = <int> (md.order/2)
    md.coeffs_sing = <double*> malloc(md.order*(pl.n_data-1)*sizeof(double))
    md.coeffs_nonsing = <double*> malloc(md.order*(pl.n_data-1)*sizeof(double))
    if pl.forward_backward == -1:    # Forward transform
        kern = &_kern_forward
        with gil:
            try:
                cffs_s_la_mv = cffs.get_coeffs('coeffs_inv_sqrt_diff_sq_lin_sing_large', order)
                cffs_ns_sqrt_sm_mv = cffs.get_coeffs('coeffs_inv_sqrt_nonsing_small', order)
                cffs_ns_sqrt_la_mv = cffs.get_coeffs('coeffs_inv_sqrt_nonsing_large', order)
                if pl.shift == 0.:
                    cffs_s_sm_mv = cffs.get_coeffs('coeffs_inv_sqrt_diff_sq_lin_sing_small', order)
                elif pl.shift == 0.5:
                    cffs_s_sm_mv = cffs.get_coeffs('coeffs_inv_sqrt_diff_sq_lin_sing_small_half_shift', order)
                else:
                    raise NotImplementedError('Method not implemented for given parameters.')
            except:
                destroy_fat_fmm_trap_end_corr(pl)
                raise
        y_cross = cffs_s_sm_mv.shape[0]
        y_large = cffs_s_la_mv.shape[0]
        for ii in range(min(y_cross,pl.n_data-1)):
            for jj in range(md.order):
                md.coeffs_sing[md.order*ii+jj] = cffs_s_sm_mv[ii,jj]
        for ii in range(y_cross, pl.n_data-1):
            y_inv_sca = pl.step_size/pl.grid[ii]*(y_cross-1)*(y_large-1)
            y_inv_sca_int = <int> mf.fmax(mf.fmin(y_inv_sca,y_large-3),1)
            for jj in range(md.order):
                md.coeffs_sing[md.order*ii+jj] = interp_cubic(y_inv_sca-y_inv_sca_int, md.order,
                                                            &cffs_s_la_mv[y_inv_sca_int-1,jj]) * \
                                                mf.sqrt(pl.grid[ii]/2./pl.step_size)
        n_cross = cffs_ns_sqrt_sm_mv.shape[0]            
        for ii in range(max(pl.n_data-1-n_cross,0),pl.n_data-1):
            for jj in range(md.order):
                md.coeffs_nonsing[md.order*ii+jj] = cffs_ns_sqrt_sm_mv[pl.n_data-2-ii,jj] * \
                                                   (pl.grid[pl.n_data-1]+(jj-order_m1_half_inner)*pl.step_size) / \
                                                   mf.sqrt((pl.grid[pl.n_data-1]+(jj-order_m1_half_inner)*pl.step_size+pl.grid[ii])*(pl.grid[pl.n_data-1]-pl.grid[ii]))
        n_large = cffs_ns_sqrt_la_mv.shape[0]      
        for ii in range(max(pl.n_data-1-n_cross,0)):
            n_inv_sca = pl.step_size/(pl.grid[pl.n_data-1]-pl.grid[ii])*n_cross*(n_large-1)
            n_inv_sca_int = <int> mf.fmax(mf.fmin(n_inv_sca,n_large-3),1)
            for jj in range(md.order):
                md.coeffs_nonsing[md.order*ii+jj] = interp_cubic(n_inv_sca-n_inv_sca_int, md.order, 
                                                               &cffs_ns_sqrt_la_mv[n_inv_sca_int-1,jj]) * \
                                                   (pl.grid[pl.n_data-1]+(jj-order_m1_half_inner)*pl.step_size) / \
                                                   mf.sqrt((pl.grid[pl.n_data-1] + \
                                                                 (jj-order_m1_half_inner)*pl.step_size+pl.grid[ii]) * \
                                                                (pl.grid[pl.n_data-1]-pl.grid[ii]))
        for ii in range(pl.n_data-1):
            md.coeffs_nonsing[md.order*ii+order_m1_half_inner] -= 0.5*kern(pl.grid[pl.n_data-1], pl.grid[ii])

    elif pl.forward_backward == 1 or pl.forward_backward == 2:    # Backward transform
        kern = &_kern_backward
        with gil:
            try:
                cffs_s_la_mv = cffs.get_coeffs('coeffs_inv_sqrt_diff_sq_sing_large', order)
                cffs_ns_sqrt_sm_mv = cffs.get_coeffs('coeffs_inv_sqrt_nonsing_small', order)
                cffs_ns_sqrt_la_mv = cffs.get_coeffs('coeffs_inv_sqrt_nonsing_large', order)
                if pl.shift == 0.:
                    cffs_s_sm_mv = cffs.get_coeffs('coeffs_inv_sqrt_diff_sq_sing_small', order)
                elif pl.shift == 0.5:
                    cffs_s_sm_mv = cffs.get_coeffs('coeffs_inv_sqrt_diff_sq_sing_small_half_shift', order)
                else:
                    raise NotImplementedError('Method not implemented for given parameters.')
            except:
                destroy_fat_fmm_trap_end_corr(pl)
                raise
        y_cross = cffs_s_sm_mv.shape[0]
        y_large = cffs_s_la_mv.shape[0]
        for ii in range(min(y_cross,pl.n_data-1)):
            for jj in range(md.order):
                md.coeffs_sing[md.order*ii+jj] = cffs_s_sm_mv[ii,jj]/pl.step_size
        for ii in range(y_cross, pl.n_data-1):
            y_inv_sca = pl.step_size/pl.grid[ii]*(y_cross-1)*(y_large-1)
            y_inv_sca_int = <int> mf.fmax(mf.fmin(y_inv_sca,y_large-3),1)
            for jj in range(md.order):
                md.coeffs_sing[md.order*ii+jj] = interp_cubic(y_inv_sca-y_inv_sca_int, md.order, 
                                                            &cffs_s_la_mv[y_inv_sca_int-1,jj]) / \
                                                mf.sqrt(pl.grid[ii]*2.*pl.step_size)
        n_cross = cffs_ns_sqrt_sm_mv.shape[0]            
        for ii in range(max(pl.n_data-1-n_cross,0),pl.n_data-1):
            for jj in range(md.order):
                md.coeffs_nonsing[md.order*ii+jj] = cffs_ns_sqrt_sm_mv[pl.n_data-2-ii,jj] / \
                                                   mf.sqrt((pl.grid[pl.n_data-1] + \
                                                            (jj-order_m1_half_inner)*pl.step_size+pl.grid[ii]) * \
                                                           (pl.grid[pl.n_data-1]-pl.grid[ii]))
        n_large = cffs_ns_sqrt_la_mv.shape[0]      
        for ii in range(max(pl.n_data-1-n_cross,0)):
            n_inv_sca = pl.step_size/(pl.grid[pl.n_data-1]-pl.grid[ii])*n_cross*(n_large-1)
            n_inv_sca_int = <int> mf.fmax(mf.fmin(n_inv_sca,n_large-3),1)
            for jj in range(md.order):
                md.coeffs_nonsing[md.order*ii+jj] = interp_cubic(n_inv_sca-n_inv_sca_int, md.order, 
                                                               &cffs_ns_sqrt_la_mv[n_inv_sca_int-1,jj]) / \
                                                   mf.sqrt((pl.grid[pl.n_data-1] + \
                                                                 (jj-order_m1_half_inner)*pl.step_size+pl.grid[ii]) * \
                                                                (pl.grid[pl.n_data-1]-pl.grid[ii]))
        for ii in range(pl.n_data-1):
            md.coeffs_nonsing[md.order*ii+order_m1_half_inner] -= 0.5*kern(pl.grid[pl.n_data-1],pl.grid[ii])

    elif pl.forward_backward == -2:    # Modified forward transform for 1/r^2 singular functions
        kern = &_kern_modified
        with gil:
            try:
                cffs_s_la_mv = cffs.get_coeffs('coeffs_inv_sqrt_diff_sq_y2_over_r2_sing_large', order)
                cffs_ns_sqrt_sm_mv = cffs.get_coeffs('coeffs_inv_sqrt_nonsing_small', order)
                cffs_ns_sqrt_la_mv = cffs.get_coeffs('coeffs_inv_sqrt_nonsing_large', order)
                if pl.shift == 0.:
                    cffs_s_sm_mv = cffs.get_coeffs('coeffs_inv_sqrt_diff_sq_y2_over_r2_sing_small', order)
                elif pl.shift == 0.5:
                    cffs_s_sm_mv = cffs.get_coeffs('coeffs_inv_sqrt_diff_sq_y2_over_r2_sing_small_half_shift', order)
                else:
                    raise NotImplementedError('Method not implemented for given parameters.')
            except:
                destroy_fat_fmm_trap_end_corr(pl)
                raise
        y_cross = cffs_s_sm_mv.shape[0]
        y_large = cffs_s_la_mv.shape[0]
        for ii in range(min(y_cross,pl.n_data-1)):
            for jj in range(md.order):
                md.coeffs_sing[md.order*ii+jj] = cffs_s_sm_mv[ii,jj]/pl.step_size
        for ii in range(y_cross, pl.n_data-1):
            y_inv_sca = pl.step_size/pl.grid[ii]*(y_cross-1)*(y_large-1)
            y_inv_sca_int = <int> mf.fmax(mf.fmin(y_inv_sca,y_large-3),1)
            for jj in range(md.order):
                md.coeffs_sing[md.order*ii+jj] = interp_cubic(y_inv_sca-y_inv_sca_int, md.order,
                                                            &cffs_s_la_mv[y_inv_sca_int-1,jj]) / \
                                                mf.sqrt(pl.grid[ii]*2.*pl.step_size)
        n_cross = cffs_ns_sqrt_sm_mv.shape[0]            
        for ii in range(max(pl.n_data-1-n_cross,0),pl.n_data-1):
            for jj in range(md.order):
                md.coeffs_nonsing[md.order*ii+jj] = cffs_ns_sqrt_sm_mv[pl.n_data-2-ii,jj] * \
                                                   (pl.grid[ii]/(pl.grid[pl.n_data-1]+(jj-order_m1_half_inner)*pl.step_size))**2 / \
                                                   mf.sqrt( (pl.grid[pl.n_data-1]+(jj-order_m1_half_inner)*pl.step_size+pl.grid[ii]) *
                                                                 (pl.grid[pl.n_data-1]-pl.grid[ii]) )
        n_large = cffs_ns_sqrt_la_mv.shape[0]      
        for ii in range(max(pl.n_data-1-n_cross,0)):
            n_inv_sca = pl.step_size/(pl.grid[pl.n_data-1]-pl.grid[ii])*n_cross*(n_large-1)
            n_inv_sca_int = <int> mf.fmax(mf.fmin(n_inv_sca,n_large-3),1)
            for jj in range(md.order):
                md.coeffs_nonsing[md.order*ii+jj] = interp_cubic(n_inv_sca-n_inv_sca_int, md.order, 
                                                               &cffs_ns_sqrt_la_mv[n_inv_sca_int-1,jj]) * \
                                                   (pl.grid[ii]/(pl.grid[pl.n_data-1]+(jj-order_m1_half_inner)*pl.step_size))**2 / \
                                                   mf.sqrt( (pl.grid[pl.n_data-1]+(jj-order_m1_half_inner)*pl.step_size+pl.grid[ii]) *
                                                                 (pl.grid[pl.n_data-1]-pl.grid[ii]) )
        for ii in range(pl.n_data-1):
            md.coeffs_nonsing[md.order*ii+order_m1_half_inner] -= 0.5*kern(pl.grid[pl.n_data-1],pl.grid[ii])

    else:
        destroy_fat_fmm_trap_end_corr(pl)
        with gil:
            raise NotImplementedError('Method not implemented for given parameters.')

    # Hierarchical decomposition
    md.pp = max(4, <int> ( -0.55*mf.log(2.*eps) + 1. ) )    # Empirical scaling, should be exponential
    md.pp1 = md.pp + 1
    md.nlevs = max(<int> ( mf.log2((pl.n_data-1.)/(2.*md.pp)) + 1. ), 2)
    md.ss = max(<int> ( (pl.n_data-1.)/2**md.nlevs + 1. ), 3)    # ss ~= 2*pp theoretical
    md.k_total = 2**(md.nlevs+1) - 1                             # Total number of intervals in all levels
    
    # Allocation of arrays FMM part
    md.cheb_roots = <double*> malloc(md.pp1*sizeof(double))
    md.kl = <int*> malloc((md.nlevs+1)*sizeof(int))
    md.kl_cum = <int*> malloc((md.nlevs+1)*sizeof(int))
    md.mtmp = <double*> malloc(md.pp1**2*sizeof(double))
    md.mtmm = <double*> malloc(md.pp1**2*sizeof(double))
    md.ltp = <double*> malloc(md.pp1*md.ss*sizeof(double))
    md.mtlk = <double*> calloc(2*md.k_total*md.pp1**2, sizeof(double))
    md.direct = <double*> calloc(2**md.nlevs*md.ss**2*2, sizeof(double))
    md.direct0 = <double*> malloc(2*md.ss*sizeof(double))

    # More hierarchical decomposition stuff
    md.kl_cum[0] = 0
    md.kl[0] = 2**md.nlevs
    for ii in range(1,md.nlevs+1):
        md.kl[ii] = 2**(md.nlevs-ii)
        md.kl_cum[ii] = md.kl_cum[ii-1] + md.kl[ii-1]
    _cheb_roots(md.pp1, md.cheb_roots)

    # Moment to moment coefficients
    for ii in range(md.pp1):
        for jj in range(md.pp1):
            md.mtmp[ii*md.pp1+jj] = _lagrange_p_int(0.5*md.cheb_roots[jj]+0.5, ii, md.cheb_roots, md.pp1)
            md.mtmm[ii*md.pp1+jj] = _lagrange_p_int(0.5*md.cheb_roots[jj]-0.5, ii, md.cheb_roots, md.pp1)    

    # Local to potential coefficients
    for ii in range(md.ss):
        temp = 2.*(ii+1)/md.ss-1.
        for jj in range(md.pp1):
            md.ltp[ii*md.pp1+jj] = _lagrange_p_int(temp, jj, md.cheb_roots, md.pp1)

    # Moment to local coefficients
    for ll in range(md.nlevs-1):
        # If even
        for kk in range(0, md.kl[ll]-2, 2):
            for ii in range(md.pp1):
                ti = (md.kl[md.nlevs-ll]*(kk+0.5+0.5*md.cheb_roots[ii])*md.ss + pl.shift)*pl.step_size
                for jj in range(md.pp1):
                    tauj = (md.kl[md.nlevs-ll]*(kk+2.5+0.5*md.cheb_roots[jj])*md.ss + pl.shift)*pl.step_size
                    md.mtlk[(2*(md.kl_cum[ll]+kk))*md.pp1**2+ii*md.pp1+jj] = kern(tauj, ti)
                    tauj = (md.kl[md.nlevs-ll]*(kk+3.5+0.5*md.cheb_roots[jj])*md.ss + pl.shift)*pl.step_size
                    md.mtlk[(2*(md.kl_cum[ll]+kk)+1)*md.pp1**2+ii*md.pp1+jj] = kern(tauj, ti)
        # If odd
        for kk in range(1, md.kl[ll]-2, 2):
            for ii in range(md.pp1):
                ti = (md.kl[md.nlevs-ll]*(kk+0.5+0.5*md.cheb_roots[ii])*md.ss + pl.shift)*pl.step_size
                for jj in range(md.pp1):
                    tauj = (md.kl[md.nlevs-ll]*(kk+2.5+0.5*md.cheb_roots[jj])*md.ss + pl.shift)*pl.step_size
                    md.mtlk[(2*(md.kl_cum[ll]+kk))*md.pp1**2+ii*md.pp1+jj] = kern(tauj, ti)

    # Direct short range coefficients
    for ii in range(1, pl.n_data):
        kk = <int> ((ii-1)/md.ss)
        ll = (ii-1) - kk*md.ss
        mm = min(pl.n_data-kk*md.ss-1, 2*md.ss)
        for jj in range(ll+1, mm):
            md.direct[kk*md.ss**2*2+md.ss*2*ll+jj] = kern(pl.grid[ii+jj-ll], pl.grid[ii])
    mm = min(pl.n_data, 2*md.ss+1)
    for jj in range(1, mm):
        md.direct0[jj-1] = kern(pl.grid[jj], pl.grid[0])    

    # Input modification filter
    if pl.forward_backward == 1:
        md.order_filter = md.order+1 + (md.order % 2)
        md.coeffs_filter = <double*> malloc(md.order_filter*sizeof(double))
        with gil:
            try:
                cffs_f_mv = cffs.get_coeffs('coeffs_deriv_smooth', md.order_filter-1)
            except:
                destroy_fat_fmm_trap_end_corr(pl)
                raise
        for ii in range(md.order_filter):
            md.coeffs_filter[ii] = cffs_f_mv[ii]*(-co.piinv)
    elif pl.forward_backward == 2 or pl.forward_backward == -1 or pl.forward_backward == -2:
        md.order_filter = 1
        md.coeffs_filter = <double*> malloc(1*sizeof(double))
        if pl.forward_backward == 2:            
            md.coeffs_filter[0] = -co.piinv*pl.step_size
        else:
            md.coeffs_filter[0] = 2.*pl.step_size

    return 0


# Execute FMM
cdef int execute_fat_fmm_trap_end_corr(abel_plan* pl, double* data_in, double* data_out, int left_boundary, 
                                    int right_boundary) except -1 nogil:

    cdef:
        method_data_fmm* md
        int ii, jj, kk, ll, mm, nn
        (double*) moments = NULL, local = NULL, data_in_temp0 = NULL, data_in_temp1 = NULL
        int order_m1_half, order_filter_m1_half, order_m1_half_inner, n_left_ext, n_right_ext

    md = <method_data_fmm*> pl.method_data
    order_m1_half = <int> ((md.order-1)/2)
    order_m1_half_inner = <int> (md.order/2)
    order_filter_m1_half = <int> ((md.order_filter-1)/2)

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
    if right_boundary == 0: # TODO maybe or right_boundary == 1 or right_boundary == 2:
        n_right_ext = order_m1_half + order_filter_m1_half
    elif right_boundary == 3:
        n_right_ext = 0
    else:
        free(data_in_temp0)
        free(data_in_temp1)
        with gil:
            raise NotImplementedError('Method not implemented for given parameters.')           
    # Copy and extend data if necessary
    nn = max(md.order, md.order_filter-1)
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
    convolve(data_in_temp0, pl.n_data+2*order_m1_half, data_in_temp1, md.order_filter, md.coeffs_filter)
    free(data_in_temp0)

    # Allocate temporary data arrays
    moments = <double*> calloc(md.k_total*md.pp1, sizeof(double))
    local = <double*> calloc(md.k_total*md.pp1, sizeof(double))

    # Finest level moment calculation
    mm = (pl.n_data-2)/md.ss - 2         # basically kl[0]-2, as first two block are not needed
    blas.dgemm('n', 'n', &md.pp1, &mm, &md.ss, &ONED, md.ltp, &md.pp1, 
               &data_in_temp1[order_m1_half+2*md.ss+1], &md.ss, &ZEROD, &moments[2*md.pp1], &md.pp1)
    mm += 2
    for ii in range(mm*md.ss+1, pl.n_data):
        nn = (ii-1) - mm*md.ss
        for jj in range(md.pp1):
            moments[mm*md.pp1+jj] += md.ltp[nn*md.pp1+jj]*data_in_temp1[order_m1_half+ii]

    # Upward Pass / Moment to moment
    mm = 2*md.pp1
    for ll in range(1, md.nlevs-1):
        blas.dgemm('t', 'n', &md.pp1, &md.kl[ll], &md.pp1, &ONED, md.mtmm, &md.pp1, 
                   &moments[md.kl_cum[ll-1]*md.pp1], &mm, &ONED, &moments[md.kl_cum[ll]*md.pp1], &md.pp1)
        blas.dgemm('t', 'n', &md.pp1, &md.kl[ll], &md.pp1, &ONED, md.mtmp, &md.pp1, 
                   &moments[(md.kl_cum[ll-1]+1)*md.pp1], &mm, &ONED, &moments[md.kl_cum[ll]*md.pp1], &md.pp1)

    # Interaction Phase / Moment to local
    for ll in range(md.nlevs-1):
        # If even
        for kk in range(0, md.kl[ll]-2, 2):
            blas.dgemv('t', &md.pp1, &md.pp1, &ONED, &md.mtlk[(2*(md.kl_cum[ll]+kk))*md.pp1**2], &md.pp1, 
                       &moments[(kk+md.kl_cum[ll]+2)*md.pp1], &ONE, &ZEROD, &local[(md.kl_cum[ll]+kk)*md.pp1], &ONE)
            blas.dgemv('t', &md.pp1, &md.pp1, &ONED, &md.mtlk[(2*(md.kl_cum[ll]+kk)+1)*md.pp1**2], &md.pp1, 
                       &moments[(kk+md.kl_cum[ll]+3)*md.pp1], &ONE, &ONED, &local[(md.kl_cum[ll]+kk)*md.pp1], &ONE)
        # If odd
        for kk in range(1, md.kl[ll]-2, 2):
            blas.dgemv('t', &md.pp1, &md.pp1, &ONED, &md.mtlk[(2*(md.kl_cum[ll]+kk))*md.pp1**2], &md.pp1, 
                       &moments[(kk+md.kl_cum[ll]+2)*md.pp1], &ONE, &ZEROD, &local[(md.kl_cum[ll]+kk)*md.pp1], &ONE)
    
    # Downward Pass / Local to local
    mm = 2*md.pp1
    for ll in range(md.nlevs-2, 0, -1):
        blas.dgemm('n', 'n', &md.pp1, &md.kl[ll], &md.pp1, &ONED, md.mtmm, &md.pp1, 
                   &local[md.kl_cum[ll]*md.pp1], &md.pp1, &ONED, &local[md.kl_cum[ll-1]*md.pp1], &mm)
        blas.dgemm('n', 'n', &md.pp1, &md.kl[ll], &md.pp1, &ONED, md.mtmp, &md.pp1, 
                   &local[md.kl_cum[ll]*md.pp1], &md.pp1, &ONED, &local[(md.kl_cum[ll-1]+1)*md.pp1], &mm)

    # Set output to zero first
    memset(data_out, 0, pl.n_data*sizeof(double))

    # Potential evaluation / local to potential
    mm = (pl.n_data-2)/md.ss - 1          # basically kl[0]-2, as last two block are not needed
    blas.dgemm('t', 'n', &md.ss, &mm, &md.pp1, &ONED, md.ltp, &md.pp1, local, &md.pp1, &ZEROD, &data_out[1], &md.ss)
    for jj in range(md.pp1):
        data_out[0] += local[jj]*_lagrange_p_int(-1., jj, md.cheb_roots, md.pp1)
    
    free(moments)
    free(local)

    # Direct short range
    ll = min((pl.n_data-2)/md.ss, md.kl[0]-1)
    mm = 2*md.ss
    for kk in range(ll):
        # Only the rows that hold data: rows beyond the data end are zero in md.direct, and reading the matching
        # input elements would run past data_in_temp1.
        nn = min(pl.n_data-kk*md.ss-1, mm)
        blas.dgemv('t', &nn, &md.ss, &ONED, &md.direct[kk*md.ss**2*2], &mm,
                   &data_in_temp1[order_m1_half+1+kk*md.ss], &ONE, &ONED, &data_out[1+kk*md.ss], &ONE)
    for ii in range(ll*md.ss+1, pl.n_data-1):
        kk = (ii-1)/md.ss
        nn = (ii-1) - kk*md.ss
        mm = min(pl.n_data-kk*md.ss-1, 2*md.ss)
        for jj in range(nn+1, mm):
            data_out[ii] += md.direct[kk*md.ss**2*2+md.ss*2*nn+jj]*data_in_temp1[order_m1_half+ii+jj-nn]
    mm = min(pl.n_data, 2*md.ss+1)
    for jj in range(1, mm):
        data_out[0] += md.direct0[jj-1]*data_in_temp1[order_m1_half+jj]

    # End correction left singular end
    # TODO maybe BLAS
    for ii in range(pl.n_data-1):
        for jj in range(md.order):
            data_out[ii] += md.coeffs_sing[md.order*ii+jj]*data_in_temp1[ii+jj]

    # End correction right nonsingular end
    mm = pl.n_data-1
    blas.dgemv('t', &md.order, &mm, &ONED, md.coeffs_nonsing, &md.order, 
               &data_in_temp1[pl.n_data-1+order_m1_half-order_m1_half_inner], &ONE, &ONED, data_out, &ONE)

    free(data_in_temp1)

    return 0


# Destroy FMM
cdef int destroy_fat_fmm_trap_end_corr(abel_plan* pl) except -1 nogil:

    cdef:
        method_data_fmm* md = <method_data_fmm*> pl.method_data

    free(md.cheb_roots)
    free(md.kl)
    free(md.kl_cum)
    free(md.mtmp)
    free(md.mtmm) 
    free(md.ltp)    
    free(md.mtlk)
    free(md.direct)
    free(md.direct0)  
    free(md.coeffs_sing)
    free(md.coeffs_nonsing)
    free(md.coeffs_filter)
    free(md)

    return 0



#############################################################################################################################################

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
    # TODO Maybe DGEMM or FFT here
    for ii in range(n_data):
        for jj in range(order):
            data_out[ii] += coeffs[jj]*data_in[ii+jj]

    return 0


cdef int _cheb_roots(int order, double* roots) nogil:
    
    cdef:
        int ii
    
    for ii in range(order):
        roots[ii] = mf.cos(0.5*co.pi*(2.*ii+1.)/order)
    
    return 0


cdef double _lagrange_p_int(double xx, unsigned int ind, double* nodes, unsigned int order) nogil:

    cdef:
        unsigned int ii
        double res
    
    if xx >= -1. and xx <= 1.:
        res = 1.
        for ii in range(ind):
            res *= (xx - nodes[ii])/(nodes[ind] - nodes[ii])
        for ii in range(ind+1, order):
            res *= (xx - nodes[ii])/(nodes[ind] - nodes[ii])
    else:
        res = 0.

    return res


cdef:
    int ONE = 1
    double ZEROD = 0.
    double ONED = 1.
    double TWOD = 2.
    double MONED = -1.
    double MTWOD = -2.
