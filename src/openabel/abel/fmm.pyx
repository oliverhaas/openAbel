from libc.stdlib cimport free
from openabel.helper cimport null_check_malloc as malloc, null_check_calloc as calloc
from libc.string cimport memset

cimport scipy.linalg.cython_blas as blas

from libc.math cimport sqrt, log, log2, cos
cimport openabel.constants as co
from openabel.abel.base cimport abel_plan
from openabel.abel.trap cimport plan_end_corr_coeffs, plan_filter_coeffs, prepare_input_end_corr


########################################################################################################################
### Fast multipole method trapezoidal with end corrections                                                           ###


ctypedef struct method_data_fmm:
    (int*) kl, kl_cum
    (double*) cheb_roots, mtmp, mtmm, ltp, mtlk, direct, direct0, coeffs_sing, coeffs_nonsing, coeffs_filter
    int pp, pp1, ss, nlevs, k_total, order, order_filter

cdef double _kern_forward(double rr, double yy) nogil:
    return rr/sqrt(rr**2-yy**2)
    
cdef double _kern_backward(double rr, double yy) nogil:
    return 1/sqrt(rr**2-yy**2)

cdef double _kern_modified(double rr, double yy) nogil:
    return (yy/rr)**2/sqrt(rr**2-yy**2)

# Plan FMM
cdef int plan_fat_fmm_trap_end_corr(abel_plan* pl, int order = 2, double eps = 1.e1*co.machine_epsilon) except -1 nogil:

    cdef:
        int ii, jj, ll, kk, mm
        method_data_fmm* md
        double ti, tauj, temp
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
    md.order = order
    pl.method_data = <void*> md

    # End correction coefficients and input modification filter
    with gil:
        try:
            plan_end_corr_coeffs(pl, order, &md.coeffs_sing, &md.coeffs_nonsing)
            plan_filter_coeffs(pl, order, &md.coeffs_filter, &md.order_filter)
        except:
            destroy_fat_fmm_trap_end_corr(pl)
            raise

    # Kernel of the transform; the shared preparation rejected every other forward_backward value
    if pl.forward_backward == -1:
        kern = &_kern_forward
    elif pl.forward_backward == -2:
        kern = &_kern_modified
    else:
        kern = &_kern_backward

    # Hierarchical decomposition
    md.pp = max(4, <int> ( -0.55*log(2.*eps) + 1. ) )    # Empirical scaling, should be exponential
    md.pp1 = md.pp + 1
    md.nlevs = max(<int> ( log2((pl.n_data-1.)/(2.*md.pp)) + 1. ), 2)
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

    return 0


# Execute FMM
cdef int execute_fat_fmm_trap_end_corr(abel_plan* pl, double* data_in, double* data_out, int left_boundary, 
                                    int right_boundary) except -1 nogil:

    cdef:
        method_data_fmm* md
        int ii, jj, kk, ll, mm, nn
        (double*) moments = NULL, local = NULL, data_in_temp1 = NULL
        int order_m1_half, order_m1_half_inner

    md = <method_data_fmm*> pl.method_data
    order_m1_half = <int> ((md.order-1)/2)
    order_m1_half_inner = <int> (md.order/2)
    data_in_temp1 = prepare_input_end_corr(pl, data_in, md.order, md.order_filter, md.coeffs_filter, left_boundary,
                                           right_boundary)

    # Allocate temporary data arrays
    moments = <double*> calloc(md.k_total*md.pp1, sizeof(double))
    local = <double*> calloc(md.k_total*md.pp1, sizeof(double))

    # Finest level moment calculation. For n_data < 3*ss+2 there is no complete block after the first two, and BLAS
    # complains about a negative dimension instead of doing nothing.
    mm = (pl.n_data-2)/md.ss - 2         # basically kl[0]-2, as first two block are not needed
    if mm > 0:
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
    if mm > 0:
        blas.dgemm('t', 'n', &md.ss, &mm, &md.pp1, &ONED, md.ltp, &md.pp1, local, &md.pp1, &ZEROD, &data_out[1],
                   &md.ss)
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



########################################################################################################################

cdef int _cheb_roots(int order, double* roots) nogil:
    
    cdef:
        int ii
    
    for ii in range(order):
        roots[ii] = cos(0.5*co.pi*(2.*ii+1.)/order)
    
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
