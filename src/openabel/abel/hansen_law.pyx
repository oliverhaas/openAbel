from libc.stdlib cimport free
from openabel.helper cimport null_check_malloc as malloc

cimport openabel.constants as constants
from libc.math cimport log
from openabel.abel.base cimport abel_plan






########################################################################################################################
### Hansen Law space state model for Abel transform                                                                  ###


ctypedef struct method_data_hansen_law:
    double* coeffs
    model_hansen_law* model


ctypedef struct model_hansen_law:
    unsigned int nk
    double* hk
    double* lamk


# Original model (factor pi different definition than original Hansen Law paper)
cdef model_hansen_law model_hansen_law_org
model_hansen_law_org.nk = 9
model_hansen_law_org.hk = [1., 0.596903, 1.09956, 2.57611, 5.65487, 12.2522, 26.0752, 61.5752, 151.739]
model_hansen_law_org.lamk = [0.0, -2.1, -6.2, -22.4, -92.5, -414.5, -1889.4, -8990.9, -47391.1]



########################################################################################################################
# Plan Hansen Law original model (9th order linear)
cdef int plan_fat_hansen_law_org_lin(abel_plan* plan) except -1 nogil:
    
    cdef:
        method_data_hansen_law* method_data = <method_data_hansen_law*> malloc(sizeof(method_data_hansen_law))

    plan.method_data = <void*> method_data
    method_data.model = &model_hansen_law_org
    method_data.coeffs = NULL

    return plan_fat_hansen_law_linear(plan)

########################################################################################################################
# Plan Hansen Law linear
cdef int plan_fat_hansen_law_linear(abel_plan* plan) except -1 nogil:
  
    cdef:
        method_data_hansen_law* md = <method_data_hansen_law*> plan.method_data
        model_hansen_law* mod = md.model
        int ii, jj, kk, ind_start
        double xjp1oxj
          
    if plan.forward_backward == -1:
        md.coeffs = <double*> malloc(3*plan.n_data*(mod.nk-1)*sizeof(double))
        if plan.shift == 0.:
            for kk in range(1,mod.nk):
                jj = 3*(kk-1)
                md.coeffs[jj] = 0.
                md.coeffs[jj+1] = 0.
                md.coeffs[jj+2] = 0.
            ind_start = 1
        else:
            ind_start = 0
        for ii in range(ind_start, plan.n_data-1):
            xjp1oxj = plan.grid[ii+1]/plan.grid[ii]
            for kk in range(1,mod.nk):
                jj = 3*(mod.nk-1)*ii + 3*(kk-1)
                md.coeffs[jj] = xjp1oxj**mod.lamk[kk]
                md.coeffs[jj+1] = 2.*mod.hk[kk]*(1.-md.coeffs[jj]*xjp1oxj)*plan.grid[ii] / \
                                          (mod.lamk[kk]+1.)
                md.coeffs[jj+2] = 2.*mod.hk[kk]*(1.-md.coeffs[jj]*xjp1oxj**2)*plan.grid[ii]**2 / \
                                          (mod.lamk[kk]+2.)

    elif plan.forward_backward == 1 or plan.forward_backward == 2:

        md.coeffs = <double*> malloc(2*plan.n_data*(mod.nk-1)*sizeof(double))

        if plan.shift == 0.:
            for kk in range(1,mod.nk):
                jj = 2*(kk-1)
                md.coeffs[jj] = 0.
                md.coeffs[jj+1] = 0.
            ind_start = 1
        else:
            ind_start = 0

        for ii in range(ind_start, plan.n_data-1):
            xjp1oxj = plan.grid[ii+1]/plan.grid[ii]
            for kk in range(1,mod.nk):
                jj = 2*(mod.nk-1)*ii + 2*(kk-1)
                md.coeffs[jj] = xjp1oxj**mod.lamk[kk]
                md.coeffs[jj+1] = constants.piinv*mod.hk[kk]*(md.coeffs[jj]-1.)/mod.lamk[kk]
    else:
        destroy_fat_hansen_law_linear(plan)
        with gil:
            raise NotImplementedError('Method not implemented for given parameters.')

    return 0


# Hansen Law with linear approximation of function
cdef int execute_fat_hansen_law_linear(abel_plan* plan, double* data_in, double* data_out) except -1 nogil:

    cdef:
        method_data_hansen_law* md = <method_data_hansen_law*> plan.method_data
        model_hansen_law* mod = md.model
        int ii, jj, kk, n_data = plan.n_data, method = plan.method
        (double*) xk, sn = [0., 0.]
        double shift = plan.shift, data_in_old

    xk = <double*> malloc(mod.nk*sizeof(double))

    if plan.forward_backward == -1:
        data_in_old = data_in[n_data-1]
        data_out[n_data-1] = 0.
        for kk in range(mod.nk):
            xk[kk] = 0.
        for ii in range(n_data-2, -1, -1):
            sn[1] = (data_in_old-data_in[ii])/plan.step_size
            sn[0] = data_in[ii] - plan.grid[ii]*sn[1]
            xk[0] += (2.*sn[0] + (plan.grid[ii+1]+plan.grid[ii])*sn[1])*plan.step_size
            data_in_old = data_in[ii]
            data_out[ii] = xk[0]
            for kk in range(1, mod.nk):
                jj = 3*(mod.nk-1)*ii + 3*(kk-1)
                xk[kk] = xk[kk]*md.coeffs[jj] - md.coeffs[jj+1]*sn[0] - md.coeffs[jj+2]*sn[1]
                data_out[ii] += xk[kk] 
    elif plan.forward_backward == 1:
        data_in_old = data_in[n_data-1]
        data_out[n_data-1] = 0.
        for kk in range(mod.nk):
            xk[kk] = 0.
        for ii in range(n_data-2, 0, -1):
            sn[1] = (data_in_old-data_in[ii])/plan.step_size
            xk[0] += -constants.piinv*log(plan.grid[ii+1]/plan.grid[ii])*sn[1]
            data_in_old = data_in[ii]
            data_out[ii] = xk[0]
            for kk in range(1, mod.nk):
                jj = 2*(mod.nk-1)*ii + 2*(kk-1)
                xk[kk] = xk[kk]*md.coeffs[jj] - md.coeffs[jj+1]*sn[1]
                data_out[ii] += xk[kk]
        data_out[0] = data_out[1]
    elif plan.forward_backward == 2:
        data_in_old = data_in[n_data-1]
        data_out[n_data-1] = 0.
        for kk in range(mod.nk):
            xk[kk] = 0.
        for ii in range(n_data-2, 0, -1):
            sn[1] = (data_in[ii]+data_in[ii-1])*0.5
            xk[0] += -constants.piinv*log(plan.grid[ii+1]/plan.grid[ii])*sn[1]
            data_in_old = data_in[ii]
            data_out[ii] = xk[0]
            for kk in range(1, mod.nk):
                jj = 2*(mod.nk-1)*ii + 2*(kk-1)
                xk[kk] = xk[kk]*md.coeffs[jj] - md.coeffs[jj+1]*sn[1]
                data_out[ii] += xk[kk]
        data_out[0] = data_out[1]
    else:
        free(xk)
        with gil:
            raise NotImplementedError('Method not implemented for given parameters.')

    free(xk)

    return 0


cdef int destroy_fat_hansen_law_linear(abel_plan* plan) except -1 nogil:

    cdef:
        method_data_hansen_law* md = <method_data_hansen_law*> plan.method_data

    free(md.coeffs)
    free(md)

    return 0
