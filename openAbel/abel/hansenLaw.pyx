

from libc.stdlib cimport malloc, free


cimport openAbel.constants as constants
cimport openAbel.mathFun as mathFun
from openAbel.abel.base cimport abel_plan






############################################################################################################################################
### Hansen Law space state model for Abel transform                                                                                      ###


ctypedef struct methodData_hansenLaw:
    double* coeffs
    model_hansenLaw* model


ctypedef struct model_hansenLaw:
    unsigned int nk
    double* hk
    double* lamk


# Original model (factor pi different definition than original Hansen Law paper)
cdef model_hansenLaw model_hansenLawOrg
model_hansenLawOrg.nk = 9
model_hansenLawOrg.hk = [1., 0.596903, 1.09956, 2.57611, 5.65487, 12.2522, 26.0752, 61.5752, 151.739]
model_hansenLawOrg.lamk = [0.0, -2.1, -6.2, -22.4, -92.5, -414.5, -1889.4, -8990.9, -47391.1]



############################################################################################################################################
# Plan Hansen Law original model (9th order linear)
cdef int plan_fat_hansenLawOrgLin(abel_plan* plan) nogil:
    
    cdef:
        methodData_hansenLaw* methodData = <methodData_hansenLaw*> malloc(sizeof(methodData_hansenLaw))

    plan.methodData = <void*> methodData
    methodData.model = &model_hansenLawOrg
    methodData.coeffs = NULL

    return plan_fat_hansenLawLinear(plan)

############################################################################################################################################
# Plan Hansen Law linear
cdef int plan_fat_hansenLawLinear(abel_plan* plan) nogil:
    
    cdef:
        unsigned int ii, jj, kk, indStart
        methodData_hansenLaw* methodData = <methodData_hansenLaw*> plan.methodData
        double xjp1oxj

    if plan.forwardBackward == -1:

        methodData.coeffs = <double*> malloc(3*plan.nData*(methodData.model.nk-1)*sizeof(double))
        if NULL == methodData.coeffs:
            free(methodData)
            with gil:
                raise MemoryError('Malloc ruturned a NULL pointer, probably not enough memory available.')

        if plan.shift == 0.:
            for kk in range(1,methodData.model.nk):
                jj = 3*(kk-1)
                methodData.coeffs[jj] = 0.
                methodData.coeffs[jj+1] = 0.
                methodData.coeffs[jj+2] = 0.
            indStart = 1
        else:
            indStart = 0
        for ii in range(indStart, plan.nData-1):
            xjp1oxj = plan.grid[ii+1]/plan.grid[ii]
            for kk in range(1,methodData.model.nk):
                jj = 3*(methodData.model.nk-1)*ii + 3*(kk-1)
                methodData.coeffs[jj] = xjp1oxj**methodData.model.lamk[kk]
                methodData.coeffs[jj+1] = 2.*methodData.model.hk[kk]*(1.-methodData.coeffs[jj]*xjp1oxj)*plan.grid[ii] / \
                                          (methodData.model.lamk[kk]+1.)
                methodData.coeffs[jj+2] = 2.*methodData.model.hk[kk]*(1.-methodData.coeffs[jj]*xjp1oxj**2)*plan.grid[ii]**2 / \
                                          (methodData.model.lamk[kk]+2.)

    elif plan.forwardBackward == 1 or plan.forwardBackward == 2:

        methodData.coeffs = <double*> malloc(2*plan.nData*(methodData.model.nk-1)*sizeof(double))
        if NULL == methodData.coeffs:
            free(methodData)
            with gil:
                raise MemoryError('Malloc ruturned a NULL pointer, probably not enough memory available.')

        if plan.shift == 0.:
            for kk in range(1,methodData.model.nk):
                jj = 2*(kk-1)
                methodData.coeffs[jj] = 0.
                methodData.coeffs[jj+1] = 0.
            indStart = 1
        else:
            indStart = 0

        for ii in range(indStart, plan.nData-1):
            xjp1oxj = plan.grid[ii+1]/plan.grid[ii]
            for kk in range(1,methodData.model.nk):
                jj = 2*(methodData.model.nk-1)*ii + 2*(kk-1)
                methodData.coeffs[jj] = xjp1oxj**methodData.model.lamk[kk]
                methodData.coeffs[jj+1] = constants.piinv*methodData.model.hk[kk]*(methodData.coeffs[jj]-1.)/methodData.model.lamk[kk]
    else:
        with gil:
            raise NotImplementedError

    return 0


# Hansen Law with linear approximation of function
cdef int execute_fat_hansenLawLinear(abel_plan* plan, double* dataIn, double* dataOut) nogil:

    cdef:
        int ii, jj, kk
        double* sn = [0., 0.]
        double dataInOld
        double* xk
        methodData_hansenLaw* methodData = <methodData_hansenLaw*> plan.methodData
        model_hansenLaw* model = methodData.model
        int nData = plan.nData
        double shift = plan.shift
        int method = plan.method

    xk = <double*> malloc(model.nk*sizeof(double))
    if NULL == xk:
        free(xk)
        with gil:
            raise MemoryError('Malloc ruturned a NULL pointer, probably not enough memory available.')

    if plan.forwardBackward == -1:
        dataInOld = dataIn[nData-1]
        dataOut[nData-1] = 0.
        for kk in range(model.nk):
            xk[kk] = 0.
        for ii in range(nData-2, -1, -1):
            sn[1] = (dataInOld-dataIn[ii])/plan.stepSize
            sn[0] = dataIn[ii] - plan.grid[ii]*sn[1]
            xk[0] += (2.*sn[0] + (plan.grid[ii+1]+plan.grid[ii])*sn[1])*plan.stepSize
            dataInOld = dataIn[ii]
            dataOut[ii] = xk[0]
            for kk in range(1, model.nk):
                jj = 3*(model.nk-1)*ii + 3*(kk-1)
                xk[kk] = xk[kk]*methodData.coeffs[jj] - methodData.coeffs[jj+1]*sn[0] - methodData.coeffs[jj+2]*sn[1]
                dataOut[ii] += xk[kk] 
    elif plan.forwardBackward == 1:
        dataInOld = dataIn[nData-1]
        dataOut[nData-1] = 0.
        for kk in range(model.nk):
            xk[kk] = 0.
        for ii in range(nData-2, 0, -1):
            sn[1] = (dataInOld-dataIn[ii])/plan.stepSize
            xk[0] += -constants.piinv*mathFun.log(plan.grid[ii+1]/plan.grid[ii])*sn[1]
            dataInOld = dataIn[ii]
            dataOut[ii] = xk[0]
            for kk in range(1, model.nk):
                jj = 2*(model.nk-1)*ii + 2*(kk-1)
                xk[kk] = xk[kk]*methodData.coeffs[jj] - methodData.coeffs[jj+1]*sn[1]
                dataOut[ii] += xk[kk]
        dataOut[0] = dataOut[1]
    elif plan.forwardBackward == 2:
        dataInOld = dataIn[nData-1]
        dataOut[nData-1] = 0.
        for kk in range(model.nk):
            xk[kk] = 0.
        for ii in range(nData-2, 0, -1):
            sn[1] = (dataIn[ii]+dataIn[ii-1])*0.5
            xk[0] += -constants.piinv*mathFun.log(plan.grid[ii+1]/plan.grid[ii])*sn[1]
            dataInOld = dataIn[ii]
            dataOut[ii] = xk[0]
            for kk in range(1, model.nk):
                jj = 2*(model.nk-1)*ii + 2*(kk-1)
                xk[kk] = xk[kk]*methodData.coeffs[jj] - methodData.coeffs[jj+1]*sn[1]
                dataOut[ii] += xk[kk]
        dataOut[0] = dataOut[1]
    else:
        with gil:
            raise NotImplementedError

    free(xk)

    return 0


cdef int destroy_fat_hansenLawLinear(abel_plan* plan) nogil:

    cdef:
        methodData_hansenLaw* methodData = <methodData_hansenLaw*> plan.methodData

    free(methodData.coeffs)
    free(methodData)

    return 0

