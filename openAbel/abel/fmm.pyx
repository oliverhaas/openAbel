

import numpy as np
import os.path as osp
import datetime

from libc.stdlib cimport free
from openAbel.helper cimport nullCheckMalloc as malloc, nullCheckCalloc as calloc
from libc.string cimport memset

cimport scipy.linalg.cython_blas as blas

import openAbel.abel.coeffs as cffs

cimport openAbel.mathFun as mf
cimport openAbel.constants as co
from openAbel.abel.base cimport abel_plan


########################################################################################################################
### Fast multipole method trapezoidal with end corrections                                                           ###


ctypedef struct methodData_FMM:
    (int*) kl, klCum
    (double*) chebRoots, mtmp, mtmm, ltp, mtlk, direct, direct0, coeffsSing, coeffsNonsing, coeffsFilter
    int pp, pp1, ss, nlevs, kTotal, order, orderFilter

cdef double _kernForward(double rr, double yy) nogil:
    return rr/mf.sqrt(rr**2-yy**2)
    
cdef double _kernBackward(double rr, double yy) nogil:
    return 1/mf.sqrt(rr**2-yy**2)

cdef double _kernModified(double rr, double yy) nogil:
    return (yy/rr)**2/mf.sqrt(rr**2-yy**2)

# Plan FMM
cdef int plan_fat_fmmTrapEndCorr(abel_plan* pl, int order = 2, double eps = co.machineEpsilon) nogil except -1:

    cdef:
        int ii, jj, ll, kk, mm
        methodData_FMM* md
        double ti, tauj, temp, nInvSca, yInvSca
        double[:,::1] cffs_s_sm_mv, cffs_s_la_mv, cffs_ns_sqrt_sm_mv, cffs_ns_sqrt_la_mv
        double[::1] cffs_f_mv
        int ordM1Hlf, ordM1HlfIn, yCross, yLarge, yInvScaInt, nCross, nLarge, nInvScaInt
        double (*kern)(double, double) nogil

    # Input check
    if NULL == pl or order <= 0 or eps < co.machineEpsilon:
        with gil:
            raise ValueError('Illegal input argument.')   

    # Main method struct
    md = <methodData_FMM*> malloc(sizeof(methodData_FMM))
    # Initialize to NULL so I can destroy properly when exception is raised
    md.chebRoots = md.kl = md.klCum = md.mtmp = md.mtmm = md.mtlk = NULL 
    md.direct = md.coeffsSing = md.coeffsNonsing = md.coeffsFilter = NULL
    pl.methodData = <void*> md

    # Small data set
    if pl.nData < order+2:
        destroy_fat_fmmTrapEndCorr(pl)
        with gil:
            raise ValueError('Not enough data points for given parameters.')

    # Load and prepare end correction coefficients
    md.order = order
    ordM1Hlf = <int> ((md.order-1)/2)
    ordM1HlfIn = <int> (md.order/2)
    md.coeffsSing = <double*> malloc(md.order*(pl.nData-1)*sizeof(double))
    md.coeffsNonsing = <double*> malloc(md.order*(pl.nData-1)*sizeof(double))
    if pl.forwardBackward == -1:    # Forward transform
        kern = &_kernForward
        with gil:
            try:
                cffs_s_la_mv = cffs.getCoeffs('coeffs_invSqrtDiffSqLin_sing_large', order)
                cffs_ns_sqrt_sm_mv = cffs.getCoeffs('coeffs_invSqrt_nonsing_small', order)
                cffs_ns_sqrt_la_mv = cffs.getCoeffs('coeffs_invSqrt_nonsing_large', order)
                if pl.shift == 0.:
                    cffs_s_sm_mv = cffs.getCoeffs('coeffs_invSqrtDiffSqLin_sing_small', order)
                elif pl.shift == 0.5:
                    cffs_s_sm_mv = cffs.getCoeffs('coeffs_invSqrtDiffSqLin_sing_small_halfShift', order)
                else:
                    raise NotImplementedError('Method not implemented for given parameters.')
            except:
                destroy_fat_fmmTrapEndCorr(pl)
                raise
        yCross = cffs_s_sm_mv.shape[0]
        yLarge = cffs_s_la_mv.shape[0]
        for ii in range(min(yCross,pl.nData-1)):
            for jj in range(md.order):
                md.coeffsSing[md.order*ii+jj] = cffs_s_sm_mv[ii,jj]
        for ii in range(yCross, pl.nData-1):
            yInvSca = pl.stepSize/pl.grid[ii]*(yCross-1)*(yLarge-1)
            yInvScaInt = <int> mf.fmax(mf.fmin(yInvSca,yLarge-3),1)
            for jj in range(md.order):
                md.coeffsSing[md.order*ii+jj] = interpCubic(yInvSca-yInvScaInt, md.order,
                                                            &cffs_s_la_mv[yInvScaInt-1,jj]) * \
                                                mf.sqrt(pl.grid[ii]/2./pl.stepSize)
        nCross = cffs_ns_sqrt_sm_mv.shape[0]            
        for ii in range(max(pl.nData-1-nCross,0),pl.nData-1):
            for jj in range(md.order):
                md.coeffsNonsing[md.order*ii+jj] = cffs_ns_sqrt_sm_mv[pl.nData-2-ii,jj] * \
                                                   (pl.grid[pl.nData-1]+(jj-ordM1HlfIn)*pl.stepSize) / \
                                                   mf.sqrt((pl.grid[pl.nData-1]+(jj-ordM1HlfIn)*pl.stepSize+pl.grid[ii])*(pl.grid[pl.nData-1]-pl.grid[ii]))
        nLarge = cffs_ns_sqrt_la_mv.shape[0]      
        for ii in range(max(pl.nData-1-nCross,0)):
            nInvSca = pl.stepSize/(pl.grid[pl.nData-1]-pl.grid[ii])*nCross*(nLarge-1)
            nInvScaInt = <int> mf.fmax(mf.fmin(nInvSca,nLarge-3),1)
            for jj in range(md.order):
                md.coeffsNonsing[md.order*ii+jj] = interpCubic(nInvSca-nInvScaInt, md.order, 
                                                               &cffs_ns_sqrt_la_mv[nInvScaInt-1,jj]) * \
                                                   (pl.grid[pl.nData-1]+(jj-ordM1HlfIn)*pl.stepSize) / \
                                                   mf.sqrt((pl.grid[pl.nData-1] + \
                                                                 (jj-ordM1HlfIn)*pl.stepSize+pl.grid[ii]) * \
                                                                (pl.grid[pl.nData-1]-pl.grid[ii]))
        for ii in range(pl.nData-1):
            md.coeffsNonsing[md.order*ii+ordM1HlfIn] -= 0.5*kern(pl.grid[pl.nData-1], pl.grid[ii])

    elif pl.forwardBackward == 1 or pl.forwardBackward == 2:    # Backward transform
        kern = &_kernBackward
        with gil:
            try:
                cffs_s_la_mv = cffs.getCoeffs('coeffs_invSqrtDiffSq_sing_large', order)
                cffs_ns_sqrt_sm_mv = cffs.getCoeffs('coeffs_invSqrt_nonsing_small', order)
                cffs_ns_sqrt_la_mv = cffs.getCoeffs('coeffs_invSqrt_nonsing_large', order)
                if pl.shift == 0.:
                    cffs_s_sm_mv = cffs.getCoeffs('coeffs_invSqrtDiffSq_sing_small', order)
                elif pl.shift == 0.5:
                    cffs_s_sm_mv = cffs.getCoeffs('coeffs_invSqrtDiffSq_sing_small_halfShift', order)
                else:
                    raise NotImplementedError('Method not implemented for given parameters.')
            except:
                destroy_fat_fmmTrapEndCorr(pl)
                raise
        yCross = cffs_s_sm_mv.shape[0]
        yLarge = cffs_s_la_mv.shape[0]
        for ii in range(min(yCross,pl.nData-1)):
            for jj in range(md.order):
                md.coeffsSing[md.order*ii+jj] = cffs_s_sm_mv[ii,jj]/pl.stepSize
        for ii in range(yCross, pl.nData-1):
            yInvSca = pl.stepSize/pl.grid[ii]*(yCross-1)*(yLarge-1)
            yInvScaInt = <int> mf.fmax(mf.fmin(yInvSca,yLarge-3),1)
            for jj in range(md.order):
                md.coeffsSing[md.order*ii+jj] = interpCubic(yInvSca-yInvScaInt, md.order, 
                                                            &cffs_s_la_mv[yInvScaInt-1,jj]) / \
                                                mf.sqrt(pl.grid[ii]*2.*pl.stepSize)
        nCross = cffs_ns_sqrt_sm_mv.shape[0]            
        for ii in range(max(pl.nData-1-nCross,0),pl.nData-1):
            for jj in range(md.order):
                md.coeffsNonsing[md.order*ii+jj] = cffs_ns_sqrt_sm_mv[pl.nData-2-ii,jj] / \
                                                   mf.sqrt((pl.grid[pl.nData-1] + \
                                                            (jj-ordM1HlfIn)*pl.stepSize+pl.grid[ii]) * \
                                                           (pl.grid[pl.nData-1]-pl.grid[ii]))
        nLarge = cffs_ns_sqrt_la_mv.shape[0]      
        for ii in range(max(pl.nData-1-nCross,0)):
            nInvSca = pl.stepSize/(pl.grid[pl.nData-1]-pl.grid[ii])*nCross*(nLarge-1)
            nInvScaInt = <int> mf.fmax(mf.fmin(nInvSca,nLarge-3),1)
            for jj in range(md.order):
                md.coeffsNonsing[md.order*ii+jj] = interpCubic(nInvSca-nInvScaInt, md.order, 
                                                               &cffs_ns_sqrt_la_mv[nInvScaInt-1,jj]) / \
                                                   mf.sqrt((pl.grid[pl.nData-1] + \
                                                                 (jj-ordM1HlfIn)*pl.stepSize+pl.grid[ii]) * \
                                                                (pl.grid[pl.nData-1]-pl.grid[ii]))
        for ii in range(pl.nData-1):
            md.coeffsNonsing[md.order*ii+ordM1HlfIn] -= 0.5*kern(pl.grid[pl.nData-1],pl.grid[ii])

    elif pl.forwardBackward == -2:    # Modified forward transform for 1/r^2 singular functions
        kern = &_kernModified
        with gil:
            try:
                cffs_s_la_mv = cffs.getCoeffs('coeffs_invSqrtDiffSqY2OR2_sing_large', order)
                cffs_ns_sqrt_sm_mv = cffs.getCoeffs('coeffs_invSqrt_nonsing_small', order)
                cffs_ns_sqrt_la_mv = cffs.getCoeffs('coeffs_invSqrt_nonsing_large', order)
                if pl.shift == 0.:
                    cffs_s_sm_mv = cffs.getCoeffs('coeffs_invSqrtDiffSqY2OR2_sing_small', order)
                elif pl.shift == 0.5:
                    cffs_s_sm_mv = cffs.getCoeffs('coeffs_invSqrtDiffSqY2OR2_sing_small_halfShift_', order)
                else:
                    raise NotImplementedError('Method not implemented for given parameters.')
            except:
                destroy_fat_fmmTrapEndCorr(pl)
                raise
        yCross = cffs_s_sm_mv.shape[0]
        yLarge = cffs_s_la_mv.shape[0]
        for ii in range(min(yCross,pl.nData-1)):
            for jj in range(md.order):
                md.coeffsSing[md.order*ii+jj] = cffs_s_sm_mv[ii,jj]/pl.stepSize
        for ii in range(yCross, pl.nData-1):
            yInvSca = pl.stepSize/pl.grid[ii]*(yCross-1)*(yLarge-1)
            yInvScaInt = <int> mf.fmax(mf.fmin(yInvSca,yLarge-3),1)
            for jj in range(md.order):
                md.coeffsSing[md.order*ii+jj] = interpCubic(yInvSca-yInvScaInt, md.order,
                                                            &cffs_s_la_mv[yInvScaInt-1,jj]) / \
                                                mf.sqrt(pl.grid[ii]*2.*pl.stepSize)
        nCross = cffs_ns_sqrt_sm_mv.shape[0]            
        for ii in range(max(pl.nData-1-nCross,0),pl.nData-1):
            for jj in range(md.order):
                md.coeffsNonsing[md.order*ii+jj] = cffs_ns_sqrt_sm_mv[pl.nData-2-ii,jj] * \
                                                   (pl.grid[ii]/(pl.grid[pl.nData-1]+(jj-ordM1HlfIn)*pl.stepSize))**2 / \
                                                   mf.sqrt( (pl.grid[pl.nData-1]+(jj-ordM1HlfIn)*pl.stepSize+pl.grid[ii]) *
                                                                 (pl.grid[pl.nData-1]-pl.grid[ii]) )
        nLarge = cffs_ns_sqrt_la_mv.shape[0]      
        for ii in range(max(pl.nData-1-nCross,0)):
            nInvSca = pl.stepSize/(pl.grid[pl.nData-1]-pl.grid[ii])*nCross*(nLarge-1)
            nInvScaInt = <int> mf.fmax(mf.fmin(nInvSca,nLarge-3),1)
            for jj in range(md.order):
                md.coeffsNonsing[md.order*ii+jj] = interpCubic(nInvSca-nInvScaInt, md.order, 
                                                               &cffs_ns_sqrt_la_mv[nInvScaInt-1,jj]) * \
                                                   (pl.grid[ii]/(pl.grid[pl.nData-1]+(jj-ordM1HlfIn)*pl.stepSize))**2 / \
                                                   mf.sqrt( (pl.grid[pl.nData-1]+(jj-ordM1HlfIn)*pl.stepSize+pl.grid[ii]) *
                                                                 (pl.grid[pl.nData-1]-pl.grid[ii]) )
        for ii in range(pl.nData-1):
            md.coeffsNonsing[md.order*ii+ordM1HlfIn] -= 0.5*kern(pl.grid[pl.nData-1],pl.grid[ii])

    else:
        destroy_fat_fmmTrapEndCorr(pl)
        with gil:
            raise NotImplementedError('Method not implemented for given parameters.')

    # Hierarchical decomposition
    md.pp = max(4, <int> ( -0.55*mf.log(2.*eps) + 1. ) )    # Empirical scaling, should be exponential
    md.pp1 = md.pp + 1
    md.nlevs = max(<int> ( mf.log2((pl.nData-1.)/(2.*md.pp)) + 1. ), 2)
    md.ss = max(<int> ( (pl.nData-1.)/2**md.nlevs + 1. ), 3)    # ss ~= 2*pp theoretical
    md.kTotal = 2**(md.nlevs+1) - 1                             # Total number of intervals in all levels
    
    # Allocation of arrays FMM part
    md.chebRoots = <double*> malloc(md.pp1*sizeof(double))
    md.kl = <int*> malloc((md.nlevs+1)*sizeof(int))
    md.klCum = <int*> malloc((md.nlevs+1)*sizeof(int))
    md.mtmp = <double*> malloc(md.pp1**2*sizeof(double))
    md.mtmm = <double*> malloc(md.pp1**2*sizeof(double))
    md.ltp = <double*> malloc(md.pp1*md.ss*sizeof(double))
    md.mtlk = <double*> calloc(2*md.kTotal*md.pp1**2, sizeof(double))
    md.direct = <double*> calloc(2**md.nlevs*md.ss**2*2, sizeof(double))
    md.direct0 = <double*> malloc(2*md.ss*sizeof(double))

    # More hierarchical decomposition stuff
    md.klCum[0] = 0
    md.kl[0] = 2**md.nlevs
    for ii in range(1,md.nlevs+1):
        md.kl[ii] = 2**(md.nlevs-ii)
        md.klCum[ii] = md.klCum[ii-1] + md.kl[ii-1]
    _chebRoots(md.pp1, md.chebRoots)

    # Moment to moment coefficients
    for ii in range(md.pp1):
        for jj in range(md.pp1):
            md.mtmp[ii*md.pp1+jj] = _lagrangePInt(0.5*md.chebRoots[jj]+0.5, ii, md.chebRoots, md.pp1)
            md.mtmm[ii*md.pp1+jj] = _lagrangePInt(0.5*md.chebRoots[jj]-0.5, ii, md.chebRoots, md.pp1)    

    # Local to potential coefficients
    for ii in range(md.ss):
        temp = 2.*(ii+1)/md.ss-1.
        for jj in range(md.pp1):
            md.ltp[ii*md.pp1+jj] = _lagrangePInt(temp, jj, md.chebRoots, md.pp1)

    # Moment to local coefficients
    for ll in range(md.nlevs-1):
        # If even
        for kk in range(0, md.kl[ll]-2, 2):
            for ii in range(md.pp1):
                ti = (md.kl[md.nlevs-ll]*(kk+0.5+0.5*md.chebRoots[ii])*md.ss + pl.shift)*pl.stepSize
                for jj in range(md.pp1):
                    tauj = (md.kl[md.nlevs-ll]*(kk+2.5+0.5*md.chebRoots[jj])*md.ss + pl.shift)*pl.stepSize
                    md.mtlk[(2*(md.klCum[ll]+kk))*md.pp1**2+ii*md.pp1+jj] = kern(tauj, ti)
                    tauj = (md.kl[md.nlevs-ll]*(kk+3.5+0.5*md.chebRoots[jj])*md.ss + pl.shift)*pl.stepSize
                    md.mtlk[(2*(md.klCum[ll]+kk)+1)*md.pp1**2+ii*md.pp1+jj] = kern(tauj, ti)
        # If odd
        for kk in range(1, md.kl[ll]-2, 2):
            for ii in range(md.pp1):
                ti = (md.kl[md.nlevs-ll]*(kk+0.5+0.5*md.chebRoots[ii])*md.ss + pl.shift)*pl.stepSize
                for jj in range(md.pp1):
                    tauj = (md.kl[md.nlevs-ll]*(kk+2.5+0.5*md.chebRoots[jj])*md.ss + pl.shift)*pl.stepSize
                    md.mtlk[(2*(md.klCum[ll]+kk))*md.pp1**2+ii*md.pp1+jj] = kern(tauj, ti)

    # Direct short range coefficients
    for ii in range(1, pl.nData):
        kk = <int> ((ii-1)/md.ss)
        ll = (ii-1) - kk*md.ss
        mm = min(pl.nData-kk*md.ss-1, 2*md.ss)
        for jj in range(ll+1, mm):
            md.direct[kk*md.ss**2*2+md.ss*2*ll+jj] = kern(pl.grid[ii+jj-ll], pl.grid[ii])
    mm = min(pl.nData, 2*md.ss+1)
    for jj in range(1, mm):
        md.direct0[jj-1] = kern(pl.grid[jj], pl.grid[0])    

    # Input modification filter
    if pl.forwardBackward == 1:
        md.orderFilter = md.order+1 + (md.order % 2)
        md.coeffsFilter = <double*> malloc(md.orderFilter*sizeof(double))
        with gil:
            try:
                cffs_f_mv = cffs.getCoeffs('coeffs_deriv_smooth', md.orderFilter-1)
            except:
                destroy_fat_fmmTrapEndCorr(pl)
                raise
        for ii in range(md.orderFilter):
            md.coeffsFilter[ii] = cffs_f_mv[ii]*(-co.piinv)
    elif pl.forwardBackward == 2 or pl.forwardBackward == -1 or pl.forwardBackward == -2:
        md.orderFilter = 1
        md.coeffsFilter = <double*> malloc(1*sizeof(double))
        if pl.forwardBackward == 2:            
            md.coeffsFilter[0] = -co.piinv*pl.stepSize
        else:
            md.coeffsFilter[0] = 2.*pl.stepSize

    return 0


# Execute FMM
cdef int execute_fat_fmmTrapEndCorr(abel_plan* pl, double* dataIn, double* dataOut, int leftBoundary, 
                                    int rightBoundary) nogil except -1:

    cdef:
        methodData_FMM* md
        int ii, jj, kk, ll, mm, nn
        (double*) moments = NULL, local = NULL, dataInTemp0 = NULL, dataInTemp1 = NULL
        int ordM1Hlf, ordFilM1Hlf, ordM1HlfIn, nLeftExt, nRightExt

    md = <methodData_FMM*> pl.methodData
    ordM1Hlf = <int> ((md.order-1)/2)
    ordM1HlfIn = <int> (md.order/2)
    ordFilM1Hlf = <int> ((md.orderFilter-1)/2)

    # Allocate temporary data arrays
    dataInTemp0 = <double*> malloc((pl.nData+md.order+md.orderFilter-2)*sizeof(double))
    dataInTemp1 = <double*> malloc((pl.nData+md.order-1)*sizeof(double))
    
    # Left boundary handling
    if leftBoundary == 0 or leftBoundary == 1 or leftBoundary == 2:
        nLeftExt = ordM1Hlf + ordFilM1Hlf
    elif leftBoundary == 3:
        nLeftExt = 0
    else:
        free(dataInTemp0)
        free(dataInTemp1)
        with gil:
            raise NotImplementedError('Method not implemented for given parameters.')
    # Right boundary handling
    if rightBoundary == 0: # TODO maybe or rightBoundary == 1 or rightBoundary == 2:
        nRightExt = ordM1Hlf + ordFilM1Hlf
    elif rightBoundary == 3:
        nRightExt = 0
    else:
        free(dataInTemp0)
        free(dataInTemp1)
        with gil:
            raise NotImplementedError('Method not implemented for given parameters.')           
    # Copy and extend data if necessary
    nn = max(md.order, md.orderFilter-1)
    for ii in range(pl.nData+md.order+md.orderFilter-2-nLeftExt-nRightExt):
        dataInTemp0[nLeftExt+ii] = dataIn[ii]
    if leftBoundary == 0:
        for ii in range(nLeftExt):
            dataInTemp0[ii] = polint(&dataInTemp0[nLeftExt], nn, ii-nLeftExt)
    elif leftBoundary == 1:
        if pl.shift == 0.:
            for ii in range(nLeftExt):
                dataInTemp0[nLeftExt-1-ii] = -dataInTemp0[nLeftExt+1+ii]
        elif pl.shift == 0.5:
            for ii in range(nLeftExt):
                dataInTemp0[nLeftExt-1-ii] = -dataInTemp0[nLeftExt+ii]
        else:
            free(dataInTemp0)
            free(dataInTemp1)
            with gil:
                raise NotImplementedError('Method not implemented for given parameters.')
    elif leftBoundary == 2:
        if pl.shift == 0.:
            for ii in range(nLeftExt):
                dataInTemp0[nLeftExt-1-ii] = dataInTemp0[nLeftExt+1+ii]
        elif pl.shift == 0.5:
            for ii in range(nLeftExt):
                dataInTemp0[nLeftExt-1-ii] = dataInTemp0[nLeftExt+ii]
        else:
            free(dataInTemp0)
            free(dataInTemp1)
            with gil:
                raise NotImplementedError('Method not implemented for given parameters.')
    elif leftBoundary == 3:
        pass
    else:
        free(dataInTemp0)
        free(dataInTemp1)
        with gil:
            raise NotImplementedError('Method not implemented for given parameters.')
    if rightBoundary == 0:
        for ii in range(nRightExt):
            dataInTemp0[pl.nData+ordM1Hlf+ordFilM1Hlf+ii] = polint(&dataInTemp0[pl.nData+ordM1Hlf+ordFilM1Hlf-nn], nn, ii+nn)
    elif rightBoundary == 3:
        pass
    else:
        free(dataInTemp0)
        free(dataInTemp1)
        with gil:
            raise NotImplementedError('Method not implemented for given parameters.')

    # Do scaling or numerical derivative
    convolve(dataInTemp0, pl.nData+md.order-1, dataInTemp1, md.orderFilter, md.coeffsFilter)
    free(dataInTemp0)

    # Allocate temporary data arrays
    moments = <double*> calloc(md.kTotal*md.pp1, sizeof(double))
    local = <double*> calloc(md.kTotal*md.pp1, sizeof(double))

    # Finest level moment calculation
    mm = (pl.nData-2)/md.ss - 2         # basically kl[0]-2, as first two block are not needed
    blas.dgemm('n', 'n', &md.pp1, &mm, &md.ss, &ONED, md.ltp, &md.pp1, 
               &dataInTemp1[ordM1Hlf+2*md.ss+1], &md.ss, &ZEROD, &moments[2*md.pp1], &md.pp1)
    mm += 2
    for ii in range(mm*md.ss+1, pl.nData):
        nn = (ii-1) - mm*md.ss
        for jj in range(md.pp1):
            moments[mm*md.pp1+jj] += md.ltp[nn*md.pp1+jj]*dataInTemp1[ordM1Hlf+ii]

    # Upward Pass / Moment to moment
    mm = 2*md.pp1
    for ll in range(1, md.nlevs-1):
        blas.dgemm('t', 'n', &md.pp1, &md.kl[ll], &md.pp1, &ONED, md.mtmm, &md.pp1, 
                   &moments[md.klCum[ll-1]*md.pp1], &mm, &ONED, &moments[md.klCum[ll]*md.pp1], &md.pp1)
        blas.dgemm('t', 'n', &md.pp1, &md.kl[ll], &md.pp1, &ONED, md.mtmp, &md.pp1, 
                   &moments[(md.klCum[ll-1]+1)*md.pp1], &mm, &ONED, &moments[md.klCum[ll]*md.pp1], &md.pp1)

    # Interaction Phase / Moment to local
    for ll in range(md.nlevs-1):
        # If even
        for kk in range(0, md.kl[ll]-2, 2):
            blas.dgemv('t', &md.pp1, &md.pp1, &ONED, &md.mtlk[(2*(md.klCum[ll]+kk))*md.pp1**2], &md.pp1, 
                       &moments[(kk+md.klCum[ll]+2)*md.pp1], &ONE, &ZEROD, &local[(md.klCum[ll]+kk)*md.pp1], &ONE)
            blas.dgemv('t', &md.pp1, &md.pp1, &ONED, &md.mtlk[(2*(md.klCum[ll]+kk)+1)*md.pp1**2], &md.pp1, 
                       &moments[(kk+md.klCum[ll]+3)*md.pp1], &ONE, &ONED, &local[(md.klCum[ll]+kk)*md.pp1], &ONE)
        # If odd
        for kk in range(1, md.kl[ll]-2, 2):
            blas.dgemv('t', &md.pp1, &md.pp1, &ONED, &md.mtlk[(2*(md.klCum[ll]+kk))*md.pp1**2], &md.pp1, 
                       &moments[(kk+md.klCum[ll]+2)*md.pp1], &ONE, &ZEROD, &local[(md.klCum[ll]+kk)*md.pp1], &ONE)
    
    # Downward Pass / Local to local
    mm = 2*md.pp1
    for ll in range(md.nlevs-2, 0, -1):
        blas.dgemm('n', 'n', &md.pp1, &md.kl[ll], &md.pp1, &ONED, md.mtmm, &md.pp1, 
                   &local[md.klCum[ll]*md.pp1], &md.pp1, &ONED, &local[md.klCum[ll-1]*md.pp1], &mm)
        blas.dgemm('n', 'n', &md.pp1, &md.kl[ll], &md.pp1, &ONED, md.mtmp, &md.pp1, 
                   &local[md.klCum[ll]*md.pp1], &md.pp1, &ONED, &local[(md.klCum[ll-1]+1)*md.pp1], &mm)

    # Set output to zero first
    memset(dataOut, 0, pl.nData*sizeof(double))

    # Potential evaluation / local to potential
    mm = (pl.nData-2)/md.ss - 1          # basically kl[0]-2, as last two block are not needed
    blas.dgemm('t', 'n', &md.ss, &mm, &md.pp1, &ONED, md.ltp, &md.pp1, local, &md.pp1, &ZEROD, &dataOut[1], &md.ss)
    for jj in range(md.pp1):
        dataOut[0] += local[jj]*_lagrangePInt(-1., jj, md.chebRoots, md.pp1)
    
    free(moments)
    free(local)

    # Direct short range
    ll = min((pl.nData-2)/md.ss, md.kl[0]-1)
    mm = 2*md.ss
    for kk in range(ll):
        blas.dgemv('t', &mm, &md.ss, &ONED, &md.direct[kk*md.ss**2*2], &mm, 
                   &dataInTemp1[ordM1Hlf+1+kk*md.ss], &ONE, &ONED, &dataOut[1+kk*md.ss], &ONE)
    for ii in range(ll*md.ss+1, pl.nData-1):
        kk = (ii-1)/md.ss
        nn = (ii-1) - kk*md.ss
        mm = min(pl.nData-kk*md.ss-1, 2*md.ss)
        for jj in range(nn+1, mm):
            dataOut[ii] += md.direct[kk*md.ss**2*2+md.ss*2*nn+jj]*dataInTemp1[ordM1Hlf+ii+jj-nn]
    mm = min(pl.nData, 2*md.ss+1)
    for jj in range(1, mm):
        dataOut[0] += md.direct0[jj-1]*dataInTemp1[ordM1Hlf+jj]

    # End correction left singular end
    # TODO maybe BLAS
    for ii in range(pl.nData-1):
        for jj in range(md.order):
            dataOut[ii] += md.coeffsSing[md.order*ii+jj]*dataInTemp1[ii+jj]

    # End correction right nonsingular end
    mm = pl.nData-1
    blas.dgemv('t', &md.order, &mm, &ONED, md.coeffsNonsing, &md.order, 
               &dataInTemp1[pl.nData-1+ordM1Hlf-ordM1HlfIn], &ONE, &ONED, dataOut, &ONE)

    free(dataInTemp1)

    return 0


# Destroy FMM
cdef int destroy_fat_fmmTrapEndCorr(abel_plan* pl) nogil except -1:

    cdef:
        methodData_FMM* md = <methodData_FMM*> pl.methodData

    free(md.chebRoots)
    free(md.kl)
    free(md.klCum)
    free(md.mtmp)
    free(md.mtmm) 
    free(md.ltp)    
    free(md.mtlk)
    free(md.direct)
    free(md.direct0)  
    free(md.coeffsSing)
    free(md.coeffsNonsing)
    free(md.coeffsFilter)
    free(md)

    return 0



#############################################################################################################################################

# Cubic interpolation
cdef inline double interpCubic(double x, int incx, double* p) nogil:

    return p[1*incx] + 0.5 * x*(p[2*incx] - p[0*incx] + 
                                x*(2.0*p[0*incx] - 5.0*p[1*incx] + 4.0*p[2*incx] - p[3*incx] + 
                                   x*(3.0*(p[1*incx] - p[2*incx]) + p[3*incx] - p[0*incx])))


# Polynomial inter-/extrapolation on equidistant grid
cdef double polint(double* data, int nData, double xx) nogil:

    cdef:
        int ii, mm, ns
        double* cc = NULL
        double* dd = NULL
        double den, res

    cc = <double*> malloc(nData*sizeof(double))
    dd = <double*> malloc(nData*sizeof(double))
    ns = <int> (xx+0.5)
    ns = min(max(ns,0),nData-1)
    for ii in range(nData):
        cc[ii] = data[ii]
        dd[ii] = data[ii]

    res = data[ns]
    ns -= 1
    for mm in range(1,nData):
        for ii in range(nData-mm):
            den = (dd[ii]-cc[ii+1])/mm
            dd[ii] = (ii+mm-xx)*den
            cc[ii]= (ii-xx)*den
        if 2*(ns+1) < nData-mm:
            res += cc[ns+1]
        else:
            res += dd[ns]
            ns -= 1
    free(cc)
    free(dd)

    return res


# Apply filter; possibly just numerical derivative
cdef int convolve(double* dataIn, int nData, double* dataOut, int order, double* coeffs) nogil:

    cdef:
        int ii, jj

    memset(dataOut, 0, nData*sizeof(double))
    # TODO Maybe DGEMM or FFT here
    for ii in range(nData):
        for jj in range(order):
            dataOut[ii] += coeffs[jj]*dataIn[ii+jj]

    return 0


cdef int _chebRoots(int order, double* roots) nogil:
    
    cdef:
        int ii
    
    for ii in range(order):
        roots[ii] = mf.cos(0.5*co.pi*(2.*ii+1.)/order)
    
    return 0


cdef double _lagrangePInt(double xx, unsigned int ind, double* nodes, unsigned int order) nogil:

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




