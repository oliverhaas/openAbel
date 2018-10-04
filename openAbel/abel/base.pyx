

from libc.stdlib cimport malloc, free

from openAbel.abel.hansenLaw cimport plan_fat_hansenLawOrgLin, execute_fat_hansenLawLinear, destroy_fat_hansenLawLinear
from openAbel.abel.trap cimport plan_fat_trapezoidalDesingConst, execute_fat_trapezoidalDesingConst, \
                                destroy_fat_trapezoidalDesingConst, \
                                plan_fat_trapezoidalEndCorr, execute_fat_trapezoidalEndCorr, \
                                destroy_fat_trapezoidalEndCorr
from openAbel.abel.fmm cimport plan_fat_fmmTrapEndCorr, execute_fat_fmmTrapEndCorr, destroy_fat_fmmTrapEndCorr



ctypedef struct abel_plan:
    int nData
    int forwardBackward
    double shift
    double stepSize
    int method
    double* grid
    void* methodData



############################################################################################################################################
### Fast Abel transforms                                                                                                                 ###
############################################################################################################################################


# Create plan for Abel transform
cdef abel_plan* plan_fat(int nData, int forwardBackward, double shift, double stepSize, 
                         int method = 3, int order = 2, double eps = 1.e-15) nogil except NULL:

    cdef:
        abel_plan* pl
        int ii

    pl = <abel_plan*> malloc(sizeof(abel_plan))
    if NULL == pl:
        with gil:
            raise MemoryError('Malloc ruturned a NULL pointer, probably not enough memory available.')
    else:
        pl.grid = NULL
        pl.methodData = NULL
    pl.nData = nData
    pl.forwardBackward = forwardBackward
    pl.shift = shift
    pl.stepSize = stepSize
    pl.method = method
    pl.grid = <double*> malloc(nData*sizeof(double))
    if NULL == pl.grid:
        free(pl)
        with gil:
            raise MemoryError('Malloc ruturned a NULL pointer, probably not enough memory available.')
    for ii in range(nData):
        pl.grid[ii] = (ii+shift)*stepSize

    with gil:
        try:
            if pl.method == 0:
                plan_fat_trapezoidalDesingConst(pl)
            elif pl.method == 1:
                plan_fat_hansenLawOrgLin(pl)
            elif pl.method == 2:
                plan_fat_trapezoidalEndCorr(pl, order = order)
            elif pl.method == 3:
                plan_fat_fmmTrapEndCorr(pl, order = order, eps = eps)
            else:
                with gil:
                    raise NotImplementedError('Method not implemented for given parameters.')
        except:
            free(pl.grid)
            free(pl)
            raise

    return pl


# Execute given plan for Abel transform
cdef int execute_fat(abel_plan* pl, double* dataIn, double* dataOut, int leftBoundary = 0, int rightBoundary = 0) nogil except -1:

    if NULL == pl:
        with gil:
            raise TypeError('Input plan is NULL.')

    if pl.method == 0:
        execute_fat_trapezoidalDesingConst(pl, dataIn, dataOut, leftBoundary, rightBoundary)
    elif pl.method == 1:
        execute_fat_hansenLawLinear(pl, dataIn, dataOut)
    elif pl.method == 2:
        execute_fat_trapezoidalEndCorr(pl, dataIn, dataOut, leftBoundary, rightBoundary)
    elif pl.method == 3:
        execute_fat_fmmTrapEndCorr(pl, dataIn, dataOut, leftBoundary, rightBoundary)
    else:
        with gil:
            raise NotImplementedError('Method not implemented for given parameters.')


# Destroy given plan for Abel transform
cdef int destroy_fat(abel_plan* pl) nogil except -1:

    if NULL == pl:
        with gil:
            raise TypeError('Input plan is NULL.')

    if pl.method == 0:
        destroy_fat_trapezoidalDesingConst(pl)
    elif pl.method == 1:
        destroy_fat_hansenLawLinear(pl)
    elif pl.method == 2:
        destroy_fat_trapezoidalEndCorr(pl)
    elif pl.method == 3:
        destroy_fat_fmmTrapEndCorr(pl)
    else:
        with gil:
            raise NotImplementedError('Method not implemented for given parameters.')

    free(pl.grid)
    free(pl)    

    return 0



