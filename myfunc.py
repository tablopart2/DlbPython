import math as math
import numpy as np


def my_sum(a, b):
    result = a + b
    return result


# static variable in ran3
ma = [0 for i in range(56)]
inext = 0
inextp = 0
iff = 0

iset = 0
gset = 0


def ran3(idum):
    mbig = 1000000000
    mseed = 161803398
    mz = 0
    fac = 1 / mbig

    global ma
    global inext
    global inextp
    global iff

    if idum < 0 or iff == 0:
        iff = 1
        mj = mseed - abs(idum)
        mj = mj % mbig
        ma[55] = mj
        mk = 1

        for i in range(1, 55):
            ii = (21 * i) % 55
            ma[ii] = mk
            mk = mj - mk

            if mk < mz:
                mk = mk + mbig

            mj = ma[ii]

        for k in range(1, 5):
            for i in range(1, 56):
                ma[i] = ma[i] - ma[1 + (i + 30) % 55]

                if ma[i] < mz:
                    ma[i] = ma[i] + mbig

        inext = 0
        inextp = 31
        idum = 1

    inext = inext + 1

    if inext == 56:
        inext = 1

    inextp = inextp + 1

    if inextp == 56:
        inextp = 1

    mj = ma[inext] - ma[inextp]

    if mj < mz:
        mj = mj + mbig

    ma[inext] = mj

    return mj * fac


def gasdev(idum):
    global iset
    global gset

    rsq = 0

    if idum < 0:
        iset = 0

    if iset == 0:
        while rsq >= 1 or rsq == 0:
            v1 = 2 * ran3(idum) - 1
            v2 = 2 * ran3(idum) - 1
            rsq = v1 * v1 + v2 * v2

        fac = math.sqrt(-2*math.log(rsq)/rsq)
        gset = v1*fac
        iset = 1
        result = v2*fac
    else:
        iset = 0
        result = gset

    return result


def determinant(correlmat, n):
    if n == 1:
        det = correlmat[0,0]
    elif n == 2:
        det = correlmat[0, 0]*correlmat[1, 1] - correlmat[0, 1]*correlmat[1, 0]
    else:
        det = 0

        for j1 in range(0, n):
            m = np.zeros((n-1, n-1))

            for i in range(1, n):
                j2 = 0

                for j in range(0, n):

                    if j == j1:
                        dettogo = 1

                    if dettogo != 1:
                        m[i-1, j2] = correlmat[i, j]
                        j2 = j2 + 1

                    dettogo = 0

            dummy = determinant(m, n - 1)
            det = det + (-1.0) ** (1.0+j1+1.0)*correlmat[0, j1]*dummy

    result = det
    return result


# def choleskydecompose(cholesky,naaray):
#
#     it_max = 100
#
#     chk_nonPD = 0
#
#     for i in range(0, naaray):







