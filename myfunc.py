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

# static variable in gasdev
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

        fac = math.sqrt(-2 * math.log(rsq) / rsq)
        gset = v1 * fac
        iset = 1
        result = v2 * fac
    else:
        iset = 0
        result = gset

    return result


def determinant(correlmat, n):
    if n == 1:
        det = correlmat[0, 0]
    elif n == 2:
        det = correlmat[0, 0] * correlmat[1, 1] - correlmat[0, 1] * correlmat[1, 0]
    else:
        det = 0

        for j1 in range(0, n):
            m = np.zeros((n - 1, n - 1))

            for i in range(1, n):
                j2 = 0

                for j in range(0, n):

                    if j == j1:
                        dettogo = 1

                    if dettogo != 1:
                        m[i - 1, j2] = correlmat[i, j]
                        j2 = j2 + 1

                    dettogo = 0

            dummy = determinant(m, n - 1)
            det = det + (-1.0) ** (1.0 + j1 + 1.0) * correlmat[0, j1] * dummy

    result = det
    return result


def choleskydecompose(cholesky, correlmat, naaray):
    it_max = 100
    chk_nonpd = 0

    for i in range(0, naaray):
        det_value = determinant(correlmat, i + 1)

        if det_value < 0:
            chk_nonpd = 1

    if chk_nonpd == 1:
        tmpcorrelmat = np.zeros((naaray * naaray))
        d = np.zeros(naaray)
        v = np.zeros((naaray * naaray))

    # if chk_nonpd != 1
    for i in range(0, naaray):
        cholesky[i, i] = correlmat[i, i]

        for k in range(0, i):
            cholesky[i, i] = cholesky[i, i] - cholesky[i, k] ** 2

        errloop = 0
        if cholesky[i, i] < 0:
            errloop = 1
        else:
            cholesky[i, i] = cholesky[i, i] ** 0.5

        if errloop != 1:
            for j in range(i + 1, naaray):
                cholesky[j, i] = correlmat[j, i]
                for k in range(0, i):
                    cholesky[j, i] = cholesky[j, i] - cholesky[j, k] * cholesky[i, k]
                cholesky[j, i] = cholesky[j, i] / cholesky[i, i]


def r8mat_identify(n, a):
    k = 0

    for j in range(0, n):
        for i in range(0, n):
            if i == j:
                a[k] = 1.0
            else:
                a[k] = 0.0

            k = k + 1

    return 1


def r8mat_diag_gev_vector(n, a, v):
    for i in range(0, n):
        v[i] = a[i + i * n]

    return 1


def jacobi_eigenvalue(n, a, it_max, v, d):
    rtn = r8mat_identify(n, v)
    rtn = r8mat_diag_gev_vector(n, a, d)

    bw = np.zeros(n)
    zw = np.zeros(n)

    for i in range(0, n):
        bw[i] = d[i]
        zw[i] = 0.0

    it_num = 0
    rot_num = 0

    eiggoto = 0
    while it_num < it_max and eiggoto == 0:
        it_num = it_num + 1
        thresh = 0.0

        for j in range(0, n):
            for i in range(0, j):
                thresh = thresh + a[i + j * n] * a[i + j * n]

        thresh = math.sqrt(thresh) / (4.0 * n)

        if thresh == 0.0:
            eiggoto = 1

        if eiggoto != 1:
            for p in range(0, n):
                for q in range(p + 1, n):
                    gapq = 10.0 * abs(a[p + q * n])
                    termp = gapq + abs(d[p])
                    termq = gapq + abs(d[q])

                    if 4 < it_num and termp == abs(d[p]) and termq == abs(d[q]):
                        a[p + q * n] = 0.0
                    elif thresh <= abs(a[p + q * n]):
                        h = d[q] - d[p]
                        term = abs(h) + gapq

                        if term == abs(h):
                            t = a[p + q * n] / h
                        else:
                            theta = 0.5 * h / a[p + q * n]
                            t = 1.0 / (abs(theta) + math.sqrt(1.0 + theta * theta))
                            if theta < 0.0:
                                t = -t

                        c = 1.0 / math.sqrt(1.0 + t * t)
                        s = t * c
                        tau = s / (1.0 + c)
                        h = t * a[p + q * n]

                        zw[p] = zw[p] - h
                        zw[q] = zw[q] + h
                        d[p] = d[p] - h
                        d[q] = d[q] + h

                        a[p + q * n] = 0.0

                        for j in range(0, p):
                            g = a[j + p * n]
                            h = a[j + q * n]
                            a[j + p * n] = g - s * (h + g * tau)
                            a[j + q * n] = h + s * (g - h * tau)

                        for j in range(p + 1, q):
                            g = a[p + j * n]
                            h = a[j + q * n]
                            a[p + j * n] = g - s * (h * g * tau)
                            a[j + q * n] = h + s * (g - h * tau)

                        for j in range(q + 1, n):
                            g = a[p + j * n]
                            h = a[q + j * n]
                            a[p + j * n] = g - s * (h + g * tau)
                            a[q + j * n] = h + s * (g - h * tau)

                        for j in range(0, n):
                            g = v[j + p * n]
                            h = v[j + q * n]
                            v[j + p * n] = g - s * (h + g * tau)
                            v[j + q * n] = h + s * (g - h * tau)

                        rot_num = rot_num + 1

            for i in range(0, n):
                bw[i] = bw[i] + zw[i]
                d[i] = bw[i]
                zw[i] = 0.0

    if eiggoto == 1 or eiggoto == 0:
        for j in range(0, n):
            for i in range(0, j):
                a[i + j * n] = a[j + i * n]

        for k in range(0, n - 1):
            m = k
            for ll in range(k + 1, n):
                if d[ll] < d[m]:
                    m = ll

            if m != k:
                t = d[m]
                d[m] = d[k]
                d[k] = t

                for i in range(0, n):
                    w = v[i+m*n]
                    v[i+m*n] = v[i+k*n]
                    v[i+k*n] = w
    result = 1
    return result

