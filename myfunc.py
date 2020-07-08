def my_sum(a, b):
    result = a + b
    return result


# static variable in ran3
ma = [0 for i in range(56)]
inext = 0
inextp = 0
iff = 0


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
    if idum < 0:
        iset = 0




