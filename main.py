import myfunc as my
import numpy as np

# test
a = 2
b = 3
c = my.my_sum(a, b)
print(c)

seed = 1

for i in range(1, 11):
    print(my.gasdev(seed))

narray = 3
correlmat = np.zeros((3, 3))
ecorr = np.zeros((3, 3))

tmpcorrelmat = np.zeros((narray * narray))
d = np.zeros(narray)
v = np.zeros((narray * narray))
tmpmatrix = np.zeros((narray, narray))
tmpmatrix2 = np.zeros((narray, narray))

correlmat[0, 0] = 1
correlmat[0, 1] = 0.903436226742437
correlmat[0, 2] = 0.903345883119763

correlmat[1, 0] = 0.903436226742437
correlmat[1, 1] = 1
correlmat[1, 2] = 0.903345883119763

correlmat[2, 0] = 0.903345883119763
correlmat[2, 1] = 0.903345883119763
correlmat[2, 2] = 1

for i in range(0, 3):
    det_value = my.determinant(correlmat, i+1)
    print(det_value)

my.choleskydecompose(ecorr, correlmat, 3)
print(ecorr)

for i in range(0, narray):
    for j in range(0, narray):
        tmpcorrelmat[i*narray + j] = correlmat[i, j]

it_max = 100
my.jacobi_eigenvalue(narray, tmpcorrelmat, it_max, v, d)

