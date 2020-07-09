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

a = np.zeros((3, 3))
print(a[1,1])

a[0, 0] = 1
a[0, 1] = 0.903436226742437
a[0, 2] = 0.903345883119763

a[1, 0] = 0.903436226742437
a[1, 1] = 1
a[1, 2] = 0.903345883119763

a[2, 0] = 0.903345883119763
a[2, 1] = 0.903345883119763
a[2, 2] = 1

for i in range(0, 3):
    det_value = my.determinant(a, i+1)


