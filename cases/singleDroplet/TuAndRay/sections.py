import numpy as np
import sys
from scipy import special

N = int(sys.argv[1])
yMin = float(sys.argv[2])
yMax = float(sys.argv[3])

dcm = float(sys.argv[4])
Z = float(sys.argv[5])
rhol = float(sys.argv[6])

y = np.power(10, np.linspace(np.log10(yMin), np.log10(yMax), N+1))
x = np.exp((np.log(y[:-1])+np.log(y[1:]))/2.0)
d = np.power(x/rhol/np.pi*6.0, 1.0/3.0)

x0 = dcm**3*np.pi/6.0*rhol

ii = -1

for i,xi in enumerate(x):

    if xi <= x0 and x[i+1] > x0:

        ii = i
        break

if ii == -1:

    raise Except('Range does not cover point')

x1 = x[ii]
x2 = x[ii+1]

d1 = d[ii]
d2 = d[ii+1]

M2 = Z/(x2+Z+x1*(dcm-d2)/(d1-dcm))
M1 = M2*(dcm-d2)/(d1-dcm)

length = str(int(np.log10(N)+1))

for i,xi in enumerate(x):

    if i == ii:
        Mi = M1

    elif i == ii+1:
        Mi = M2

    else:
        Mi = 0.0

    print(("%0"+length+"i")%(i+1), '%16.12e'%xi, '%16.12e'%Mi)
