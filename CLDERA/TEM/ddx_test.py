import numpy as np
import pdb

x       = np.array([1, 2, 4, 7, 8, 10])
xdiff   = np.diff(x)
ii      = np.arange(len(x) - 2) + 1
xsqd    = np.diff(xdiff**2)
f       = x**2
fdiff   = np.diff(f)

y       = np.arange(len(x)) # "index" coordinate
ydiff   = np.diff(y)
ysqd    = np.diff(ydiff**2)


dfdx1 = np.gradient(f, x)                # numpy (uneven finite diff)
dfdx2 = np.gradient(f) / np.gradient(x)  # numpy (assume even finite diff)
dfdx3 = lambda f,x,xdiff: [(f[1]-f[0])/xdiff[0], 
                    [(f[i+1] - f[i-1])/(xdiff[i]+xdiff[i-1]) for i in ii], 
                    (f[-1]-f[-2])/xdiff[-1]]              # even finite diff

dfdx4 = lambda f,x,xdiff,xsqd: [(f[1]-f[0])/xdiff[0], 
                     [(xdiff[i-1]**2*f[i+1] + (xsqd[i-1]*f[i]) - xdiff[i]**2*f[i-1])/
                     (xdiff[i-1]*xdiff[i]*(xdiff[i-1]+xdiff[i])) for i in ii], 
                     (f[-1]-f[-2])/xdiff[-1]]              # uneven finite diff

pdb.set_trace()
   
