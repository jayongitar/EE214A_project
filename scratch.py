from scipy.interpolate import RegularGridInterpolator
from numpy import linspace, zeros, array
x = linspace(1,4,11)
y = linspace(1,7,22)
z = linspace(1,9,33)
V = zeros((11,22,33))
for i in range(11):
    for j in range(22):
        for k in range(33):
            V[i,j,k] = 100*x[i] + 10*y[j] + z[k]
fn = RegularGridInterpolator((x,y,z), V)
pts = array([[1,1,1],[3,5,7]])
print(fn(pts))