# -*- coding: utf-8 -*-

from fenics import *
from matplotlib import pyplot

parameters ["plotting_backend"] = "matplotlib"

# Units
cm = 1e-2
um = 1e-4 * cm
dyn = 1
pa = 10 * dyn/cm**2

# Scaled variables
r0 = 20*um
R = r0/r0
Le = 10*R
W = 0.2*R
Disp = 4*um/r0
mu0 = 1.13e5*pa
Mu = mu0/mu0
Lam = 5.54e6*pa/mu0

# Create mesh and define function space
nr = 10
nz = 1000 
mesh = RectangleMesh(Point(0, 0), Point(W, Le), nr, nz)
V = VectorFunctionSpace(mesh, 'P', 1)

# Define boundaries
def astr_boundary(x, on_boundary):
    return near(x[0], W) and near(x[1], Le/2)

def fixed_boundary(x, on_boundary):
    return near(x[1], 0) or near(x[1], Le)

#astr_boundary = 'near(x[0], W) && near(x[1], Le/2)'
#fixed_boundary = 'near(x[1], 0) || near(x[1], Le)'
bc1 = DirichletBC(V, Constant((Disp, 0)), astr_boundary, method='pointwise')
bc2 = DirichletBC(V, Constant((0, 0)), fixed_boundary)
bcs = [bc1, bc2]

# Stress and strain

def epsilon(u):
    return 0.5*(nabla_grad(u) + nabla_grad(u).T)

def sigma(u):
    return 0.5*Lam*nabla_div(u)*Identity(d) + 2*Mu*epsilon(u)

# Define variational problem
u = TrialFunction(V)
d = u.geometric_dimension()
v = TestFunction(V)
f = Constant((0, 0))
T = Constant((0, 0))
a = inner(sigma(u), epsilon(v))*dx
L = dot(f, v)*dx + dot(T, v)*ds

# Compute solution
u = Function(V)
solve(a == L, u, bcs)

# Plot solution
u = u
p1 = plot(u, title='Displacement', mode='displacement')
pyplot.colorbar(p1)
pyplot.show()

# Plot stress
W = TensorFunctionSpace(mesh, "Lagrange", 2)
sig = project(sigma(u), W)
[s11, s12, s21, s22] = sig.split(True)
#s_gd = s11.dx(1)
#p_bm = -s_gd*mu0*um/pa
#p2 = plot(s11, title='Stress (Pa/um)')
#pyplot.colorbar(p2)
#pyplot.show()

#print(Le/50)


#
# Save solution to file in VTK format
File('data/sigr_1.pvd') << s11
#File('data/sigz.pvd') << s22

# Hold plot
interactive()