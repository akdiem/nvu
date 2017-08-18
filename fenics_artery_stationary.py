# -*- coding: utf-8 -*-

from fenics import *
from dolfin import *
from mshr import *

from matplotlib import pyplot as plt
import numpy as np
from scipy.interpolate import interp1d
from subprocess import run

parameters["plotting_backend"] = "matplotlib"
prm = parameters['krylov_solver']
prm['absolute_tolerance'] = 1e-7
prm['relative_tolerance'] = 1e-4
prm['maximum_iterations'] = 1000

#list_linear_solver_methods()

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
ast = 5*um/r0

Y = 1.0e6 * pa
nu = 0.49
lam = Y*nu/((1+nu)*(1-2*nu))
mu = Y/(2*(1+nu))
Mu = mu/mu
Lam = lam/mu

# Create mesh and define function space
N = 512
deg = 2
elem = "Lagrange"
geom = Rectangle(Point(0, 0), Point(W, Le))
mesh = generate_mesh(geom, N)
print(mesh)

V = VectorFunctionSpace(mesh, elem, deg)

# Define boundaries
def astr_boundary(x, on_boundary):
    return near(x[0], W) and x[1] < Le/2+ast/2 and x[1] > Le/2-ast/2


def fixed_boundary(x, on_boundary):
    return near(x[1], 0) or near(x[1], Le)

disp_bc = Expression(('d', '0'), d=Disp, degree=deg)
bc1 = DirichletBC(V, disp_bc, astr_boundary, method='pointwise')
bc2 = DirichletBC(V, Constant((0, 0)), fixed_boundary)
bcs = [bc1, bc2]

# Stress and strain
def epsilon(u):
    return 0.5*(nabla_grad(u) + nabla_grad(u).T)

def sigma(u):
    return Lam*nabla_div(u)*Identity(d) + 2*Mu*epsilon(u)

# Define variational problem
u = TrialFunction(V)
d = u.geometric_dimension()
v = TestFunction(V)
f = Constant((0, 0))
a = inner(sigma(u), epsilon(v))*dx
L = dot(f, v)*dx

# Create VTK file for saving solution
ufile = File('data/mesh/u_%d.pvd' % N)
mfile = File('data/mesh/mises_%d.pvd' % N)

# Compute solution
u = Function(V)
solve(a == L, u, bcs)
    
# Calculate stress
W = TensorFunctionSpace(mesh, elem, deg)
sig = sigma(u) - (1./3)*tr(sigma(u))*Identity(d)  # deviatoric stress
#[s11, s12, s21, s22] = sig.split(True)
    
# von Mises stress
von_Mises = sqrt(3./2*inner(sig, sig))
V = FunctionSpace(mesh, elem, deg)
von_Mises = project(von_Mises, V)
    
# Save to file and plot solution
ufile << u
mfile << von_Mises