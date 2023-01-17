

# Solve VWF model
# Author: Edwina Yeo, yeo@maths.ox.ac.uk
# Date: 11/22

import matplotlib.pyplot as plt
from fenics import *
import numpy as np
import os
from scipy.optimize import fsolve
from time import time
import math as math
# update pyplot settings for all plots.
plt.rcParams.update({
    "text.usetex": True,
    "font.size": 16, "font.family": "serif",
    "font.serif": ["Palatino"]})
import math
lc1 = 0.5  # mesh factor elsewhere
lc2 = 0.5  # mesh precision near the stenosis
l1 = 1e-3  # length of flat bit in stenosis
l2 = 1.5 # length of sloped section in stenosis
h = 0.3  # height of stenosis
# mesh grid points
l = 12  # half length of domain
nz =10*l  # in straight region of pipe
nsten = 20  # along the stenosis
nx = 20


set_log_level(30)  # output for fenics solver, only if warnings come out.
degree = 1  # degree of lagrange elements
EPS = 1e-4
geomfoldername = ""  # location to store data/figures in

# VWF parameters which are fitted to lippok et al.
xi   =     0.0626
kappa = 0.001 # diffusivity
hatb =      1.7545e-04
gamma_star=    1.9057e+04  #shear stress at which vwf unfolds
delta=      9.9697e-04
LL =       22.6


def write_rect_mesh(l,nx,ny):

    mesh = RectangleMesh(Point(-10, 0.0), Point(l, 1), nx, ny)

    class Noslip(SubDomain):
        def inside(self, x, on_boundary):
            return on_boundary

    # Sub domain for inflow (right)
    class Inflow(SubDomain):
        def inside(self, x, on_boundary):
            return x[0] <= -10 + DOLFIN_EPS and on_boundary

    # Sub domain for outflow
    class Outflow(SubDomain):
        def inside(self, x, on_boundary):
            return x[0] >= l - DOLFIN_EPS and on_boundary

    # # Sub domain for centre
    class Centre(SubDomain):
        def inside(self, x, on_boundary):
            return x[1] <= 0 + DOLFIN_EPS and on_boundary



    # Create mesh functions over the cell facets
    sub_domains = MeshFunction("size_t", mesh, mesh.topology().dim() - 1)

    # Mark all facets as sub domain 3
    sub_domains.set_all(3)

    # Mark no-slip facets as sub domain 0, 0.0
    noslip = Noslip()
    noslip.mark(sub_domains, 0)

    # Mark inflow as sub domain 1, 01
    inflow = Inflow()
    inflow.mark(sub_domains, 1)

    # Mark outflow as sub domain 2, 0.2, True
    outflow = Outflow()
    outflow.mark(sub_domains, 2)

    # Mark no-slip facets as sub domain 0, 0.0
    centre = Centre()
    centre.mark(sub_domains, 4)

    return sub_domains, mesh

def geo_write(l1, l2, h, l, geomfoldername,lc,lc2):

    k = open(geomfoldername + "/mesh.geo", "w")
    

    L=['//+\n SetFactory("OpenCASCADE");\n'+ \
       ' \nl=' + str(l) + ';\nl1=' + str(l1) + ';\nh=' + str(h) + ';\nl2=' + str( l2) +';\nlc=' + str( lc) +';\nlc2=' + str( lc2) +''+ \
       ';\n//+\n//+\nPoint(5) = {0, 1-h, 0, lc};'+ \
       '\n//+\nPoint(6) = {0, 1-2*h, 0, lc};'+ \
       '\n//+\nPoint(7) = {l2, 1, 0, lc};'+ \
       '\n//+\nPoint(8) = {l2+1, 1, 0, lc};'+ \
       '\n//+\nPoint(9) = {-(l2+1), 1, 0, lc};'+ \
       '\n//+\nPoint(10) = {-l2, 1, 0, lc};'+ \
       '\n//+\nPoint(11) = {-l1, 1-h, 0, lc};'+ \
       '\n//+\nPoint(12) = {l1, 1-h, 0, lc};'+ \
       '\n//+\nBSpline(5) = {9, 10, 6, 7, 8};'+ \
       '\n//+\nTranslate {-l1, 0, 0} {\n'+ \
       ' Duplicata { Curve{5}; }\n}'+ \
       '\n//+\nTranslate {l1, 0, 0} {\n    Curve{5};\n}'+ \
       '\n//+\nLine(7) = {11, 12};'+ \

       '\n//+\nPoint(17) = {l, 1, 0, lc2};'+ \
       '\n//+\nPoint(18) = {l, 0, 0, lc2};'+ \
       '\n//+\nPoint(19) = {-10, 0, 0, lc2};'+ \
       '\n//+\nPoint(20) = {-10, 1, 0, lc2};'+ \
       '\n//+\nLine(8) = {20, 19};'+ \
       '\n//+\nLine(9) = {19, 18};'+ \
       '\n//+\nLine(10) = {18, 17};'+ \
       '\n//+\nLine(11) = {17, 16};'+ \
       '\n//+\nLine(12) = {16, 13};'+ \
       '\n//+\nLine(13) = {13, 20};'+ \
       '\n//+\nBooleanFragments{ Curve{6}; Delete; }{ Curve{7}; Delete; }'+ \
       '\n//+\nBooleanFragments{ Curve{5};Delete; }{ Curve{7}; Delete; }'+ \
       '\n//+\nRecursive Delete {\n    Curve{16};\n}'+ \
       '\n//+\nRecursive Delete {\n    Curve{15};\n}'+ \
       '\n//+\nRecursive Delete {\n    Curve{12};\n}'+ \
       '\n//+\nCurve Loop(1) = {13, 8, 9, 10, 11, -17, -7, -14};'+ \
       '\n//+\nCurve Loop(2) = {13, 8, 9, 10, 11, -17, -7, -14};'+ \
       '\n//+\nSurface(1) = {2};\n'+ \
       '//+\nPhysical Point(1) = {20};'+ \
       '//+\nPhysical Point(2) = {19};'+ \
       '//+\nPhysical Point(3) = {18};'+ \
       '//+\nPhysical Point(4) = {17};'+ \
       '//+\nPhysical Curve(5) = {8};'+ \
       '//+\nPhysical Curve(6) = {9};'+ \
       '//+\nPhysical Curve(7) = {10};'+ \
       '//+\nPhysical Curve(8) = {11, 17, 14, 13};'+ \
       '//+\nPhysical Surface(9) = {1};'+\
       '\n//+\n Field[1] = Distance;'+\
       '\n//+\nField[1].NodesList = {11,14};'+\
       '\n//+\nField[1].EdgesList = {7};'+\
       '\n//+\nField[2] = Threshold;'+\
       '\n//+\nField[2].IField = 1;'+\
       '\n//+\nField[2].LcMin = lc / 2;'+\
       '\n//+\nField[2].LcMax = lc;'+\
       '\n//+\nField[2].DistMin = 0.15;'+\
       '\n//+\nField[2].DistMax = 0.5;'+\
       '\n//+\nField[6] = Box;'+\
       '\n//+\nField[6].VIn = lc /3;'+\
       '\n//+\nField[6].VOut = lc;'+\
       '\n//+\nField[6].XMin = -l1-l2/2;'+\
       '\n//+\nField[6].XMax = -l1+0.2;'+\
       '\n//+\nField[6].YMin =0;'+\
       '\n//+\nField[6].YMax =h+0.1;'+\
       '\n//+\nField[6].Thickness = 0.3;'+\
       '\n//+\nField[8] = Box;'+\
       '\n//+\nField[8].VIn = lc /3;'+\
       '\n//+\nField[8].VOut = lc;'+\
       '\n//+\nField[8].XMin = l1+l2;'+\
       '\n//+\nField[8].XMax = l1-0.2;'+\
       '\n//+\nField[8].YMin =0;'+\
       '\n//+\nField[8].YMax =h+0.1;'+\
       '\n//+\nField[8].Thickness = 0.2;'+\

       '\n//+\nField[7] = Min;'+\
       '\n//+\nField[7].FieldsList = {2,6,8};'+\
       '\n//+\nBackground Field = 7;']
    k.writelines(L)
    k.close()# Takes .geo file and writes .msh file, then constructs fenics mesh .xml with labeled subdomains

# uses mesh.geo to write mesh.msh for use by fenics
def write_mesh(l1, l2, h, geomfoldername, l,lc):
    # takes in geometry options l1, l2, h and mesh.geo file in geomfilename and writes subspace using lengh l.
    
    # write mesh
    geo_write(l1, l2, h, l, geomfoldername,lc,lc*2)

    # create mesh and import to Fenics format
    # note mesh must be in gmsh 2 format. We use GMSH 3.0.6.
    # command = "/home/yeo/software/gmsh-3.0.6-Linux64/bin/./gmsh -2 " + geomfoldername + "mesh.geo"



    command = "gmsh -2 " + geomfoldername + "/mesh.geo -format msh2 -v 0"
    os.system(command)
    # os.system(' dolfin.cpp.common.set_log_level(30)')
    os.system("dolfin-convert " + geomfoldername + "/mesh.msh " + geomfoldername + "/mesh.xml")
    mesh = Mesh(geomfoldername + "/mesh.xml")
#    for i in range(0,refinement_level):
#        mesh=refine(mesh)
    # # Sub domain for no-slip (mark whole boundary, inflow and outflow will overwrite)
    class Noslip(SubDomain):
        def inside(self, x, on_boundary):
            return on_boundary

    # Sub domain for inflow (right)
    class Inflow(SubDomain):
        def inside(self, x, on_boundary):
            return x[0] <= -10 + DOLFIN_EPS and on_boundary

    # Sub domain for outflow
    class Outflow(SubDomain):
        def inside(self, x, on_boundary):
            return x[0] >= l - DOLFIN_EPS and on_boundary

    # # Sub domain for centre
    class Centre(SubDomain):
        def inside(self, x, on_boundary):
            return x[1] <= 0 + DOLFIN_EPS and on_boundary

    # Create mesh functions over the cell facets
    sub_domains = MeshFunction("size_t", mesh, mesh.topology().dim() - 1)

    # Mark all facets as sub domain 3
    sub_domains.set_all(3)

    # Mark no-slip facets as sub domain 0, 0.0
    noslip = Noslip()
    noslip.mark(sub_domains, 0)

    # Mark inflow as sub domain 1, 01
    inflow = Inflow()
    inflow.mark(sub_domains, 1)

    # Mark outflow as sub domain 2, 0.2, True
    outflow = Outflow()
    outflow.mark(sub_domains, 2)

    # Mark no-slip facets as sub domain 0, 0.0
    centre = Centre()
    centre.mark(sub_domains, 4)

    return sub_domains, mesh
# uses mesh.geo to write mesh.msh for use by fenics

def get_base_dofs(mesh, subdomains, boundary_index):
    # access base coordinates, and dof of mesh to access the solution at these points.
    # Note arrays need to be reordered for plotting.
    V = FunctionSpace(mesh, 'CG', 1)

    x_base = interpolate(Expression('x[0]', degree=1), V)
    y_base = interpolate(Expression('x[1]', degree=1), V)

    # Unique vertices of the line (by indices)
    vs = list(set(sum((f.entities(0).tolist() for f in SubsetIterator(subdomains, boundary_index)), [])))
    # Get degrees of freedom associated to vertices (by indices)
    v2d = vertex_to_dof_map(V)
    d = v2d[vs]  # list of dofs to access solution at later
    # Dof values
    x_base = x_base.vector().get_local()[d]
    y_base = y_base.vector().get_local()[d]

    return d, x_base, y_base

# access all the mesh coordinates
def get_mesh_dofs(mesh, V):
    # access base coordinates, and dof of mesh to access the solution at these points.
    # Note arrays need to be reordered for plotting.
    x_base = interpolate(Expression('x[0]', degree=1), V)
    y_base = interpolate(Expression('x[1]', degree=1), V)

    # Get degrees of freedom associated to vertices (by indices)
    v2d = vertex_to_dof_map(V)
    # Dof values
    x_mesh = x_base.vector().get_local()[v2d]
    y_mesh = y_base.vector().get_local()[v2d]

    return x_mesh, y_mesh

# extract value of variable on base of mesh
def get_base_val(var, d):
    # d is a list of dofs of the base coordinates
    val=var.vector().get_local()[d]
    val=np.array(val,dtype=np.float64)
    return val

def A_parabolic(L, Re, r):
    # algebraic equation for Azz is aAzz^3+bAzz-c=0
    # L,De,r=(i for i in data)
    sr=2*r#shear rate at the inlet
    tau=Re*xi*((np.tanh(hatb*Re*(sr-gamma_star/Re))+1)/2+delta)

    # print(tau)
    # print(tau)
    wr = -2 * r

    def func(x):
        f=L*L/(L*L-x[2]-x[0])
        a=L*L/(L*L-1)
        B=[f*x[0]-a,-x[0]*wr+1/tau*(f*x[1]),-2*x[1]*wr+1/tau*(f*x[2]-a)]

        return B

    A = fsolve(func, [1, 0,1])
    return A
# inlet values for the FENE-p equations
def inlet_A(mesh, Re, L):
    # print(Re)

    V = FunctionSpace(mesh, "CG", 1)
    class Azz_inlet(UserExpression):
        def eval(self, values, x):
            A = A_parabolic(L, Re, x[1] )
            values[0] = A[2]

    class Arz_inlet(UserExpression):
        def eval(self, values, x):
            A = A_parabolic(L, Re, x[1] )
            values[0] = A[1]

    class Arr_inlet(UserExpression):
        def eval(self, values, x):
            A = A_parabolic(L, Re, x[1])
            values[0] = A[0]

    Arz_in = Arz_inlet(element=V.ufl_element())
    Azz_in = Azz_inlet(element=V.ufl_element())
    Arr_in = Arr_inlet(element=V.ufl_element())

    return Arr_in, Arz_in, Azz_in

# solve steady FENE-p equation
def fene_init(mesh, u, w, ur, uz, wr, wz, sub_domains, Re_val,sr):
    # Takes input mesh.xml file (string name given)
    # foldername is a string where the data should be saved
    # input velocity is a string name of text file with velocity components evaluated at the mesh vertices. The columns
    # are u,w,ur,uz,wr,wz
    # longthin is bool, true means the system is long and thin and use the

    #  Create velocity function.
    V_u = VectorFunctionSpace(mesh, 'CG', degree + 1)  # space for the velocity field.
    V = FunctionSpace(mesh, 'CG', 1)
    vel = interpolate(Expression(('w', 'u'), degree=3, w=w, u=u), V_u)

    # define mixed function space for the three tensor components
    CG = FiniteElement('CG', triangle, degree)  # discontinuous galerkin for phi
    element = MixedElement([CG, CG, CG])
    VV = FunctionSpace(mesh, element)

    # define functions and test functions (since this is nonlinear)
    w = Function(VV)
    J = TestFunction(VV)

    # split to get components
    Arr, Arz, Azz = split(w)  # trial functions
    v2, v3, v4 = split(J)  # test functions

    #  Solution function.

    r = Expression('x[1]', degree=2)

    Re = Expression('val2',degree=2, val2=Re_val)

    # weak forms for each tensor component

    tau_gamma=Expression('hata * ((tanh(hatb * Re_ * (sr -gamma_star/ Re_)) + 1) / 2 + delta)',degree=2,delta=delta,gamma_star=gamma_star,sr=sr,Re_=Re_val,hatb=hatb,hata=xi)
    TAU= interpolate(tau_gamma,V)
#    file_u = File(path + "/Tau.pvd")
#    file_u << TAU
    # plot_save(TAU,'relax time sr',path)
    # base_plot(x_base,get_base_val(TAU,d),'base relax time',path)
    G = LL *LL/ (LL *LL - Arr - Azz)
    a=LL*LL/(LL*LL-2)
    np.savetxt(path + '/base_tau.txt', get_base_val(TAU, d), fmt='%1.9f')


    Frr = Re * (tau_gamma) * (dot(vel, grad(Arr)) * v2 + kappa * dot(grad(Arr), grad(v2))
                              - 2 * (Arr * ur + Arz * uz) * v2) * r * dx
    Frz = Re * (tau_gamma) * (dot(vel, grad(Arz)) * v3 + kappa * dot(grad(Arz), grad(v3)) \
                              - (Azz * uz +  Arr * wr) * v3) * r * dx
    Fzz = Re * (tau_gamma) * (dot(vel, grad(Azz)) * v4 + kappa * dot(grad(Azz), grad(v4)) \
                              - 2 * (Azz * wz + Arz * wr) * v4) * r * dx

    a2 = Frr + Frz + Fzz
    L2 = (  -( G*Arr-a)* v2 * r * dx -G* Arz * v3 * r * dx - (
            G*Azz-a)* v4 * r * dx)
    F2 = a2 - L2

    return F2,tau_gamma,Re,VV,V,w
# solve steady FENE-p equation
def fene_solve( F2,tau_gamma,Re,sub_domains, Re_val,sr,VV,V,w,steps):
    # Takes input mesh.xml file (string name given)
    # foldername is a string where the data should be saved
    # input velocity is a string name of text file with velocity components evaluated at the mesh vertices. The columns
    # are u,w,ur,uz,wr,wz
    # longthin is bool, true means the system is long and thin and use the
    V_1 = FunctionSpace(mesh, 'CG', 1)

    #  Create velocity function.

    Re_fene_iter = np.linspace(1e-5, 1,steps)  # deborah scales for continuation De=xi*Re
    # De = Expression('val*val2', degree=2, val=xi_iter[0], val2=Re_val)
    # continue up to the correct deborah number
    tau_gamma.sr = sr
    if steps==0:
        tau_gamma.Re_ =  Re_val
        Re.val =  Re_val

        TAU = interpolate(tau_gamma, V)

        Arr_in, Arz_in, Azz_in = inlet_A(mesh,  Re_val, LL)
        bc_rr = DirichletBC(VV.sub(0), Arr_in, sub_domains, 1, "geometric")
        bc_rz = DirichletBC(VV.sub(1), Arz_in, sub_domains, 1, "geometric")
        bc_zz = DirichletBC(VV.sub(2), Azz_in, sub_domains, 1, "geometric")
        bc_old = [bc_rr, bc_rz, bc_zz]
        solve(F2 == 0, w, bc_old,solver_parameters = {"newton_solver":{ "linear_solver" : "mumps"}})

        Arr_h, Arz_h, Azz_h = w.split()
        ARR = interpolate(Arr_h, V_1)
        ARZ = interpolate(Arz_h, V_1)
        AZZ = interpolate(Azz_h, V_1)


        r = Expression('x[1]', degree=2)

        s = TestFunction(V)
        du = TrialFunction(V)
        A = inner(du, s) * r * dx
        S = (Arr_h + Azz_h - 2) * s * r * dx  # +Arz*vel[1]*s*dx
        Ext = Function(V)
        solve(A == S, Ext)
        V_u = VectorFunctionSpace(mesh, 'CG', degree + 1)  # space for the velocity field.


    for k in range(len(Re_fene_iter)):
        tau_gamma.Re_=Re_fene_iter[k]*Re_val
        Re.val = Re_fene_iter[k]*Re_val
#        print('solving fene for ',str( Re_fene_iter[k]*Re_val))
        TAU = interpolate(tau_gamma, V)


        Arr_in, Arz_in, Azz_in = inlet_A(mesh, Re_fene_iter[k]*Re_val , LL)
        bc_rr = DirichletBC(VV.sub(0), Arr_in, sub_domains, 1, "geometric")
        bc_rz = DirichletBC(VV.sub(1), Arz_in, sub_domains, 1, "geometric")
        bc_zz = DirichletBC(VV.sub(2), Azz_in, sub_domains, 1, "geometric")
        bc_old = [bc_rr, bc_rz, bc_zz]
        solve(F2 == 0, w, bc_old,solver_parameters = {"newton_solver":{ "linear_solver" : "mumps"}})

        Arr_h, Arz_h, Azz_h = w.split()


        ARR=interpolate(Arr_h,V_1)
        ARZ=interpolate(Arz_h,V_1)
        AZZ=interpolate(Azz_h,V_1)


        r = Expression('x[1]', degree=2)

        s = TestFunction(V)
        du = TrialFunction(V)
        A = inner(du, s) * r * dx
        S =(Arr_h+Azz_h-2)/2 * s * r * dx  # +Arz*vel[1]*s*dx
        Ext = Function(V)
        solve(A == S, Ext)


        V_u = VectorFunctionSpace(mesh, 'CG', degree + 1)  # space for the velocity field.
    Abar=interpolate(Expression(('(AZZ-1)/sqrt((ARR-1)*(ARR-1)+(AZZ-1)*(AZZ-1))','(ARR-1)/sqrt((ARR-1)*(ARR-1)+(AZZ-1)*(AZZ-1))'),degree=1,ARR=ARR,AZZ=AZZ),V_u)
#    file_ext = File(path+"/Azz.pvd")
#    file_ext << AZZ

#    file_ext = File(path+"/Arr.pvd")
#    file_ext << ARR
#
#    file_ext = File(path+"/Arz.pvd")
#    file_ext << ARZ

    write_txt(mesh,Ext,'extsol',path)
    write_txt(mesh,ARR,'rrsol',path)
    write_txt(mesh,AZZ,'zzsol',path)
    write_txt(mesh,ARZ,'rzsol',path)


    return Ext

def create_level_mesh(mesh_level):
    # mesh is coarse mesh, psi is streamfunction solution from NS equations in P1. mesh level determines the overall
    # fineness of the mesh


    command = "/home/yeo/software/gmsh-3.0.6-Linux64/bin/./gmsh -2 " + geomfoldername + "mesh.geo -bgm level"+str(mesh_level)+'bgm.pos'

    # command = "gmsh -2 " + geomfoldername + "mesh.geo"

    os.system(command)
    os.system("dolfin-convert " + geomfoldername + "mesh.msh " + geomfoldername + "mesh.xml")
    mesh = Mesh(geomfoldername + "mesh.xml")

    # # Sub domain for no-slip (mark whole boundary, inflow and outflow will overwrite)
    class Noslip(SubDomain):
        def inside(self, x, on_boundary):
            return on_boundary

    # Sub domain for inflow (right)
    class Inflow(SubDomain):
        def inside(self, x, on_boundary):
            return x[0] <= -l + DOLFIN_EPS and on_boundary

    # Sub domain for outflow
    class Outflow(SubDomain):
        def inside(self, x, on_boundary):
            return x[0] >= l - DOLFIN_EPS and on_boundary

    # # Sub domain for centre
    class Centre(SubDomain):
        def inside(self, x, on_boundary):
            return x[1] <= 0 + DOLFIN_EPS and on_boundary

    # Create mesh functions over the cell facets
    sub_domains = MeshFunction("size_t", mesh, mesh.topology().dim() - 1)

    # Mark all facets as sub domain 3
    sub_domains.set_all(3)

    # Mark no-slip facets as sub domain 0, 0.0
    noslip = Noslip()
    noslip.mark(sub_domains, 0)

    # Mark inflow as sub domain 1, 01
    inflow = Inflow()
    inflow.mark(sub_domains, 1)

    # Mark outflow as sub domain 2, 0.2, True
    outflow = Outflow()
    outflow.mark(sub_domains, 2)

    # Mark no-slip facets as sub domain 0, 0.0
    centre = Centre()
    centre.mark(sub_domains, 4)

    return sub_domains, mesh

def find_stag_points(x_base,val_base):
    a = np.column_stack((x_base, val_base) ) # reorder points
    a = a[a[:, 0].argsort()]
    x_base = a[:, 0]
    val_ordered = a[:, 1]

    zero_crossings = np.where(np.diff(np.sign(val_ordered)))
    #
    # plt.plot(x_base, val_ordered,'o')
    # plt.plot(x_base[zero_crossings[0]],val_ordered[zero_crossings[0]],'o')


    # plt.show()
    stags=val_ordered[zero_crossings[0]]
    xs=x_base[zero_crossings[0]]
    return stags,xs

def NS(mesh, re,sub_domains):
    # a, Re, ww, bcs = NS_init(mesh, sub_domains,re)

    # Compute solution
    # update the renolds number
    V2 = VectorElement("CG", triangle, 2)
    QQ = FiniteElement("CG", triangle, 1)
    W_space = FunctionSpace(mesh, MixedElement([V2, QQ]))

    # Define Function and TestFunction(s)
    ww = Function(W_space)
    (u, p) = split(ww)
    (v, q) = split(TestFunction(W_space))

    # Define viscosity and bcs
    # p0 = Expression("-4*x[0]+4*l", degree=2, l=L) #inlet pressure
    Uin = Expression(('(1-x[1]*x[1])', '0.0'), degree=2)  # inlet veloctiy field

    bcu_noslip_b = DirichletBC(W_space.sub(0), Constant((0, 0)), sub_domains, 0)  # both velocity components are zero on base
    bcout = DirichletBC(W_space.sub(0).sub(1), 0.0, sub_domains, 2)  # unidirectional flow out
    bc_centre_v = DirichletBC(W_space.sub(0).sub(1), Constant(0), sub_domains, 4)  # centre v=0
    bcu_inflow = DirichletBC(W_space.sub(0), Uin, sub_domains, 1)  # inlet velocity profile
    # bcp_inflow = DirichletBC(W.sub(1), Expression('8*L', degree=1, L=L), sub_domains, 1)  # inlet pressure
    bcp_outflow = DirichletBC(W_space.sub(1), Constant(0.0), sub_domains, 2)  # outlet pressure=0
    bcs = [bcu_noslip_b, bcout, bcu_inflow, bcp_outflow, bc_centre_v]

    # Define variational form
    r = Expression('x[1]', degree=2)
    er = Constant((0.0, 1.0))
    epsilon = sym(nabla_grad(u))
    Re = Expression('val', degree=2, val=re)
    a = +Re * inner(grad(u) * u, v) * r * dx + (
            2 * inner(epsilon, nabla_grad(v)) - div(u) * q - div(v) * p) * r * dx - (
                dot(u, er) * q + dot(v, er) * p) * dx  # 2d cartesian
    Re.val = re
    print('Solving NS equations for Re=', re)
    solve(a == 0, ww, bcs,solver_parameters = {"newton_solver":{ "linear_solver" : "mumps"}})
    # J = derivative(a, ww)

    (uu, p) = ww.split()
    w, u = uu.split()  # split into velocity components
    meshsize=mesh.num_vertices()
    max_w=np.max(w.vector().get_local())
    # print(max_w)
    V_1 = FunctionSpace(mesh, "P", 1)  # space for the shear component vectors


    (uu, p) = ww.split()
#    file_u = File(path + "/u.pvd")
#    file_u << uu
#    file_u = File(path + "1/p.pvd")
#    file_u << p
    W, U = uu.split()  # split into velocity components

    # plotting(uu,'u inside error calc')

    K=VectorFunctionSpace(mesh,'CG',1)
    s = TestFunction(K)
    du = TrialFunction(K)
    A = inner(du, s) * r * dx
    S = dot(grad(U), s) * r * dx
    du = Function(K)
    solve(A == S, du)
    (uz, ur) = du.split()

    dw = TrialFunction(K)
    A = inner(dw, s) * r * dx
    S = dot(grad(W), s) * r * dx
    dw = Function(K)
    solve(A == S, dw)

    (wz, wr) = dw.split()


    gamma = sqrt( (2 * ur * ur) + (2 * wz * wz) + (uz + wr) * (uz + wr))

    k = TrialFunction(V_1)
    s = TestFunction(V_1)

    a = k * s * dx
    L = gamma * s * dx
    A = assemble(a, keep_diagonal=True)
    b = assemble(L)
    A.ident_zeros()

    shear_rate = Function(V_1)
    solve(A, shear_rate.vector(), b)
#    file_u = File(path + "/sr.pvd")
#    file_u << shear_rate
    base_vort=get_base_val(shear_rate,d)
    np.savetxt(path+'/sr.txt',base_vort)
    # +base_plot(x_base,base_vort,'sr',path)
    # # solve for wall shear stress (sigma -(sigma.n)n)=WSS
    n = FacetNormal(mesh)
    tangent = as_vector([n[1], -n[0]])
    # Define the stress tensor and vector
    sigma = 2 * sym(nabla_grad(uu))
    T = sigma * n

    # Define the wall shear stress
    Tt = T - inner(T, n) * n

    k = TrialFunction(K)
    s = TestFunction(K)

    a = inner(k, s) * ds(0)
    L = inner(Tt, s) * ds(0)
    A = assemble(a, keep_diagonal=True)
    b = assemble(L)
    A.ident_zeros()

    shear_stress = Function(K)
    solve(A, shear_stress.vector(), b)

    ss_1, ss_2 = shear_stress.split()
    # ss_1,ss_2=1,2

    V_1=FunctionSpace(mesh,'CG',1)

    # Compute vorticity by L2 projection
    k = TrialFunction(V_1)
    s = TestFunction(V_1)
    r = Expression('x[1]', degree=2)
    x = SpatialCoordinate(mesh)

    a = k * s * dx
    L = (x[1] * Dx(uu[0], 1) + uu[0] - x[1] * Dx((uu[1]), 0)) * s * dx
    vort = Function(V_1)
    solve(a == L, vort)
    # File("vorticity.pvd", "compressed") << vort

#    file_ext = File(path+"/vort.pvd")
#    file_ext << vort
#    base_vort=get_base_val(vort,d)
#    np.savetxt(path+'/basevort.txt',base_vort)
    # base_plot(x_base,base_vort,'base vort',path)
    # Compute stream function on wall
    # Laplace(psi) = -vort
    a = inner(grad(k), grad(s)) * dx
    L = vort * s * dx
    psi = Function(V_1)
    wall = DirichletBC(V_1, 0, sub_domains, 0)
    solve(a == L, psi, wall)
    
    
    L = vort * s * dx
    psi = Function(V_1)
    wall = DirichletBC(V_1, 0, sub_domains, 0)
    solve(a == L, psi, wall)

        # Compute stream function in domain 
    # Laplace(psi) = -vort

    # Compute vorticity by L2 projection
    k = TrialFunction(V_1)
    s = TestFunction(V_1)


    a = inner(grad(s), grad(k))*dx
    L = vort* s*dx
    psi_tot=Function(V_1)
    solve(a == L, psi_tot)
    write_txt(mesh,psi_tot,'psi',path)


 
    # Compute vorticity by L2 projection
    k = TrialFunction(V_1)
    s = TestFunction(V_1)


    a = k * s * ds(0)
    L = (-shear_stress[0]*n[1]+shear_stress[1]*n[0])  * s * ds(0)
    shear_stress_sign = Function(V_1)
    solve(a == L, shear_stress_sign)
    # # Compute u.grad(sr)
    # k = TrialFunction(V_1)
    # s = TestFunction(V_1)
    # file_u = File(path + "/ur.pvd")
    # file_u << ur
    # file_u = File(path + "/wr.pvd")
    # file_u << wr
    # file_u = File(path + "/wz.pvd")
    # file_u << wz
    # file_u = File(path + "/uz.pvd")
    # file_u << uz
    # plot(ur)
    # plt.show()
    # plot(uz)
    # plt.show()
    # plot(wr)
    # plt.show()
    # plot(wz)
    # plt.show()
    return max_w, meshsize,ss_1,ss_2,psi,shear_stress_sign,u, w ,ur,uz,wr,wz,shear_stress,shear_rate

def write_txt(mesh,var,name,path1):
    V = FunctionSpace(mesh, "CG", 1)
    u = interpolate(var, V)
    coords = V.tabulate_dof_coordinates()
    vec = u.vector().get_local()
    outfile = open(path1+'/'+name+'.txt', "w")
    for coord, val in zip(coords, vec):
        print(coord[0], coord[1], val, file=outfile)
    return


h=0.5
l2=2
l1=1.5

l=30+l2+l1
#mesh_level=3
# rs=np.linspace(0,1,100)
# Ars=np.zeros((100,3))
# for ii in range(len(rs)):
#     Ars[ii,:]=A_parabolic(LL,400,rs[ii])
#
#
# plt.plot(Ars[:,1],rs)
# plt.show()
#
# plt.plot(rs,Ars[:,2])
#
# plt.show()
#
#
# l2s=[1]
# for ii in l2s:
#     l2=ii
#re=500
#lcs=[0.4,0.3,0.2,0.1,0.07,0.05,0.04,0.03,0.02,0.01]
#lcs=[0.025]
#
##
## Numerical Validation ----------------------------------------
##
##1) mesh convergence
#
## # # #create sequence of finer meshes
#steps=re
#max_ext=np.zeros((len(lcs),1))
#min_ext=np.zeros((len(lcs),1))
#
#max_ws=np.zeros((len(lcs),1))
#mesh_vertices=np.zeros((len(lcs),1))
#k=len(lcs)-1
#for j in lcs:
#     lc=j
#
#     path1='data/Convergence'
#     if not os.path.exists(path1):
#         os.mkdir(path1)
#     path = path1 + '/Re_' + str(re)+'h-'+str(h)+'l1-'+str(l1)+'l2-'+str(l2)+'lc'+str(lc)
#     if not os.path.exists(path):
#         os.mkdir(path)
#    
#     print('Solving system for mesh lc', lc)
#         
#     sub_domains, mesh = write_mesh(l1, l2, h, path1, l,lc)  # write mesh and domains for solving
#    # sub_domains, mesh = write_rect_mesh(l,100,100)  # write mesh and domains for solving
#    
#     ds = ds(subdomain_data=sub_domains)  # assign surface integration measure
#     d, x_base, y_base = get_base_dofs(mesh, sub_domains, 0)  # get dofs of base and their z values
#     
#
#     (max_w, meshsize,ss_1,ss_2,psi,shear_stress_sign,u, w ,ur,uz,wr,wz,shear_stress,sr) =NS(mesh,re,sub_domains)
#
#
#     write_txt(mesh,u,'usol',path)
#     write_txt(mesh,w,'wsol',path)
#    
#     V_1 = FunctionSpace(mesh, 'CG',  1)  # continous galerkin space for the velocity field.
#     wss = interpolate(Expression('sqrt(ss_1* ss_1+ss_2 * ss_2)',degree=2,ss_1=ss_1,ss_2=ss_2),V_1)
#    # base_plot(x_base,get_base_val(wss,d),'wss base', path)
#     eig_E= interpolate(Expression('(uz+wr)*(uz+wr)/4-ur*wz',degree=2,uz=uz,wr=wr,wz=wz,ur=ur),V_1)
#     eig_W= interpolate(Expression('(uz-wr)*(uz-wr)/4',degree=2,uz=uz,wr=wr),V_1)
#     rot_rate= interpolate(Expression('uz-wr',degree=2,uz=uz,wr=wr,ur=ur,wz=wz),V_1)
#     total_rate= interpolate(Expression('abs((uz-wr))-sr',degree=2,uz=uz,wr=wr,ur=ur,wz=wz,sr=sr),V_1)
#    # rot_rate= interpolate(Expression('sqrt(1/2*((uz-wr)^2+(uz-wr)^2))',degree=2,uz=uz,wr=wr),V_1)
#    
#    
#     write_txt(mesh,rot_rate,'rot',path)
#     write_txt(mesh,total_rate,'tot',path)
#     write_txt(mesh,sr,'sr',path)
#    
#     #  Create velocity function.
#     (F2, tau_gamma, Re, VV, V,W) = fene_init(mesh, u, w, ur, uz, wr, wz, sub_domains, re, sr)
#
#     Ext=fene_solve(F2,tau_gamma,Re,sub_domains,re,sr,VV,V,W,steps)
#
#    
#     max_ws[k,:]=max_w
#     max_ext[k]=max(Ext.vector())
#     min_ext[k]=min(Ext.vector())
#
#     print(min_ext)
#
#     mesh_vertices[k]=meshsize
#     foldername="data/Convergence/"
#     np.savetxt(foldername + 'meshverts4.txt', mesh_vertices)
#     np.savetxt(foldername + 'max_w4.txt',max_ws)
#     np.savetxt(foldername + 'max_e4.txt',max_ext)
#
#     np.savetxt(foldername + 'min_e4.txt',min_ext)
#     
#     k=k-1
#foldername="data/Convergence/"
#np.savetxt(foldername + 'meshverts4.txt', mesh_vertices)
#np.savetxt(foldername + 'max_w4.txt',max_ws)
#np.savetxt(foldername + 'max_e4.txt',max_ext)


##2) continuation convergence 
##
##
##
#
##stepss=[re*1.5,re*2]
##stepss=[re*4,re*0.75]
#stepss=[re*0.2]
##steps=re
#max_ext=np.zeros((len(lcs),1))
#min_ext=np.zeros((len(lcs),1))
#lc=0.02
#
##max_ws=np.zeros((len(lcs),1))
##mesh_vertices=np.zeros((len(lcs),1))
#k=len(lcs)-1
#for j in stepss:
#     steps=int(np.floor(j))
#     path1='data/Convergence'
#     if not os.path.exists(path1):
#         os.mkdir(path1)
#     path = path1 + '/Re_' + str(re)+'h-'+str(h)+'l1-'+str(l1)+'l2-'+str(l2)+'steps'+str(steps)
#     if not os.path.exists(path):
#         os.mkdir(path)
#     print(path)
#     print('Solving system for steps', steps)
#         
#     sub_domains, mesh = write_mesh(l1, l2, h, path1, l,lc)  # write mesh and domains for solving
#    # sub_domains, mesh = write_rect_mesh(l,100,100)  # write mesh and domains for solving
#    
#     ds = ds(subdomain_data=sub_domains)  # assign surface integration measure
#     d, x_base, y_base = get_base_dofs(mesh, sub_domains, 0)  # get dofs of base and their z values
#     
#
#     (max_w, meshsize,ss_1,ss_2,psi,shear_stress_sign,u, w ,ur,uz,wr,wz,shear_stress,sr) =NS(mesh,re,sub_domains)
##
##
#     write_txt(mesh,u,'usol',path)
#     write_txt(mesh,w,'wsol',path)
#    
#     V_1 = FunctionSpace(mesh, 'CG',  1)  # continous galerkin space for the velocity field.
#     wss = interpolate(Expression('sqrt(ss_1* ss_1+ss_2 * ss_2)',degree=2,ss_1=ss_1,ss_2=ss_2),V_1)
#    # base_plot(x_base,get_base_val(wss,d),'wss base', path)
#     eig_E= interpolate(Expression('(uz+wr)*(uz+wr)/4-ur*wz',degree=2,uz=uz,wr=wr,wz=wz,ur=ur),V_1)
#     eig_W= interpolate(Expression('(uz-wr)*(uz-wr)/4',degree=2,uz=uz,wr=wr),V_1)
#     rot_rate= interpolate(Expression('uz-wr',degree=2,uz=uz,wr=wr,ur=ur,wz=wz),V_1)
#     total_rate= interpolate(Expression('abs((uz-wr))-sr',degree=2,uz=uz,wr=wr,ur=ur,wz=wz,sr=sr),V_1)
#    # rot_rate= interpolate(Expression('sqrt(1/2*((uz-wr)^2+(uz-wr)^2))',degree=2,uz=uz,wr=wr),V_1)
#    
#    
#     write_txt(mesh,rot_rate,'rot',path)
#     write_txt(mesh,total_rate,'tot',path)
#     write_txt(mesh,sr,'sr',path)
#    
#     #  Create velocity function.
#     (F2, tau_gamma, Re, VV, V,W) = fene_init(mesh, u, w, ur, uz, wr, wz, sub_domains, re, sr)
#
#     Ext=fene_solve(F2,tau_gamma,Re,sub_domains,re,sr,VV,V,W,steps)
#
#    
##     max_ws[k,:]=max_w
#     max_ext[k]=max(Ext.vector())
#     min_ext[k]=min(Ext.vector())
#
#     print(min_ext)
#
##     mesh_vertices[k]=meshsize
#     foldername="data/Convergence/"
#     np.savetxt(path+ 'max_e_step.txt',max_ext)
#     np.savetxt(path + 'min_e_step.txt',min_ext)
#     
#     k=k-1
##foldername="data/Convergence/"
##np.savetxt(foldername + 'meshverts3.txt', mesh_vertices)
##np.savetxt(foldername + 'max_w3.txt',max_ws)
#np.savetxt(foldername + 'max_e3.txt',max_ext)



#
#
#path1='data'
#
#lc=0.025
#Res=[500]
##Res=[500]
#for i in range(len(Res)):
#    re=Res[i]
#    steps=int(np.floor(re*0.2))
#
#    if not os.path.exists(path1):
#        os.mkdir(path1)
#    path = path1 + '/Re_' + str(re)+'h-'+str(h)+'l1-'+str(l1)+'l2-'+str(l2)
#    if not os.path.exists(path):
#        os.mkdir(path)
#    
#    sub_domains, mesh = write_mesh(l1, l2, h, path1, l,lc)  # write mesh and domains for solving
#        # sub_domains, mesh = write_rect_mesh(l,100,100)  # write mesh and domains for solving
#    
#    ds = ds(subdomain_data=sub_domains)  # assign surface integration measure
#    d, x_base, y_base = get_base_dofs(mesh, sub_domains, 0)  # get dofs of base and their z values
#    
#    
#    (max_w, meshsize,ss_1,ss_2,psi,shear_stress_sign,u, w ,ur,uz,wr,wz,shear_stress,sr) =NS(mesh,re,sub_domains)
#    
#    write_txt(mesh,u,'usol',path)
#    write_txt(mesh,w,'wsol',path)
#    
#    V_1 = FunctionSpace(mesh, 'CG',  1)  # continous galerkin space for the velocity field.
#    wss = interpolate(Expression('sqrt(ss_1* ss_1+ss_2 * ss_2)',degree=2,ss_1=ss_1,ss_2=ss_2),V_1)
#    rot_rate= interpolate(Expression('uz-wr',degree=2,uz=uz,wr=wr,ur=ur,wz=wz),V_1)
#    total_rate= interpolate(Expression('abs((uz-wr))-sr',degree=2,uz=uz,wr=wr,ur=ur,wz=wz,sr=sr),V_1)
#    
#    
#    write_txt(mesh,rot_rate,'rot',path)
#    write_txt(mesh,total_rate,'tot',path)
#    write_txt(mesh,sr,'sr',path)
#    
#
#    np.savetxt(path + '/ss_sign.txt', get_base_val(shear_stress_sign, d), fmt='%1.9f')
#    np.savetxt(path + '/wss.txt', get_base_val(wss, d), fmt='%1.9f')
#    np.savetxt(path + '/x_base.txt', x_base, fmt='%1.9f')
#    #
#    # #  Create velocity function.
#    (F2, tau_gamma, Re, VV, V,W) = fene_init(mesh, u, w, ur, uz, wr, wz, sub_domains, re, sr)
#    Ext=fene_solve(F2,tau_gamma,Re,sub_domains,re,sr,VV,V,W,steps)
#    base_ext=get_base_val(Ext,d)
#    np.savetxt(path + '/base_ext.txt', base_ext, fmt='%1.9f')

# path1='data'
# #
# lc=0.025
# #Res=[200,300,400]
# l2s=[2,3,4,5]
# #Res=[500]
# re=500
# for i in range(len(l2s)):
#    l2=l2s[i]
#    steps=int(np.floor(re*0.2))

#    if not os.path.exists(path1):
#        os.mkdir(path1)
#    path = path1 + '/Re_' + str(re)+'h-'+str(h)+'l1-'+str(l1)+'l2-'+str(l2)
#    if not os.path.exists(path):
#        os.mkdir(path)
   
#    sub_domains, mesh = write_mesh(l1, l2, h, path, l,lc)  # write mesh and domains for solving
#        # sub_domains, mesh = write_rect_mesh(l,100,100)  # write mesh and domains for solving
   
#    ds = ds(subdomain_data=sub_domains)  # assign surface integration measure
#    d, x_base, y_base = get_base_dofs(mesh, sub_domains, 0)  # get dofs of base and their z values    
# #    (max_w, meshsize,ss_1,ss_2,psi,shear_stress_sign,u, w ,ur,uz,wr,wz,shear_stress,sr) =NS(mesh,re,sub_domains)
#    np.savetxt(path+'/y_base.txt',y_base)

#    write_txt(mesh,u,'usol',path)
#    write_txt(mesh,w,'wsol',path)
#    
#    V_1 = FunctionSpace(mesh, 'CG',  1)  # continous galerkin space for the velocity field.
#    wss = interpolate(Expression('sqrt(ss_1* ss_1+ss_2 * ss_2)',degree=2,ss_1=ss_1,ss_2=ss_2),V_1)
#    rot_rate= interpolate(Expression('uz-wr',degree=2,uz=uz,wr=wr,ur=ur,wz=wz),V_1)
#    total_rate= interpolate(Expression('abs((uz-wr))-sr',degree=2,uz=uz,wr=wr,ur=ur,wz=wz,sr=sr),V_1)
#    
#    write_txt(mesh,rot_rate,'rot',path)
#    write_txt(mesh,total_rate,'tot',path)
#    write_txt(mesh,sr,'sr',path)
#
#    np.savetxt(path + '/ss_sign.txt', get_base_val(shear_stress_sign, d), fmt='%1.9f')
#    np.savetxt(path + '/wss.txt', get_base_val(wss, d), fmt='%1.9f')
#    np.savetxt(path + '/x_base.txt', x_base, fmt='%1.9f')
#    #
#    # #  Create velocity function.
#    (F2, tau_gamma, Re, VV, V,W) = fene_init(mesh, u, w, ur, uz, wr, wz, sub_domains, re, sr)
#    Ext=fene_solve(F2,tau_gamma,Re,sub_domains,re,sr,VV,V,W,steps)
#    base_ext=get_base_val(Ext,d)
#    np.savetxt(path + '/base_ext.txt', base_ext, fmt='%1.9f')



lc=0.025
#Res=[200,300,400]
# hs=[0.4,0.3,0.2]
hs=[0.5]
#Res=[500]
re=400
for i in range(len(hs)):
    h=hs[i]
    steps=int(np.floor(re*0.2))

    if not os.path.exists(path1):
        os.mkdir(path1)
    path = path1 + '/Reynolds_' + str(re)+'h-'+str(h)+'l1-'+str(l1)+'l2-'+str(l2)
    if not os.path.exists(path):
        os.mkdir(path)
    print(path)
    
    sub_domains, mesh = write_mesh(l1, l2, h, path, l,lc)  # write mesh and domains for solving
       # sub_domains, mesh = write_rect_mesh(l,100,100)  # write mesh and domains for solving
   
    ds = ds(subdomain_data=sub_domains)  # assign surface integration measure
    d, x_base, y_base = get_base_dofs(mesh, sub_domains, 0)  # get dofs of base and their z values    
    (max_w, meshsize,ss_1,ss_2,psi,shear_stress_sign,u, w ,ur,uz,wr,wz,shear_stress,sr) =NS(mesh,re,sub_domains)

    write_txt(mesh,u,'usol',path)
    write_txt(mesh,w,'wsol',path)

    V_1 = FunctionSpace(mesh, 'CG',  1)  # continous galerkin space for the velocity field.
    wss = interpolate(Expression('sqrt(ss_1* ss_1+ss_2 * ss_2)',degree=2,ss_1=ss_1,ss_2=ss_2),V_1)
    rot_rate= interpolate(Expression('uz-wr',degree=2,uz=uz,wr=wr,ur=ur,wz=wz),V_1)
    total_rate= interpolate(Expression('abs((uz-wr))-sr',degree=2,uz=uz,wr=wr,ur=ur,wz=wz,sr=sr),V_1)

    write_txt(mesh,rot_rate,'rot',path)
    write_txt(mesh,total_rate,'tot',path)
    write_txt(mesh,sr,'sr',path)

    np.savetxt(path + '/ss_sign.txt', get_base_val(shear_stress_sign, d), fmt='%1.9f')
    np.savetxt(path + '/wss.txt', get_base_val(wss, d), fmt='%1.9f')
    np.savetxt(path + '/x_base.txt', x_base, fmt='%1.9f')
    #
    # #  Create velocity function.
    (F2, tau_gamma, Re, VV, V,W) = fene_init(mesh, u, w, ur, uz, wr, wz, sub_domains, re, sr)
    Ext=fene_solve(F2,tau_gamma,Re,sub_domains,re,sr,VV,V,W,steps)
    base_ext=get_base_val(Ext,d)
    np.savetxt(path + '/base_ext.txt', base_ext, fmt='%1.9f')


