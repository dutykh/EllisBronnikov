# Matrices assembling for Ellis-Bronnikov phantom wormhole's quasi-normal modes computation
# Scalar perturbation case
MatrixAssembler := proc (
  d::integer, # d : number of digits used in computations
  n::integer, # n : number of Tchebyshev modes
  L::integer, # L : angular momentum
  c::numeric, # c : phantom mass
  p::string   # p : string containing the path where we save the assembled matrices
  )
  local L00::function, L01::function, L02::function, L10::function, L11::function, L12::function, L20::function, L21::function, L22::function, M0::Matrix, M1::Matrix, M2::Matrix, i::integer, j::integer, expr0::algebraic, expr1::algebraic, expr2::algebraic, xi::numeric, path::string, nstr::string:
  with(LinearAlgebra):
  Digits := d:
  # Definition of the 2nd order ODE coefficients:
  L00 := y -> -1/4*cos(1/2*Pi*y)^2*((-c^2 + 4)*cos(1/2*Pi*y)^2 + 2*c*sin(Pi*y) + 4*L*(L+1))/(-y^2 + 1)^2:
  L01 := y -> 2*cos(1/2*Pi*y)^2*(y + 1)*(-c*cos(1/2*Pi*y)^2*(Pi*y - Pi) + Pi*(y - 1)*sin(Pi*y))/(Pi^2*(-y^2 + 1)^3):
  L02 := y -> 4*cos(1/2*Pi*y)^4/(Pi^2*(-y^2 + 1)^2):
  L10 := y -> -1/(4*Pi^2*(-y^2 + 1)^4)*(-8*c*exp(c*Pi)*(-y + 1)^2*cos(1/2*Pi*y)^3*(-2*Pi*(y + 1)*sin(1/2*Pi*y) + (c*Pi*(y + 1) - 2)*cos(1/2*Pi*y)) + 4*Pi*(-y^2 + 1)^2*cos(1/2*Pi*y)^2*exp(-1/2*c*Pi*(y - 1))*(2*(Pi*c*y - 2)*cos(1/2*Pi*y)^2 + c*sin(Pi*y) - Pi*c*y) + 8*c*(y + 1)^2*cos(1/2*Pi*y)^3*(-2*Pi*(y - 1)*sin(1/2*Pi*y) + (c*Pi*(y - 1) - 2)*cos(1/2*Pi*y))):
  L11 := y -> 2*cos(1/2*Pi*y)^2*((y^2 - 1)*exp(-1/2*c*Pi*(y - 1))*(sin(Pi*y)*(Pi*c*y - 2) - 2*Pi*y) + 4*exp(c*Pi)*c*(-y + 1)*cos(1/2*Pi*y)^2 + 4*(y + 1)*c*cos(1/2*Pi*y)^2)/(Pi^2*(-y^2 + 1)^3):
  L12 := y -> 0:
  L20 := y -> -1/(4*Pi^2*(-y^2 + 1)^4)*(16*c^2*cos(1/2*Pi*y)^4*((y - 1)^2*exp(2*c*Pi) - 2*(y^2 - 1)*exp(c*Pi) + (y + 1)^2) + 16*c*(y^2 - 1)*(y + 1)*exp(-1/2*c*Pi*(y - 1))*cos(1/2*Pi*y)^2*((1/2*Pi*c*y - 1)*sin(Pi*y) - Pi*y) - 16*c*(y^2 - 1)*(y - 1)*exp(-1/2*c*Pi*(y - 3))*cos(1/2*Pi*y)^2*((1/2*Pi*c*y - 1)*sin(Pi*y) - Pi*y) - 4*(y - 1)^2*(y + 1)^2*exp(-c*Pi*(y - 1))*((Pi*c*y - 2)^2*cos(1/2*Pi*y)^4 - (Pi*c*y - 2)^2*cos(1/2*Pi*y)^2 + Pi*y*(Pi*c*y - 2)*sin(Pi*y) - Pi^2*(y^2 - 1))):
  L21 := y -> 0:
  L22 := y -> 0:
  F := y -> simplify(add(a[j]*ChebyshevT(j, y), j=0..n-1)):
  M0 := Matrix(n):
  M1 := Matrix(n):
  M2 := Matrix(n):
  for i from 1 to n do
    xi := cos((2.0*i-1.0)*Pi/(2.0*n)); # Chebyshev roots collocation points
    expr0 := evalf(L00(xi)*F(xi) + L01(xi)*subs(x=xi, diff(F(x),x)) + L02(xi)*subs(x=xi, diff(F(x),x$2))):
    expr1 := evalf(L10(xi)*F(xi) + L11(xi)*subs(x=xi, diff(F(x),x)) + L12(xi)*subs(x=xi, diff(F(x),x$2))):
    expr2 := evalf(L20(xi)*F(xi) + L21(xi)*subs(x=xi, diff(F(x),x)) + L22(xi)*subs(x=xi, diff(F(x),x$2))):
    for j from 1 to n do
      M0[i,j] := coeff(expr0, a[j-1]):
      M1[i,j] := coeff(expr1, a[j-1]):
      M2[i,j] := coeff(expr2, a[j-1]):
    end do:
  end do:
  # We finally export the data from Maple and save in files:
  path := cat(p, "/data/"):
  nstr := convert(n, string);
  ExportMatrix(cat(path, "M0_", nstr, ".mat"), M0, target=MATLAB, mode=ascii):
  ExportMatrix(cat(path, "M1_", nstr, ".mat"), M1, target=MATLAB, mode=ascii):
  ExportMatrix(cat(path, "M2_", nstr, ".mat"), M2, target=MATLAB, mode=ascii):
end proc: