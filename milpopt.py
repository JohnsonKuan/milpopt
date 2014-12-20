"""

MILPOpt - Python Package for Mixed Integer Linear Programming (MILP) Optimization
==========================================================================================

Author: Johnson Kuan
Version: 1.0.0
Created 12/19/2014

Algorithms:
  LP: Simplex
  ILP: Gomory-Chvatal Cutting-Plane

"""

import math, random

#transpose
def T(x):
  return list(map(lambda *x: list(x), *x))

#divide
def D(a, b):
  try:
    a / b
  except ZeroDivisionError:
    return float('inf')
  else:
    return a / b

class milpopt:

  def __init__(self, c, A, b, t = 10**(-10), I = ()):
    self.c = c
    self.A = A
    self.b = b
    #floating point tolerance for integer constraints
    self.t = t
    self.I = I

  #validate inputs
  def validate(self):

    #lists
    if not isinstance(self.c, list) \
    or not isinstance(self.A, list) \
    or not isinstance(self.b, list):
      return 'Inputs must be lists.'
      quit()
    #dimensions
    elif len(self.c) != len(self.A[0]):
      return 'Dimension of objective coefficients (%s) and matrix A column coefficients (%s) do not match.' % (len(self.c), len(self.A[0]))
      quit()
    elif len(self.b) != len(T(self.A)[0]):
      return 'Dimension of basic coefficients (%s) and matrix A row coefficients (%s) do not match.' % (len(self.c), len(T(self.A)[0]))
      quit()
    #integer constraint
    elif len(self.c) != len(self.I):
      return 'Dimension of objective coefficients (%s) and integer constraint tuple (%s) do not match.' % (len(self.I), len(self.c))
      quit()
    else:
      return True

  #convert canonical form to dictionary
  def dict(self):

    #validate inputs
    if milpopt.validate(self) != True:
      return milpopt.validate(self)
    else:
      NBI = [i for i in range(1, len(self.c) + 1)]
      BI = [i for i in range(len(self.c) + 1, len(self.c) + len(self.b) + 1)]
      B = self.b
      A = [list(map(lambda x: (-1) * x, a)) for a in self.A]
      Z = [0] + self.c
      return (BI, NBI, B, A, Z)

  #lp simplex
  def lp_opt(self, BI, NBI, B, A, Z):

    #choose entering variable
    I = [NBI[i] for i, j in enumerate(Z[1:]) if j > self.t]

    #quit if no entering variable
    try:
      type(max(I))
    except ValueError:
      print('no entering variable')
      return (BI, NBI, B, A, Z)
      quit()

    #bland's rule
    #index of entering variable
    ev = min(I)
    evi = NBI.index(ev)

    #choose leaving variable
    I = [True for j in A if j[evi] < (-1) * self.t]

    #quit if unbounded
    try:
      type(max(I))
    except ValueError:
      print('unbounded')
      quit()

    I = [((-1) * D(B[i], j[evi])) for i, j in enumerate(A)]
    T = [j for j in I if j != float('inf') and j != float('-inf') and j == j]

    temp = max(T) + 1

    for i, j in enumerate(I):
      if j < temp and j >= (-1) * self.t and A[i][evi] < (-1) * self.t:
        temp = j
        lv = BI[i]
        lvi = i
      elif j == temp:
        #bland's rule
        lv = min(lv,BI[i])
        lvi = BI.index(lv)

    #rearrange dictionary
    R = A[lvi][:evi] + [-1] + A[lvi][evi + 1:]
    R = [(-1) * (j / A[lvi][evi]) for j in R]

    #update basic
    B = B[:lvi] + [(-1) * B[lvi] / A[lvi][evi]] + B[lvi + 1:]
    for i in range(len(B)):
      if i != lvi:
        B[i] += (A[i][evi] * B[lvi])

    #update A
    A[lvi] = R
    for i in range(len(A)):
      if i != lvi:
        I1 = [A[i][evi] * j for j in R]
        I2 = A[i][:evi] + [0] + A[i][evi + 1:]
        A[i] = [sum(x) for x in zip(I1, I2)]

    #update basic indices
    BI = BI[:lvi] + [ev] + BI[lvi + 1:]

    #update non-basic indices
    NBI = NBI[:evi] + [lv] + NBI[evi + 1:]

    #update objective value and coefficients
    I1 = [Z[1:][evi] * j for j in [B[lvi]] + R]
    I2 = [Z[0]] + Z[1:][:evi] + [0] + Z[1:][evi + 1:]
    Z = [sum(x) for x in zip(I1,I2)]

    if max(Z[1:]) <= self.t:
      #end simplex optimization phase and return final dictionary
      return (BI, NBI, B, A, Z)
    else:
      return milpopt.lp_opt(self, BI, NBI, B, A, Z)

  #lp initialize
  def lp_init(self, BI, NBI, B, A, Z):

    #check sign of basic coefficients
    I = [(0,1)[j < (-1) * self.t] for j in B]
    if 1 not in I:
      print('Initial primal feasible')
      return (BI, NBI, B, A, Z)
      quit()

    #convert to dual problem
    B_ = [1.0] *(len(Z) - 1)
    A_ = [list(map(lambda x: (-1) * x, a)) for a in T(A)]
    Z_ = [(-1.0) * j for j in [Z[0]] + B]
    NBI_ = BI
    BI_ = NBI

    #optimize dual problem
    sol = milpopt.lp_opt(self, BI_, NBI_, B_, A_, Z_)

    #convert back to primal
    B_ = [(-1.0) * j for j in sol[4][1:]]
    A_ = [list(map(lambda x: (-1) * x, a)) for a in T(sol[3])]
    Z_ = [(-1.0) * j for j in [sol[4][0]] + sol[2]]
    NBI_ = sol[0]
    BI_ = sol[1]

    #restore original objective
    I = [0.0] * (len(NBI) + 1)

    for i, j in enumerate(NBI):
      if j in BI_:
        k = BI_.index(j)
        I1 = [Z[1:][i] * j for j in [B_[k]] + A_[k]]
        I = [sum(x) for x in zip(I, I1)]

    for i, j in enumerate(NBI):
      if j not in BI_:
        k = NBI_.index(j)
        I1 = [0] + [(0, Z[1:][i])[n == k] for n in range(len(NBI_))]
        I = [sum(x) for x in zip(I, I1)]

    Z_ = I

    return (BI_, NBI_, B_, A_, Z_)

  #ilp gomory-chvatal cutting-plane
  def ilp_cut_plane(self, BI, NBI, B, A, Z):

    #check for integral solution
    I = [BI[i] for i, j in enumerate(B) if (j - math.floor(j)) > self.t and (math.ceil(j) - j) > self.t and BI[i] <= len(self.I)]
    try:
      type(max(I))
    except ValueError:
      return (BI, NBI, B, A, Z)
      quit()

    #start ilp
    #identify row with fractional constant
    I = [BI[i] for i, j in enumerate(B) if (j - math.floor(j)) > self.t and (math.ceil(j) - j) > self.t and BI[i] <= len(self.I)]

    #add cuts
    for b in I:
      if self.I[b - 1] == 1:
        bi = BI.index(b)

        #derive cutting plane
        B_ = (-1.0) * (B[bi] - math.floor(B[bi]))
        A_ = [((-1.0) * j) - math.floor(((-1.0) * j)) for j in A[bi]]

        #add cutting plane constraint
        B = B + [B_]
        A.append(A_)
        BI = BI + [len(BI + NBI) + 1]

    #convert to dual problem
    B_ = [(-1.0) * j for j in Z[1:]]
    A_ = [list(map(lambda x: (-1) * x, a)) for a in T(A)]
    Z_ = [(-1.0) * j for j in [Z[0]] + B]
    NBI_ = BI
    BI_ = NBI

    #optimize
    BI, NBI, B, A, Z = milpopt.lp_opt(self, BI_, NBI_, B_, A_, Z_)

    #convert to primal problem
    B_ = [(-1.0) * j for j in Z[1:]]
    A_ = [list(map(lambda x: (-1) * x, a)) for a in T(A)]
    Z_ = [(-1.0) * j for j in [Z[0]] + B]
    NBI_ = BI
    BI_ = NBI

    return milpopt.ilp_cut_plane(self, BI_, NBI_, B_, A_, Z_)

  #solve
  def solver(self):

    #primal dictionary
    BI, NBI, B, A, Z = milpopt.dict(self)

    #initialize
    BI, NBI, B, A, Z = milpopt.lp_init(self, BI, NBI, B, A, Z)

    if 1 in self.I:

      BI, NBI, B, A, Z = milpopt.lp_opt(self, BI, NBI, B, A, Z)

      return milpopt.ilp_cut_plane(self, BI, NBI, B, A, Z)

    else:

      return milpopt.lp_opt(self, BI, NBI, B, A, Z)
