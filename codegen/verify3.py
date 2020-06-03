#   Copyright (c) 2014-2015, Sandia Corporation
#   All rights reserved.
#
#   This file is part of fast-matmul and is under the BSD 2-Clause License, 
#   which can be found in the LICENSE file in the root directory, or at 
#   http://opensource.org/licenses/BSD-2-Clause.

import sys
from collections import defaultdict
from fractions import Fraction

# holds a polynomial in epsilon over integers (could be changed to floats pretty easily)
# internal storage is a dict from powers of epsilon to coefficients
class Number:
    def __init__( self, arg1, arg2=None ):
        self.orig = arg1
        self.orig2 = None
        self.val = defaultdict(lambda:0)
        co = 0
        if type(arg1) == type( "1" ):
            # one argument, string constructor
            try:
                self.val[0] = Fraction(arg1)
            except:
                if arg1 == "x":
                    self.val[1] = 1
                    co = 1
                elif arg1 == "xi":
                    self.val[-1] = 1
                    co = -1
                elif arg1 == "-x":
                    self.val[1] = -1
                    co = 1
                elif arg1 == "-xi":
                    self.val[-1] = -1
                    co = -1
                elif arg1 == "x2":
                    self.val[2] = 1
                    co = 2
                elif arg1 == "-x2":
                    self.val[2] = -1
                    co = 2
                elif arg1 == "x2i":
                    self.val[-2] = 1
                    co = -2
                elif arg1 == "-x2i":
                    self.val[-2] = -1
                    co = -2
                elif arg1 == "-2x2":
                    self.val[2] = -2
                    co = 2
                elif arg1 == "x3":
                    self.val[3] = 1
                    co = 3
                elif arg1 == "-x3":
                    self.val[3] = -1
                    co = 3
                elif arg1 == "2x3":
                    self.val[3] = 2
                    co = 3
                elif arg1 == "-2x3":
                    self.val[3] = -2
                    co = 3
                elif arg1 == "x3i":
                    self.val[-3] = 1
                    co = -3
                elif arg1 == "-x3i":
                    self.val[-3] = -1
                    co = -3
                elif arg1 == "x4":
                    self.val[4] = 1
                    co = 4
                elif arg1 == "(1+-x3)":
                    self.val[0] = 1
                    self.val[3] = -1
                elif arg1 == "(1+x2)":
                    self.val[0] = 1
                    self.val[2] = 1
                elif arg1 == "(x+x2)":
                    self.val[1] = 1
                    self.val[2] = 1
                elif arg1 == "(x2+x3)":
                    self.val[2] = 1
                    self.val[3] = 1
                elif arg1 == "(x+-x2)":
                    self.val[1] = 1
                    self.val[2] = -1
                elif arg1 == "(-x2+-x3)":
                    self.val[2] = -1
                    self.val[3] = -1
                elif arg1 == "(-x+x2)":
                    self.val[1] = -1
                    self.val[2] = 1
                elif arg1 == "(x+-x4)":
                    self.val[1] = 1
                    self.val[4] = -1
                else:
                    print( "Not supported", arg1 )
                    raise NotImplemented
        elif type(arg1) == type(defaultdict(lambda:0)):
            # one argument, dict constructor
            self.val = arg1
        else:
            raise NotImplemented
        if arg2 != None:
            # two argument constructor
            #assert( type(arg1) == type(1) )
            #assert( type(arg2) == type(1) )
            #self.val[arg2]=arg1
            self.orig2 = arg2
            self.val[co] = Fraction(arg2)
    def __add__( self, other ):
        newVal = defaultdict(lambda:0)
        for p in self.val:
            newVal[p] += self.val[p]
        for p in other.val:
            newVal[p] += other.val[p]
        return Number(newVal)
    def __mul__( self, other ):
        newVal = defaultdict(lambda:0)
        for p1 in self.val:
            for p2 in other.val:
                newVal[p1+p2] += self.val[p1] * other.val[p2]
        return Number(newVal)
    def __eq__( self, other ):
        for p in self.val:
            if p <= 0:
                if self.val[p] != other.val[p]:
                    return False
        for p in other.val:
            if p <= 0:
                if self.val[p] != other.val[p]:
                    return False
        return True
    def __ne__( self, other ):
        return not (self == other)
    def exacteq( self, other ):
        for p in self.val:
            if self.val[p] != other.val[p]:
                return False
        for p in other.val:
            if self.val[p] != other.val[p]:
                return False
        return True

    # find minimum index corresponding to nonzero
    def findmin(self):
        minnum = 5
        for p in self.val:
            if self.val[p] != 0 and p < minnum:
                minnum = p
        return minnum

    def __repr__(self):
        # assume it is simple
        #if self.exacteq(Number("0")):
            #return "0"
        #if self.exacteq(Number("1")):
            #return "1"
        #if self.exacteq(Number("-1")):
            #return "-1"
        #if self.exacteq(Number("x")):
            #return "x"
        #if self.exacteq(Number("-x")):
            #return "-x"
        #if self.exacteq(Number("xi")):
            #return "xi"
        #if self.exacteq(Number("-xi")):
            #return "-xi"
        #if self.exacteq(Number("x2")):
            #return "x2"
        #if self.exacteq(Number("-x2")):
            #return "-x2"
        #if self.exacteq(Number("x2i")):
            #return "x2i"
        #if self.exacteq(Number("-x2i")):
            #return "-x2i"
        if self.orig2 != None:
            return str(self.orig2) + str(self.orig)
        return str(self.orig)


def main():        
	# the sizes; we check mxnxk matmul with q products, assuming row-major ordering.
	m = int(sys.argv[1])
	n = int(sys.argv[2])
	k = int(sys.argv[3])
	q = int(sys.argv[4])
	
	# read in U,V,W
	U = []
	V = []
	W = []
	#for i in range(m*n):
	#    U.append([])
	#    V.append([])
	#    W.append([])
	
	nnz = 0
	
	for j in range(m*n):
	    line = sys.stdin.readline().split()
	    while line[0] == "#":
	        line = sys.stdin.readline().split()
	    U.append([])
	    for i in range(q):
	        p = 0
	        coef = ''
	        while p < len(line[i]):
	            if line[i][p] != "x":
	                coef = coef + line[i][p]
	                p += 1
	            else:
	                break
	        if p == 0 or p == len(line[i]) or coef == '-' or coef == '(' or coef == '(-' or coef == '(1+-' or coef == '(1+':
	            U[j].append(Number(line[i]))
	        else:
	            line[i] = line[i][p:]
	            U[j].append(Number(line[i],coef))
	        if line[i] != "0":
	            nnz += 1
	
	for j in range(n*k):
	    line = sys.stdin.readline().split()
	    while line[0] == "#":
	        line = sys.stdin.readline().split()
	    V.append([])
	    for i in range(q):
	        p = 0
	        coef = ''
	        while p < len(line[i]):
	            if line[i][p] != "x":
	                coef = coef + line[i][p]
	                p += 1
	            else:
	                break
	        if p == 0 or p == len(line[i]) or coef == '-' or coef == '(' or coef == '(-' or coef == '(1+-' or coef == '(1+':
	            V[j].append(Number(line[i]))
	        else:
	            line[i] = line[i][p:]
	            V[j].append(Number(line[i],coef))
	        if line[i] != "0":
	            nnz += 1
	
	for j in range(m*k):
	    line = sys.stdin.readline().split()
	    while line[0] == "#":
	        line = sys.stdin.readline().split()
	    W.append([])
	    for i in range(q):
	        p = 0
	        coef = ''
	        while p < len(line[i]):
	            if line[i][p] != "x":
	                coef = coef + line[i][p]
	                p += 1
	            else:
	                break
	        if p == 0 or p == len(line[i]) or coef == '-' or coef == '(' or coef == '(-' or coef == '(1+-' or coef == '(1+':
	            W[j].append(Number(line[i]))
	        else:
	            line[i] = line[i][p:]
	            W[j].append(Number(line[i],coef))
	        if line[i] != "0":
	            nnz += 1
	
	print( "U" )
	for r in U:
	    print( r )
	print( "V" )
	for r in V:
	    print( r )
	print( "W" )
	for r in W:
	    print( r )
	print( "Naive number of additions:", nnz-len(U[0])-len(V[0])-len(W) )
	print("Number of nonzero entries:", nnz)
	
	def isConnected(M):
	    def add(s, r):
	        for i in range(len(r)):
	            if not r[i].exacteq(Number("0")):
	                s.add(i)
	    def overlap(s, r):
	        for i in range(len(r)):
	            if (not r[i].exacteq(Number("0"))) and i in s:
	                return True
	        return False
	    current = set()
	    add(current, M[0])
	    while True:
	        sOld = len(current)
	        for j in range(len(M)):
	            if overlap(current,M[j]):
	                add(current,M[j])
	        if len(current) == sOld:
	            break
	    if len(current) == len(M[0]):
	        return True
	    return False
	
	if not isConnected(U):
	    print( "U is disconnected" )
	if not isConnected(V):
	    print( "V is disconnected" )
	if not isConnected(W):
	    print( "W is disconnected" )

        
  # check algorithm and compute error parameter sigma
	sigma = 5
	#UVW = [0,0,0,0,0,0]
	for a in range(m):
	    for b in range(n):
	        # a and b describe row of U
	        rU = a*n+b
	        for c in range(n):
	            for d in range(k):
	                rV = c*k+d
	                for e in range(m):
	                    for f in range(k):
	                        rW = e*k+f
	                        # compute the contribution
	                        sum = Number("0")
	                        for i in range(q):
	                            sum += U[rU][i]*V[rV][i]*W[rW][i]
	                        if a == e and b == c and d == f:
	                            # should be a 1
	                            if sum != Number("1"):
	                                print( "Trouble at", a, b, c, d, e, f, "sum should be 1, is ", sum )
	                            else:
	                                # subtract 1 for sigma computation below
	                                sum += Number("-1")
	                        else:
	                            if sum != Number("0"):
	                                print( "Trouble at", a, b, c, d, e, f, "sum should be 0, is ", sum )
	                        # if largest error term has smaller exponent, update smallest exponent
	                        tmin = sum.findmin()
	                        if tmin < sigma:
	                            sigma = tmin
	                            UVW = [a, b, c, d, e, f]
	print ("sigma: ", sigma)
	print ("index of matrix when sigma satisfied: A[",UVW[0], "][",UVW[1], "], B[",UVW[2], "][",UVW[3], "], and C[",UVW[4], "][",UVW[5], "]")
  
	# calculate parameter phi
	phi = 0

	Tmin = 5
	UVWr = [0,0,0,0]
	for a in range(q):
	    minU, minV, minW = 5, 5, 5
	    Ub, Vb, Wb = 0, 0, 0
	    for b in range(m*n):
	        if U[b][a].findmin() < minU:
	            minU = U[b][a].findmin()
	            Ub = b
	    for b in range(n*k):
	        if V[b][a].findmin() < minV:
	            minV = V[b][a].findmin()
	            Vb = b
	    for b in range(k*m):
	        if W[b][a].findmin() < minW:
	            minW = W[b][a].findmin()
	            Wb = b
	    if minU+minV+minW < Tmin:
	        Tmin = minU+minV+minW
	        UVWr = [a, Ub, Vb, Wb]
	phi =  -Tmin
	print ("phi: ", phi)
	print ("index of UVW matrix when phi satisfied: U[", UVWr[1], "][", UVWr[0], "], V[", UVWr[2], "][", UVWr[0], "], and W[", UVWr[3], "][", UVWr[0], "]")
        

if __name__ == '__main__':
    main()
