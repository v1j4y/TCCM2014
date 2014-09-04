#!/usr/bin/env python
# parse mathematica script output on a terminal with
# sympy pprint to get awesome termiinal output
# 
# Presently supports printing 
# 1. Matrices (1 or many)
# 2. Equations
#
# TODO:
#   1. add support for plots maybe
#   2. better printing of equations 
#   3. avoid using eval()

from sympy import *
import numpy as np
import scipy.linalg
import re
from StringIO import StringIO


def main():
    file=open("inp","r")
    count=0
    for line in file.readlines():

        # read output file(inp) for different types of input
        #   1. check if multiple matrices
        if "{" in line:
            line=re.sub("{",'',line)
            line=re.sub(",",'\t',line)
            line=re.sub("}",'\n',line)
            c=StringIO(line)
            mat=np.loadtxt(c)
            A=Matrix(mat)
            pprint(A)
            count=1
        
        #   2. check if string (equation)
        num=line.split('\t',1)
        num=''.join(num[0])
        try:
            num=float(num)
        except ValueError:
            a,b,c, d, e, p, t, i=symbols('a b c d e p t i')
            x, y, z, h,m, Vo, k, r=symbols('x y z h m Vo k r')
            count=1
            func=line
            func=func.lower()
            func=re.sub('\^','**',func)
            func=re.sub('\[','(',func)
            func=re.sub('\]',')',func)
            func=re.sub('\(i\*','(I*',func)
            func=re.sub('arctan','atan',func)
            func=re.sub('arccos','acos',func)
            func=re.sub('arcsin','asin',func)
            pprint(eval(func))

        #   3. if nothing then assume a single matrix
    if count !=1:
            mat=np.loadtxt("inp")
            A=Matrix(mat)
            pprint(A)

if __name__ == "__main__":
    main()
