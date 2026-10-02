from sage.all import *
import time
R=PolynomialRing(QQ,"x"); x=R.gen()
K=NumberField(x**7-18,"a")
E=EllipticCurve(K,[1,0,0,-1,-1])
print("GRH-assisted exploratory Simon descent; only explicit points will be certified",flush=True)
r=E.simon_two_descent(verbose=2,limbigprime=20,lim1=10,lim3=10,limtriv=3,maxprob=20)
print("SIMON_OUTPUT",r,flush=True)
print("COMPLETED",flush=True)
