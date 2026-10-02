"""Exact bounded point search in selected division subfields, with modular sieving.

sage search_box.py --radius 3 --seconds 300
An unsuccessful bounded search gives no upper bound on the rank.
"""
import argparse
import json
import time
import numpy as np
from sage.all import EllipticCurve, GF, NumberField, PolynomialRing, QQ, ZZ, gcd, prime_range

p = argparse.ArgumentParser(description=__doc__)
p.add_argument('--radius', type=int, default=3)
p.add_argument('--seconds', type=int, default=300)
p.add_argument('--basis', choices=['integral','reduced'], default='reduced')
p.add_argument('--case', choices=['7a','5c','5d'], default='7a')
args = p.parse_args()
assert 1 <= args.radius <= 5
started = time.monotonic()
R = PolynomialRing(QQ, 'x'); x = R.gen()
polynomials = {'7a':x**7-18,'5c':x**10-18*x**5-324,'5d':x**10-54*x**5+324}
f = polynomials[args.case]
ainv = [1,0,0,-1,-1] if args.case == '7a' else [1,0,0,-3,-3]
K = NumberField(f, 'a'); a = K.gen()
E = EllipticCurve(K, ainv)
degree = K.degree()
assert (2*args.radius+1)**degree <= 20000000
basis = list(K.reduced_basis(prec=160) if args.basis == 'reduced' else K.integral_basis())
print('CASE',args.case,'FIELD',f,'\nBASIS_TYPE',args.basis,'BASIS', basis, flush=True)
# Prime-to-residue-characteristic torsion injects under good reduction.
good = []
for ell in prime_range(5,200):
    if EllipticCurve(QQ,ainv).discriminant() % ell == 0:
        continue
    e = EllipticCurve(GF(ell),ainv)
    if e.discriminant() and f.change_ring(GF(ell)).roots():
        good.append((ell,e.cardinality()))
    if len(good)==3:
        break
torsion_bound = gcd([n for ell,n in good])
for ell,n in good:
    v = min(ZZ(m).valuation(ell) for q,m in good if q!=ell)
    torsion_bound *= ell**max(0,v-ZZ(torsion_bound).valuation(ell))
print('TORSION_BOUND',torsion_bound,'FROM_REDUCTIONS',good,flush=True)
width = 2*args.radius+1
count = width**degree
grid = np.indices((width,)*degree, dtype=np.int16).reshape(degree,-1).T - args.radius
primes = list(prime_range(11,200))
conditions = []
for ell in primes:
    F = GF(ell)
    squares = np.array([bool(F(z).is_square()) for z in range(ell)])
    if any(any(ZZ(c.denominator()) % ell == 0 for c in b.list()) for b in basis):
        continue
    for r, _ in f.change_ring(F).roots():
        weights = np.array([int(R(b.polynomial()).change_ring(F)(r)) for b in basis],dtype=np.int64)
        residues = np.asarray((grid @ weights) % ell, dtype=np.int16)
        conditions.append((ell, squares, residues))
    if len(conditions) >= 16:
        break
print('RADIUS',args.radius,'NUMERATORS',count,'SIEVE_CONDITIONS',len(conditions),flush=True)
tested = 0
for denominator in [1,2,3,4,5,6,7,8,9,10,12,14,21,28]:
    live = np.arange(count)
    for ell, squares, residues in conditions:
        if denominator % ell == 0:
            continue
        z = (residues[live].astype(np.int64)*pow(denominator,-1,ell)) % ell
        fz = (4*z**3+z*z+4*ainv[3]*z+4*ainv[4]) % ell
        live = live[squares[fz]]
    print('DENOMINATOR',denominator,'SURVIVORS',len(live),flush=True)
    for index in live:
        if time.monotonic()-started > args.seconds:
            print('TIME_LIMIT exact square tests',tested,flush=True)
            raise SystemExit(0)
        xx = sum(int(c)*b for c,b in zip(grid[index],basis))/denominator
        d = 4*xx**3+xx**2+4*ainv[3]*xx+4*ainv[4]
        tested += 1
        if not d.is_square():
            continue
        P = E(xx,(-xx+d.sqrt())/2)
        if torsion_bound*P != E(0):
            print('PROVED_INFINITE_ORDER',P,flush=True)
            print('COORDINATES_JSON',json.dumps([[str(c) for c in q.list()] for q in [P[0],P[1]]]),flush=True)
            print('CERTIFICATE: the torsion bound times P is nonzero.',flush=True)
            raise SystemExit(0)
    print('ELAPSED',time.monotonic()-started,'EXACT_TESTS',tested,flush=True)
print('COMPLETED: no infinite-order point in the stated finite box.',flush=True)
