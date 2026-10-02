"""Exact field checks for Coates--Sujatha Examples 4.7 and 4.8.

Run with SageMath: sage verify_fields.py
This proves the four field inclusions stated in the report. It does not
compute ranks over the infinite towers.
"""
from sage.all import EllipticCurve, NumberField, PolynomialRing, QQ, pari, version

R = PolynomialRing(QQ, 'x')
x = R.gen()
cases = {
    5: ([1,0,0,-3,-3], [1,0,0,-28,272],
        [x**5-x**4+4*x**3+4*x**2-x+13, x**5+5*x**3+5*x-1]),
    7: ([1,0,0,-1,-1], [1,0,0,-141,657],
        [x**7-18, x**7-14*x**4+7*x**3+28*x**2+21*x+8]),
}
print(version(), flush=True)
for p, (ainv, binv, polys) in cases.items():
    E, Ep = EllipticCurve(ainv), EllipticCurve(binv)
    assert any(phi.codomain().is_isomorphic(Ep) for phi in E.isogenies_prime_degree(p))
    print('\nPRIME', p, 'CURVES', E.cremona_label(), Ep.cremona_label(), flush=True)
    print('VERIFIED rational isogeny of degree', p, flush=True)
    for idx, g in enumerate(polys):
        A = [E, Ep][idx]
        factors = A.division_polynomial(p).factor()
        if idx == 0:
            f = next(f for f, e in factors if f.degree() == p)
            h, b = f, x
        else:
            f = max((f for f,e in factors), key=lambda f:f.degree())
            subfields = pari(f).nfsubfields(p)
            assert len(subfields) == 1
            h, b = map(R, subfields[0])
        assert h(b) % f == 0
        iso = pari(h).nfisisom(g)
        assert iso != 0
        for c in iso:
            assert h(R(c)) % g == 0
        K = NumberField(g, 'a')
        assert g.is_irreducible()
        gal = pari(g).polgalois()
        assert int(gal[0]) == p*(p-1)
        print('CASE', str(p)+'ab'[idx], flush=True)
        print('POLYNOMIAL', g, flush=True)
        print('DIVISION_FACTOR', f, flush=True)
        print('SUBFIELD_POLYNOMIAL', h, flush=True)
        print('SUBFIELD_EMBEDDING', b, flush=True)
        print('ISOMORPHISM_ROOTS', iso, flush=True)
        print('FIELD_DISCRIMINANT', K.discriminant(), flush=True)
        print('GALOIS_GROUP', gal, flush=True)
        print('CONTAINED_IN original E[%s] division field' % (p**(idx+1)), flush=True)
g5 = cases[5][2][0]
pure = x**5-18
isopure = pari(g5).nfisisom(pure)
assert isopure != 0
for c in isopure:
    assert g5(R(c)) % pure == 0
print('\nPURE_QUINTIC isomorphism roots', isopure, flush=True)

# Let beta^5 = u, where u^2-u-1=0. Then b=beta-beta^(-1)
# satisfies the polynomial of case 5b. If zeta is a primitive fifth root,
# b_j=zeta^j*beta-zeta^(-j)/beta are all five roots of that polynomial.
# Thus beta=(b_1-zeta^(-1)*b_0)/(zeta-zeta^(-1)) lies in its splitting
# field after adjoining zeta, hence in Q(E'[5]) and Q(E[25]).
B = NumberField(x**10-x**5-1, 'beta')
beta = B.gen()
assert cases[5][2][1](beta-1/beta) == 0
assert (beta**5)**2-beta**5-1 == 0
print('VERIFIED beta-beta^(-1) is a root of case 5b', flush=True)
for k, g in [(1, x**10-18*x**5-324), (2, x**10-54*x**5+324)]:
    assert g.is_irreducible()
    # c^5=18*beta^(5k), so check its quadratic polynomial exactly.
    t = 18*beta**(5*k)
    assert t**2-18*([1,3][k-1])*t+(-1)**k*324 == 0
    print('MIXED_FIELD', g, 'CONTAINED_IN E[25] division field', flush=True)
    print('FIELD_DISCRIMINANT', pari(g).nfdisc(), flush=True)
print('\nCOMPLETED: all field relations checked exactly.', flush=True)
