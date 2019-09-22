pari.allocatemem(10^9)

import hashlib
babyjubjubr = 2736030358979909402780800718157159386076813972158567259200215660948447373041
babyjubjubq = 21888242871839275222246405745257275088548364400416034343698204186575808495617

Fq = GF(babyjubjubq)
Fr = GF(babyjubjubr)

a = Fq(168700)
d = Fq(168696)
base_x = Fq(16540640123574156134436876038791482806971768689494387082833631921987005038935)
base_y = Fq(20819045374670962167435360035096875258406992893633759881276124905556507972311)

# convert edwards to montgomery
A = 2*(a+d)/(a-d)
B = 4/(a-d)
base_u = (1+base_y)/(1-base_y)
base_v = base_u/base_x

# convert montgomery to weierstrass
wa = (3-A^2)/(3*B^2)
wb = (2*A^3-9*A)/(27*B^3)
x0, y0 = (base_u+A/3)/B, base_v/B

# twisted edwards curve arithmetic

def ed_add(P, Q):
  [x1, y1] = P
  [x2, y2] = Q
  x3 = (x1*y2+y1*x2)/(1+d*x1*x2*y1*y2)
  y3 = (y1*y2-a*x1*x2)/(1-d*x1*x2*y1*y2)
  return [x3, y3]

def ed_double(P):
  [x, y] = P
  x3 = (x*y+y*x)/(1+d*x*x*y*y)
  y3 = (y*y-a*x*x)/(1-d*x*x*y*y)
  return [x3, y3]

def ed_negate(P):
  [x, y] = P
  return [-x, y]

# montgomery to weierstrass
def m_to_w(P):
  [x, y] = P
  return [(x+A/3)/B, y/B]

# weierstrass to montgomery
def w_to_m(P):
  [x, y] = P
  return [B*x-(A/3), B*y]

# twisted edwards to montgomery
def e_to_m(P):
  [x, y] = P
  return [(1+y)/(1-y), (1+y)/(x-x*y)]

# montgomery to twisted edwards
def m_to_e(P):
  [x, y] = P
  return [x/y, (x-1)/(x+1)]

# weierstrass to twisted edwards
def w_to_e(P):
  return m_to_e(w_to_m(P))

# twisted edwards to weierstrass
def e_to_w(P):
  return m_to_w(e_to_m(P))

assert m_to_w([base_u, base_v]) == [x0, y0]
assert w_to_m([x0, y0]) == [base_u, base_v]
assert e_to_m([base_x, base_y]) == [base_u, base_v]
assert m_to_e([base_u, base_v]) == [base_x, base_y]
assert w_to_e([x0, y0]) == [base_x, base_y]
assert e_to_w([base_x, base_y]) == [x0, y0]

Ebjj = EllipticCurve(Fq, [wa, wb])
a = Ebjj.random_point()
b = Ebjj.random_point()

g = Ebjj(x0, y0)
g0 = Ebjj(10480227264716662452755474067665177302533350463766626583754719465282143096156,
8369223176222797522123539354678750082493875221914526607111902450577289923070)
g1 = Ebjj(11880246950178646564845070742914370409150886704747897303831567178374129810211,
11568178448482280388943253237660832271540088086331062617007232034122551287432)
h0 = Ebjj(3168195724794061035067437460306519652998969886362841281644734115844445886185,
8489771060550565563483479806652696626788502212157814537990951108069267515955)
h1 = Ebjj(12239960479043393909262572194446271981202432681542270714195285876953060316412,
16901143140203092244252603065645010058104514960352419130166756310402509981799)

# hard_coded_points = [g, g0, g1, h0, h1]
# for p in hard_coded_points:
#   print(w_to_e([p[0], p[1]]))

Scalars = GF(Ebjj.cardinality() / 8)
assert Ebjj.cardinality() / 8 == babyjubjubr

def random_value(prefix, i):
  return int(hashlib.sha256('%s%d' % (prefix, i)).hexdigest(), 16)

def prove_prime(w, s):
  [x0, x1, x2, x3, x4, x5, x6, b0, b1] = w
  [y, g0, g1, h0, h1, c0, c1, c2, c3, c4] = s

  A0 = c0
  A1 = c1
  A2 = c2
  A3 = c0 - c1 + c2
  A4 = c1 - c2 + c4
  A5 = c3
  A6 = c1 - c2 + c3 + c4

  r0 = Integer(Scalars.random_element())
  r1 = Integer(Scalars.random_element())
  r2 = Integer(Scalars.random_element())
  r3 = Integer(Scalars.random_element())

  s0 = Integer(Scalars.random_element())
  s1 = Integer(Scalars.random_element())
  s2 = Integer(Scalars.random_element())
  s3 = Integer(Scalars.random_element())
  s4 = Integer(Scalars.random_element())
  s5 = Integer(Scalars.random_element())

  R0 = r0*g0
  R1 = r1*(g0 - h0 + y)
  R2 = r2*h1
  R3 = r3*(g0 + h1)
  S0 = s0*g1 + s1*h0
  T0 = (s0*b0)*g1 + s2*h0
  S1 = s3*g1 + s4*y
  T1 = (s3*b1)*g1 + s5*y

  a0 = random_value('1', 1)
  a1 = random_value('1', 1)
  a2 = random_value('1', 1)
  a3 = random_value('1', 1)
  a4 = random_value('1', 1)
  a5 = random_value('1', 1)
  a6 = random_value('1', 1)

  d0 = r0 + a0*x0 + a5*x5
  u0 = s0 + a1*b0
  v0 = s1 + a1*x1
  w0 = s2 + x1*(a1 - u0)
  u1 = s3 + a2*b1
  v1 = s4 + a2*x2
  w1 = s5 + x2*(a2 - u1)
  d1 = r1 + a3*x3
  d2 = r2 + a4*x4
  d3 = r3 + a6*x6

  assert d0*g0 == R0 + a0*A0 + a5*A5
  assert u0*g1 + v0*h0 == S0 + a1*A1
  assert w0*h0 == T0 + (a1 - u0)*A1
  assert u1*g1 + v1*y == S1 + a2*A2
  assert w1*y == T1 + (a2 - u1)*A2
  assert d1*(g0 - h0 + y) == R1 + a3*A3
  assert d2*h1 == R2 + a4*A4
  assert d3*(g0 + h1) == R3 + a6*A6

  print "a"
  print(a[0])

  fstmsg = [R0, R1, R2, R3, S0, S1, T0, T1]
  for msg in fstmsg:
    emsg = w_to_e([msg[0], msg[1]])
    print(emsg[0])
    print(emsg[1])

  sndmsg = [d0, d1, d2, d3, u0, v0, w0, u1, v1, w1, (a1 - u0), (a2 - u1)]
  for msg in sndmsg:
    print(Fr(msg))

Ebjj = EllipticCurve(Fq, [wa, wb])
a = Ebjj.random_point()
b = Ebjj.random_point()
Scalars = GF(Ebjj.cardinality() / 8)
x = random_value('hello tcr', 1)
x0 = x
x1 = x
x2 = x
x3 = x
x4 = x
x5 = x
x6 = x
xx = random_value('hello y', 1)
y = xx*g
# print "y"
# print(w_to_e([y[0], y[1]]))

b = Integer(1)
b0 = b
b1 = b
c0 = x*g0
c1 = b*g1+x*h0
c2 = b*g1+x*y
c3 = x*g0
c4 = x*(y-h0)+x*h1

# commitments = [c0, c1, c2, c3, c4]
# for c in commitments:
#   print(w_to_e([c[0], c[1]]))

"""
# A0 = x0 g0 = c0
A0 = x0*g0
# A1 = b0 g1 + x1 h0 = c1
A1 = b0*g1+x1*h0
# A2 = b1 g1 + x2 y = c2
A2 = b1*g1+x2*y
# A3 = x3 (g0 - h0 + y) = c0 - c1 + c2
A3 = x3*(g0-h0+y)
# A4 = x4 h1 = c0 - c2 + c4
A4 = x4*h1
# A5 = x5 g0 = c3
A5 = x5*g0
# A6 = x6 (g0 + h1) = c1 - c2 + c3 + c4
A6 = x6*(g0+h1)
"""
s = [y, g0, g1, h0, h1, c0, c1, c2, c3, c4]
w = [x, x, x, x, x, x, x, b, b]

prove_prime(w, s)

# print(hashlib.sha256(c0 + c1 + c2 + c3))
