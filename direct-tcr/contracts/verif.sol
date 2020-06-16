pragma solidity ^0.5.0;
pragma experimental ABIEncoderV2;
contract TCR {
  struct Point {
    uint x;
    uint y;
  }

  function P1() internal returns (Point memory) {
    return Point(1, 2);
  }

  function negate(Point memory p) internal returns (Point memory) {
    uint q = 21888242871839275222246405745257275088696311157297823662689037894645226208583;
    if (p.x == 0 && p.y == 0)
      return Point(0, 0);
    return Point(p.x, q - (p.y % q));
  }

function add(Point memory a, Point memory b) public view returns (Point memory r) {
 uint[4] memory input;
 input[0] = a.x;
 input[1] = a.y;
 input[2] = b.x;
 input[3] = b.y;

 assembly {
   if iszero(staticcall(gas, 0x06, input, 0x80, r, 0x40)) {
       revert(0,0)
   }
 }
 return r;
}
    
function mul(Point memory a, uint k) public view returns (Point memory p) {
 uint[3] memory input;
 input[0] = a.x;
 input[1] = a.y;
 input[2] = k;

 assembly {
   if iszero(staticcall(gas, 0x07, input, 0x60, p, 0x40)) {
       revert(0,0)
   }
 }
 return p;
}


  struct statement {
    Point g0;
    Point g1;
    Point h0;
    Point h1;
    Point y;
  }

  struct commitments {
    Point c0;
    Point c1;
    Point c2;
    Point c3;
    Point c4;
  }

  struct computed {
    Point A0;
    Point A1;
    Point A2;
    Point A3;
    Point A4;
    Point A5;
  }

  struct proof {
    Point R0;
    Point R1;
    Point R2;
    Point R3;
    Point S0;
    Point S1;
    Point T;
    uint d0;
    uint d1;
    uint u;
    uint v;
    uint w;
    uint a0;
    uint a1;
    uint a2;
  }

  function verify(proof memory p) public returns (uint) {
    Point memory g = Point(1, 2);
    statement memory s = statement(g, g, g, g, g);
    computed memory a = computed(g, g, g, g, g, g);
    // compute a0, a1, a2
    uint success;
    // g0 * d0 == R0 + A0 * a0
    if (mul(s.g0, p.d0).x == add(p.R0, (mul(a.A0, p.a0))).x) {
      success++;
    }
    // (h0 - y) * d0 == R1 + A3 * a0
    if (mul(add(s.h0, Point(s.y.x, -s.y.y)), p.d0).x == add(p.R1, (mul(a.A3, p.a0))).x) {
      success++;
    }
    // h1 * d1 == R2 + A4 * a1
    if (mul(s.h1, p.d1).x == add(p.R2, mul(a.A4, p.a1)).x) {
      success++;
    }
    // h1 * d1 == R2 + A4 * a1
    if (mul(s.h1, p.d1).x == add(p.R2, mul(a.A4, p.a1)).x) {
      success++;
    }
    // g0 * d1 == R3 + A5 * a1
    if (mul(s.g0, p.d1).x == add(p.R3, mul(a.A5, p.a1)).x) {
      success++;
    }
    // g1 * u + h0 * v == S0 + A1 * a2
    if (add(mul(s.g1, p.u), mul(s.h0, p.v)).x == add(p.S0, mul(a.A1, p.a2)).x) {
      success++;
    }
    // g1 * u + y * v == S1 + A2 * a2
    if (add(mul(s.g1, p.u), mul(s.y, p.v)).x == add(p.S1, mul(a.A2, p.a2)).x) {
      success++;
    }
    // h0 * w == T + A1 * (a2 - u) (do we have a scalar minus function? :/)
    if (mul(s.h0, p.w).x == add(p.T, mul(a.A1, (p.a2 - p.u))).x) {
      success++;
    }
    return success;
  }
}
