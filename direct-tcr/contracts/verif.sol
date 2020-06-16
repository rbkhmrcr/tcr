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

  function negate(Point memory p) public view returns (Point memory) {
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
  
  struct pok {
      Point R;
      uint d;
  }
  
  uint[] amts;
    
  function deposit(uint amt, pok memory p, Point memory g, Point memory y) public returns (Point memory) {
        uint256 c = uint256(sha256(abi.encode(p.R.x)));
        Point memory ycr = add(mul(y, c), p.R);
        Point memory rhs = negate(mul(g, p.d));
        amts.push(amt);
        return add(ycr, rhs);
  }
  
   
  function update(uint amt, uint i, pok memory p, Point memory g, Point memory y) public returns (Point memory) {
        uint256 c = uint256(sha256(abi.encode(p.R.x)));
        Point memory ycr = add(mul(y, c), p.R);
        Point memory rhs = negate(mul(g, p.d));
        amts[i] = amt;
        return add(ycr, rhs);
  }
  
  function verify_vote1(pok memory p, Point memory g, Point memory y) public returns (Point memory) {
        uint256 c = uint256(sha256(abi.encode(p.R.x)));
        Point memory ycr = add(mul(y, c), p.R);
        Point memory rhs = negate(mul(g, p.d));
        return add(ycr, rhs);
  }

  function verify_vote2(proof memory p, computed memory com) public returns (uint) {
    Point memory g = Point(1, 2);
    Point memory g0 = Point(1368015179489954701390400359078579693043519447331113978918064868415326638035,9918110051302171585080402603319702774565515993150576347155970296011118125764);
    Point memory g1 = Point(3010198690406615200373504922352659861758983907867017329644089018310584441462,4027184618003122424972590350825261965929648733675738730716654005365300998076);
    Point memory h0 = Point(3932705576657793550893430333273221375907985235130430286685735064194643946083,18813763293032256545937756946359266117037834559191913266454084342712532869153);
    Point memory h1 = Point(10835225521862395592687560951453385602895512958032257955899877380493200080708,2623520004791921319615054428233368525468155544765295675952919303096698181037);
    statement memory s = statement(g, g0, g1, h0, h1);
    // compute a0, a1, a2
    uint success;
    // g0 * d0 == R0 + A0 * a0
    if (mul(s.g0, p.d0).x == add(p.R0, (mul(com.A0, p.a0))).x) {
      success++;
    }
    // (h0 - y) * d0 == R1 + A3 * a0
    if (mul(add(s.h0, Point(s.y.x, -s.y.y)), p.d0).x == add(p.R1, (mul(com.A3, p.a0))).x) {
      success++;
    }
    // h1 * d1 == R2 + A4 * a1
    if (mul(s.h1, p.d1).x == add(p.R2, mul(com.A4, p.a1)).x) {
      success++;
    }
    // h1 * d1 == R2 + A4 * a1
    if (mul(s.h1, p.d1).x == add(p.R2, mul(com.A4, p.a1)).x) {
      success++;
    }
    // g0 * d1 == R3 + A5 * a1
    if (mul(s.g0, p.d1).x == add(p.R3, mul(com.A5, p.a1)).x) {
      success++;
    }
    // g1 * u + h0 * v == S0 + A1 * a2
    if (add(mul(s.g1, p.u), mul(s.h0, p.v)).x == add(p.S0, mul(com.A1, p.a2)).x) {
      success++;
    }
    // g1 * u + y * v == S1 + A2 * a2
    if (add(mul(s.g1, p.u), mul(s.y, p.v)).x == add(p.S1, mul(com.A2, p.a2)).x) {
      success++;
    }
    // h0 * w == T + A1 * (a2 - u) (do we have a scalar minus function? :/)
    if (mul(s.h0, p.w).x == add(p.T, mul(com.A1, (p.a2 - p.u))).x) {
      success++;
    }
    return success;
  }
}
