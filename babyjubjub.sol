// This file is MIT Licensed.
// PRs welcome @ github.com/rbkhmrcr/tcr

pragma solidity ^0.5.0;

library babyjubjub {
    // Points of the twisted edwards curve 168700x^2 + y^2 = 1 + 168696x^2y^2
    // uint aliases uint256
    struct bjjPoint {
        uint x;
        uint y;
    }
    // A = 168700
    uint public constant A = 0x292FC;
    // D = 168696 
    uint public constant D = 0x292F8;
    // field modulus = 21888242871839275222246405745257275088548364400416034343698204186575808495617
    uint public constant modulus = 0x30644E72E131A029B85045B68181585D2833E84879B9709143E1F593F0000001;


     // p + q = ((x1y2 + y1x2) / (1 + dx1x2y1y2), (y1y2 - ax1x2) / (1 - dx1x2y1y2))
    function pointadd(bjjPoint memory p, bjjPoint memory q) internal view returns (bjjPoint memory r) {
        if (p.x == 0 && p.y == 0) {
            return q;
        }

        if (q.x == 0 && q.y == 0) {
            return p;
        }

        uint dprod = mulmod(D, mulmod(mulmod(p.x, q.x, modulus), mulmod(p.y, q.y, modulus), modulus), modulus);
        // x3 = (x1y2 + y1x2) / (1 + dx1x2y1y2)
        uint x3num = addmod(mulmod(p.x, q.y, modulus), mulmod(p.y, q.x, modulus), modulus);
        r.x = mulmod(x3num, invmod(addmod(1, dprod, modulus)), modulus);
        
        // y3 = (y1y2 - ax1x2) / (1 - dx1x2y1y2)
        uint ax1x2 = mulmod(A, mulmod(p.x, q.x, modulus), modulus);
        uint y3num = submod(mulmod(p.y, q.y, modulus), ax1x2, modulus);
        r.y = mulmod(y3num, invmod(submod(1, dprod, modulus)), modulus);
        return r;
    }

    function pointdbl(bjjPoint memory p) internal view returns (bjjPoint memory r) {
        return pointadd(p, p);
    }

    function scalarmul(bjjPoint memory p, uint k) internal view returns (bjjPoint memory r) {
        uint rem = k;
        bjjPoint memory pp = bjjPoint(p.x, p.y);
        r = bjjPoint(0, 0);

        while (rem != 0) {
            if ((rem & 1) != 0) {
                r = pointadd(r, pp);
            }
            pp = pointdbl(r);
            rem = rem / 2;
        }
        return r;
    }

    function isoncurve(bjjPoint memory p) internal pure returns (bool) {
        uint x2 = mulmod(p.x, p.x, modulus);
        uint y2 = mulmod(p.y, p.y, modulus);
        uint lhs = addmod(mulmod(A, x2, modulus), y2, modulus);
        uint rhs = addmod(1, mulmod(mulmod(D, x2, modulus), y2, modulus), modulus);
        return submod(lhs, rhs, modulus) == 0;
    }

    function submod(uint a, uint b, uint mod) internal pure returns (uint) {
        uint pa = a;
        if (a <= b) {
            pa += mod;
        }
        return addmod(pa - b, 0, mod);
    }

    function invmod(uint a) internal view returns (uint) {
        return expmod(a, modulus - 2, modulus);
    }

    function expmod(uint base, uint e, uint m) internal view returns (uint o) {
  
        assembly {
            // define pointer
            let p := mload(0x40)
            // store data assembly-favouring ways
            mstore(p, 0x20)             // Length of Base
            mstore(add(p, 0x20), 0x20)  // Length of Exponent
            mstore(add(p, 0x40), 0x20)  // Length of Modulus
            mstore(add(p, 0x60), base)  // Base
            mstore(add(p, 0x80), e)     // Exponent
            mstore(add(p, 0xa0), m)     // Modulus
            // call modexp precompile! -- old school gas handling
            let success := staticcall(gas, 0x05, p, 0xc0, p, 0x20)
            // gas fiddling
            switch success case 0 {
                revert(0, 0)
            }
            // data
            o := mload(p)
        }
    }
}

