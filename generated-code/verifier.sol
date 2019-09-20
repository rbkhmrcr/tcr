// This file is LGPL3 Licensed

/**
 * @title Elliptic curve operations on twist points for alt_bn128
 * @author Mustafa Al-Bassam (mus@musalbas.com)
 * @dev Homepage: https://github.com/musalbas/solidity-BN256G2
 */

library BN256G2 {
    uint256 internal constant FIELD_MODULUS = 0x30644e72e131a029b85045b68181585d97816a916871ca8d3c208c16d87cfd47;
    uint256 internal constant TWISTBX = 0x2b149d40ceb8aaae81be18991be06ac3b5b4c5e559dbefa33267e6dc24a138e5;
    uint256 internal constant TWISTBY = 0x9713b03af0fed4cd2cafadeed8fdf4a74fa084e52d1852e4a2bd0685c315d2;
    uint internal constant PTXX = 0;
    uint internal constant PTXY = 1;
    uint internal constant PTYX = 2;
    uint internal constant PTYY = 3;
    uint internal constant PTZX = 4;
    uint internal constant PTZY = 5;

    /**
     * @notice Add two twist points
     * @param pt1xx Coefficient 1 of x on point 1
     * @param pt1xy Coefficient 2 of x on point 1
     * @param pt1yx Coefficient 1 of y on point 1
     * @param pt1yy Coefficient 2 of y on point 1
     * @param pt2xx Coefficient 1 of x on point 2
     * @param pt2xy Coefficient 2 of x on point 2
     * @param pt2yx Coefficient 1 of y on point 2
     * @param pt2yy Coefficient 2 of y on point 2
     * @return (pt3xx, pt3xy, pt3yx, pt3yy)
     */
    function ECTwistAdd(
        uint256 pt1xx, uint256 pt1xy,
        uint256 pt1yx, uint256 pt1yy,
        uint256 pt2xx, uint256 pt2xy,
        uint256 pt2yx, uint256 pt2yy
    ) public view returns (
        uint256, uint256,
        uint256, uint256
    ) {
        if (
            pt1xx == 0 && pt1xy == 0 &&
            pt1yx == 0 && pt1yy == 0
        ) {
            if (!(
                pt2xx == 0 && pt2xy == 0 &&
                pt2yx == 0 && pt2yy == 0
            )) {
                assert(_isOnCurve(
                    pt2xx, pt2xy,
                    pt2yx, pt2yy
                ));
            }
            return (
                pt2xx, pt2xy,
                pt2yx, pt2yy
            );
        } else if (
            pt2xx == 0 && pt2xy == 0 &&
            pt2yx == 0 && pt2yy == 0
        ) {
            assert(_isOnCurve(
                pt1xx, pt1xy,
                pt1yx, pt1yy
            ));
            return (
                pt1xx, pt1xy,
                pt1yx, pt1yy
            );
        }

        assert(_isOnCurve(
            pt1xx, pt1xy,
            pt1yx, pt1yy
        ));
        assert(_isOnCurve(
            pt2xx, pt2xy,
            pt2yx, pt2yy
        ));

        uint256[6] memory pt3 = _ECTwistAddJacobian(
            pt1xx, pt1xy,
            pt1yx, pt1yy,
            1,     0,
            pt2xx, pt2xy,
            pt2yx, pt2yy,
            1,     0
        );

        return _fromJacobian(
            pt3[PTXX], pt3[PTXY],
            pt3[PTYX], pt3[PTYY],
            pt3[PTZX], pt3[PTZY]
        );
    }

    /**
     * @notice Multiply a twist point by a scalar
     * @param s     Scalar to multiply by
     * @param pt1xx Coefficient 1 of x
     * @param pt1xy Coefficient 2 of x
     * @param pt1yx Coefficient 1 of y
     * @param pt1yy Coefficient 2 of y
     * @return (pt2xx, pt2xy, pt2yx, pt2yy)
     */
    function ECTwistMul(
        uint256 s,
        uint256 pt1xx, uint256 pt1xy,
        uint256 pt1yx, uint256 pt1yy
    ) public view returns (
        uint256, uint256,
        uint256, uint256
    ) {
        uint256 pt1zx = 1;
        if (
            pt1xx == 0 && pt1xy == 0 &&
            pt1yx == 0 && pt1yy == 0
        ) {
            pt1xx = 1;
            pt1yx = 1;
            pt1zx = 0;
        } else {
            assert(_isOnCurve(
                pt1xx, pt1xy,
                pt1yx, pt1yy
            ));
        }

        uint256[6] memory pt2 = _ECTwistMulJacobian(
            s,
            pt1xx, pt1xy,
            pt1yx, pt1yy,
            pt1zx, 0
        );

        return _fromJacobian(
            pt2[PTXX], pt2[PTXY],
            pt2[PTYX], pt2[PTYY],
            pt2[PTZX], pt2[PTZY]
        );
    }

    /**
     * @notice Get the field modulus
     * @return The field modulus
     */
    function GetFieldModulus() public pure returns (uint256) {
        return FIELD_MODULUS;
    }

    function submod(uint256 a, uint256 b, uint256 n) internal pure returns (uint256) {
        return addmod(a, n - b, n);
    }

    function _FQ2Mul(
        uint256 xx, uint256 xy,
        uint256 yx, uint256 yy
    ) internal pure returns (uint256, uint256) {
        return (
            submod(mulmod(xx, yx, FIELD_MODULUS), mulmod(xy, yy, FIELD_MODULUS), FIELD_MODULUS),
            addmod(mulmod(xx, yy, FIELD_MODULUS), mulmod(xy, yx, FIELD_MODULUS), FIELD_MODULUS)
        );
    }

    function _FQ2Muc(
        uint256 xx, uint256 xy,
        uint256 c
    ) internal pure returns (uint256, uint256) {
        return (
            mulmod(xx, c, FIELD_MODULUS),
            mulmod(xy, c, FIELD_MODULUS)
        );
    }

    function _FQ2Add(
        uint256 xx, uint256 xy,
        uint256 yx, uint256 yy
    ) internal pure returns (uint256, uint256) {
        return (
            addmod(xx, yx, FIELD_MODULUS),
            addmod(xy, yy, FIELD_MODULUS)
        );
    }

    function _FQ2Sub(
        uint256 xx, uint256 xy,
        uint256 yx, uint256 yy
    ) internal pure returns (uint256 rx, uint256 ry) {
        return (
            submod(xx, yx, FIELD_MODULUS),
            submod(xy, yy, FIELD_MODULUS)
        );
    }

    function _FQ2Div(
        uint256 xx, uint256 xy,
        uint256 yx, uint256 yy
    ) internal view returns (uint256, uint256) {
        (yx, yy) = _FQ2Inv(yx, yy);
        return _FQ2Mul(xx, xy, yx, yy);
    }

    function _FQ2Inv(uint256 x, uint256 y) internal view returns (uint256, uint256) {
        uint256 inv = _modInv(addmod(mulmod(y, y, FIELD_MODULUS), mulmod(x, x, FIELD_MODULUS), FIELD_MODULUS), FIELD_MODULUS);
        return (
            mulmod(x, inv, FIELD_MODULUS),
            FIELD_MODULUS - mulmod(y, inv, FIELD_MODULUS)
        );
    }

    function _isOnCurve(
        uint256 xx, uint256 xy,
        uint256 yx, uint256 yy
    ) internal pure returns (bool) {
        uint256 yyx;
        uint256 yyy;
        uint256 xxxx;
        uint256 xxxy;
        (yyx, yyy) = _FQ2Mul(yx, yy, yx, yy);
        (xxxx, xxxy) = _FQ2Mul(xx, xy, xx, xy);
        (xxxx, xxxy) = _FQ2Mul(xxxx, xxxy, xx, xy);
        (yyx, yyy) = _FQ2Sub(yyx, yyy, xxxx, xxxy);
        (yyx, yyy) = _FQ2Sub(yyx, yyy, TWISTBX, TWISTBY);
        return yyx == 0 && yyy == 0;
    }

    function _modInv(uint256 a, uint256 n) internal view returns (uint256 result) {
        bool success;
        assembly {
            let freemem := mload(0x40)
            mstore(freemem, 0x20)
            mstore(add(freemem,0x20), 0x20)
            mstore(add(freemem,0x40), 0x20)
            mstore(add(freemem,0x60), a)
            mstore(add(freemem,0x80), sub(n, 2))
            mstore(add(freemem,0xA0), n)
            success := staticcall(sub(gas, 2000), 5, freemem, 0xC0, freemem, 0x20)
            result := mload(freemem)
        }
        require(success);
    }

    function _fromJacobian(
        uint256 pt1xx, uint256 pt1xy,
        uint256 pt1yx, uint256 pt1yy,
        uint256 pt1zx, uint256 pt1zy
    ) internal view returns (
        uint256 pt2xx, uint256 pt2xy,
        uint256 pt2yx, uint256 pt2yy
    ) {
        uint256 invzx;
        uint256 invzy;
        (invzx, invzy) = _FQ2Inv(pt1zx, pt1zy);
        (pt2xx, pt2xy) = _FQ2Mul(pt1xx, pt1xy, invzx, invzy);
        (pt2yx, pt2yy) = _FQ2Mul(pt1yx, pt1yy, invzx, invzy);
    }

    function _ECTwistAddJacobian(
        uint256 pt1xx, uint256 pt1xy,
        uint256 pt1yx, uint256 pt1yy,
        uint256 pt1zx, uint256 pt1zy,
        uint256 pt2xx, uint256 pt2xy,
        uint256 pt2yx, uint256 pt2yy,
        uint256 pt2zx, uint256 pt2zy) internal pure returns (uint256[6] memory pt3) {
            if (pt1zx == 0 && pt1zy == 0) {
                (
                    pt3[PTXX], pt3[PTXY],
                    pt3[PTYX], pt3[PTYY],
                    pt3[PTZX], pt3[PTZY]
                ) = (
                    pt2xx, pt2xy,
                    pt2yx, pt2yy,
                    pt2zx, pt2zy
                );
                return pt3;
            } else if (pt2zx == 0 && pt2zy == 0) {
                (
                    pt3[PTXX], pt3[PTXY],
                    pt3[PTYX], pt3[PTYY],
                    pt3[PTZX], pt3[PTZY]
                ) = (
                    pt1xx, pt1xy,
                    pt1yx, pt1yy,
                    pt1zx, pt1zy
                );
                return pt3;
            }

            (pt2yx,     pt2yy)     = _FQ2Mul(pt2yx, pt2yy, pt1zx, pt1zy); // U1 = y2 * z1
            (pt3[PTYX], pt3[PTYY]) = _FQ2Mul(pt1yx, pt1yy, pt2zx, pt2zy); // U2 = y1 * z2
            (pt2xx,     pt2xy)     = _FQ2Mul(pt2xx, pt2xy, pt1zx, pt1zy); // V1 = x2 * z1
            (pt3[PTZX], pt3[PTZY]) = _FQ2Mul(pt1xx, pt1xy, pt2zx, pt2zy); // V2 = x1 * z2

            if (pt2xx == pt3[PTZX] && pt2xy == pt3[PTZY]) {
                if (pt2yx == pt3[PTYX] && pt2yy == pt3[PTYY]) {
                    (
                        pt3[PTXX], pt3[PTXY],
                        pt3[PTYX], pt3[PTYY],
                        pt3[PTZX], pt3[PTZY]
                    ) = _ECTwistDoubleJacobian(pt1xx, pt1xy, pt1yx, pt1yy, pt1zx, pt1zy);
                    return pt3;
                }
                (
                    pt3[PTXX], pt3[PTXY],
                    pt3[PTYX], pt3[PTYY],
                    pt3[PTZX], pt3[PTZY]
                ) = (
                    1, 0,
                    1, 0,
                    0, 0
                );
                return pt3;
            }

            (pt2zx,     pt2zy)     = _FQ2Mul(pt1zx, pt1zy, pt2zx,     pt2zy);     // W = z1 * z2
            (pt1xx,     pt1xy)     = _FQ2Sub(pt2yx, pt2yy, pt3[PTYX], pt3[PTYY]); // U = U1 - U2
            (pt1yx,     pt1yy)     = _FQ2Sub(pt2xx, pt2xy, pt3[PTZX], pt3[PTZY]); // V = V1 - V2
            (pt1zx,     pt1zy)     = _FQ2Mul(pt1yx, pt1yy, pt1yx,     pt1yy);     // V_squared = V * V
            (pt2yx,     pt2yy)     = _FQ2Mul(pt1zx, pt1zy, pt3[PTZX], pt3[PTZY]); // V_squared_times_V2 = V_squared * V2
            (pt1zx,     pt1zy)     = _FQ2Mul(pt1zx, pt1zy, pt1yx,     pt1yy);     // V_cubed = V * V_squared
            (pt3[PTZX], pt3[PTZY]) = _FQ2Mul(pt1zx, pt1zy, pt2zx,     pt2zy);     // newz = V_cubed * W
            (pt2xx,     pt2xy)     = _FQ2Mul(pt1xx, pt1xy, pt1xx,     pt1xy);     // U * U
            (pt2xx,     pt2xy)     = _FQ2Mul(pt2xx, pt2xy, pt2zx,     pt2zy);     // U * U * W
            (pt2xx,     pt2xy)     = _FQ2Sub(pt2xx, pt2xy, pt1zx,     pt1zy);     // U * U * W - V_cubed
            (pt2zx,     pt2zy)     = _FQ2Muc(pt2yx, pt2yy, 2);                    // 2 * V_squared_times_V2
            (pt2xx,     pt2xy)     = _FQ2Sub(pt2xx, pt2xy, pt2zx,     pt2zy);     // A = U * U * W - V_cubed - 2 * V_squared_times_V2
            (pt3[PTXX], pt3[PTXY]) = _FQ2Mul(pt1yx, pt1yy, pt2xx,     pt2xy);     // newx = V * A
            (pt1yx,     pt1yy)     = _FQ2Sub(pt2yx, pt2yy, pt2xx,     pt2xy);     // V_squared_times_V2 - A
            (pt1yx,     pt1yy)     = _FQ2Mul(pt1xx, pt1xy, pt1yx,     pt1yy);     // U * (V_squared_times_V2 - A)
            (pt1xx,     pt1xy)     = _FQ2Mul(pt1zx, pt1zy, pt3[PTYX], pt3[PTYY]); // V_cubed * U2
            (pt3[PTYX], pt3[PTYY]) = _FQ2Sub(pt1yx, pt1yy, pt1xx,     pt1xy);     // newy = U * (V_squared_times_V2 - A) - V_cubed * U2
    }

    function _ECTwistDoubleJacobian(
        uint256 pt1xx, uint256 pt1xy,
        uint256 pt1yx, uint256 pt1yy,
        uint256 pt1zx, uint256 pt1zy
    ) internal pure returns (
        uint256 pt2xx, uint256 pt2xy,
        uint256 pt2yx, uint256 pt2yy,
        uint256 pt2zx, uint256 pt2zy
    ) {
        (pt2xx, pt2xy) = _FQ2Muc(pt1xx, pt1xy, 3);            // 3 * x
        (pt2xx, pt2xy) = _FQ2Mul(pt2xx, pt2xy, pt1xx, pt1xy); // W = 3 * x * x
        (pt1zx, pt1zy) = _FQ2Mul(pt1yx, pt1yy, pt1zx, pt1zy); // S = y * z
        (pt2yx, pt2yy) = _FQ2Mul(pt1xx, pt1xy, pt1yx, pt1yy); // x * y
        (pt2yx, pt2yy) = _FQ2Mul(pt2yx, pt2yy, pt1zx, pt1zy); // B = x * y * S
        (pt1xx, pt1xy) = _FQ2Mul(pt2xx, pt2xy, pt2xx, pt2xy); // W * W
        (pt2zx, pt2zy) = _FQ2Muc(pt2yx, pt2yy, 8);            // 8 * B
        (pt1xx, pt1xy) = _FQ2Sub(pt1xx, pt1xy, pt2zx, pt2zy); // H = W * W - 8 * B
        (pt2zx, pt2zy) = _FQ2Mul(pt1zx, pt1zy, pt1zx, pt1zy); // S_squared = S * S
        (pt2yx, pt2yy) = _FQ2Muc(pt2yx, pt2yy, 4);            // 4 * B
        (pt2yx, pt2yy) = _FQ2Sub(pt2yx, pt2yy, pt1xx, pt1xy); // 4 * B - H
        (pt2yx, pt2yy) = _FQ2Mul(pt2yx, pt2yy, pt2xx, pt2xy); // W * (4 * B - H)
        (pt2xx, pt2xy) = _FQ2Muc(pt1yx, pt1yy, 8);            // 8 * y
        (pt2xx, pt2xy) = _FQ2Mul(pt2xx, pt2xy, pt1yx, pt1yy); // 8 * y * y
        (pt2xx, pt2xy) = _FQ2Mul(pt2xx, pt2xy, pt2zx, pt2zy); // 8 * y * y * S_squared
        (pt2yx, pt2yy) = _FQ2Sub(pt2yx, pt2yy, pt2xx, pt2xy); // newy = W * (4 * B - H) - 8 * y * y * S_squared
        (pt2xx, pt2xy) = _FQ2Muc(pt1xx, pt1xy, 2);            // 2 * H
        (pt2xx, pt2xy) = _FQ2Mul(pt2xx, pt2xy, pt1zx, pt1zy); // newx = 2 * H * S
        (pt2zx, pt2zy) = _FQ2Mul(pt1zx, pt1zy, pt2zx, pt2zy); // S * S_squared
        (pt2zx, pt2zy) = _FQ2Muc(pt2zx, pt2zy, 8);            // newz = 8 * S * S_squared
    }

    function _ECTwistMulJacobian(
        uint256 d,
        uint256 pt1xx, uint256 pt1xy,
        uint256 pt1yx, uint256 pt1yy,
        uint256 pt1zx, uint256 pt1zy
    ) internal pure returns (uint256[6] memory pt2) {
        while (d != 0) {
            if ((d & 1) != 0) {
                pt2 = _ECTwistAddJacobian(
                    pt2[PTXX], pt2[PTXY],
                    pt2[PTYX], pt2[PTYY],
                    pt2[PTZX], pt2[PTZY],
                    pt1xx, pt1xy,
                    pt1yx, pt1yy,
                    pt1zx, pt1zy);
            }
            (
                pt1xx, pt1xy,
                pt1yx, pt1yy,
                pt1zx, pt1zy
            ) = _ECTwistDoubleJacobian(
                pt1xx, pt1xy,
                pt1yx, pt1yy,
                pt1zx, pt1zy
            );

            d = d / 2;
        }
    }
}
// This file is MIT Licensed.
//
// Copyright 2017 Christian Reitwiessner
// Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
// The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
pragma solidity ^0.5.0;
library Pairing {
    struct G1Point {
        uint X;
        uint Y;
    }
    // Encoding of field elements is: X[0] * z + X[1]
    struct G2Point {
        uint[2] X;
        uint[2] Y;
    }
    /// @return the generator of G1
    function P1() pure internal returns (G1Point memory) {
        return G1Point(1, 2);
    }
    /// @return the generator of G2
    function P2() pure internal returns (G2Point memory) {
        return G2Point(
            [11559732032986387107991004021392285783925812861821192530917403151452391805634,
             10857046999023057135944570762232829481370756359578518086990519993285655852781],
            [4082367875863433681332203403145435568316851327593401208105741076214120093531,
             8495653923123431417604973247489272438418190587263600148770280649306958101930]
        );
    }
    /// @return the negation of p, i.e. p.addition(p.negate()) should be zero.
    function negate(G1Point memory p) pure internal returns (G1Point memory) {
        // The prime q in the base field F_q for G1
        uint q = 21888242871839275222246405745257275088696311157297823662689037894645226208583;
        if (p.X == 0 && p.Y == 0)
            return G1Point(0, 0);
        return G1Point(p.X, q - (p.Y % q));
    }
    /// @return the sum of two points of G1
    function addition(G1Point memory p1, G1Point memory p2) internal returns (G1Point memory r) {
        uint[4] memory input;
        input[0] = p1.X;
        input[1] = p1.Y;
        input[2] = p2.X;
        input[3] = p2.Y;
        bool success;
        assembly {
            success := call(sub(gas, 2000), 6, 0, input, 0xc0, r, 0x60)
            // Use "invalid" to make gas estimation work
            switch success case 0 { invalid() }
        }
        require(success);
    }
    /// @return the sum of two points of G2
    function addition(G2Point memory p1, G2Point memory p2) internal returns (G2Point memory r) {
        (r.X[1], r.X[0], r.Y[1], r.Y[0]) = BN256G2.ECTwistAdd(p1.X[1],p1.X[0],p1.Y[1],p1.Y[0],p2.X[1],p2.X[0],p2.Y[1],p2.Y[0]);
    }
    /// @return the product of a point on G1 and a scalar, i.e.
    /// p == p.scalar_mul(1) and p.addition(p) == p.scalar_mul(2) for all points p.
    function scalar_mul(G1Point memory p, uint s) internal returns (G1Point memory r) {
        uint[3] memory input;
        input[0] = p.X;
        input[1] = p.Y;
        input[2] = s;
        bool success;
        assembly {
            success := call(sub(gas, 2000), 7, 0, input, 0x80, r, 0x60)
            // Use "invalid" to make gas estimation work
            switch success case 0 { invalid() }
        }
        require (success);
    }
    /// @return the result of computing the pairing check
    /// e(p1[0], p2[0]) *  .... * e(p1[n], p2[n]) == 1
    /// For example pairing([P1(), P1().negate()], [P2(), P2()]) should
    /// return true.
    function pairing(G1Point[] memory p1, G2Point[] memory p2) internal returns (bool) {
        require(p1.length == p2.length);
        uint elements = p1.length;
        uint inputSize = elements * 6;
        uint[] memory input = new uint[](inputSize);
        for (uint i = 0; i < elements; i++)
        {
            input[i * 6 + 0] = p1[i].X;
            input[i * 6 + 1] = p1[i].Y;
            input[i * 6 + 2] = p2[i].X[0];
            input[i * 6 + 3] = p2[i].X[1];
            input[i * 6 + 4] = p2[i].Y[0];
            input[i * 6 + 5] = p2[i].Y[1];
        }
        uint[1] memory out;
        bool success;
        assembly {
            success := call(sub(gas, 2000), 8, 0, add(input, 0x20), mul(inputSize, 0x20), out, 0x20)
            // Use "invalid" to make gas estimation work
            switch success case 0 { invalid() }
        }
        require(success);
        return out[0] != 0;
    }
    /// Convenience method for a pairing check for two pairs.
    function pairingProd2(G1Point memory a1, G2Point memory a2, G1Point memory b1, G2Point memory b2) internal returns (bool) {
        G1Point[] memory p1 = new G1Point[](2);
        G2Point[] memory p2 = new G2Point[](2);
        p1[0] = a1;
        p1[1] = b1;
        p2[0] = a2;
        p2[1] = b2;
        return pairing(p1, p2);
    }
    /// Convenience method for a pairing check for three pairs.
    function pairingProd3(
            G1Point memory a1, G2Point memory a2,
            G1Point memory b1, G2Point memory b2,
            G1Point memory c1, G2Point memory c2
    ) internal returns (bool) {
        G1Point[] memory p1 = new G1Point[](3);
        G2Point[] memory p2 = new G2Point[](3);
        p1[0] = a1;
        p1[1] = b1;
        p1[2] = c1;
        p2[0] = a2;
        p2[1] = b2;
        p2[2] = c2;
        return pairing(p1, p2);
    }
    /// Convenience method for a pairing check for four pairs.
    function pairingProd4(
            G1Point memory a1, G2Point memory a2,
            G1Point memory b1, G2Point memory b2,
            G1Point memory c1, G2Point memory c2,
            G1Point memory d1, G2Point memory d2
    ) internal returns (bool) {
        G1Point[] memory p1 = new G1Point[](4);
        G2Point[] memory p2 = new G2Point[](4);
        p1[0] = a1;
        p1[1] = b1;
        p1[2] = c1;
        p1[3] = d1;
        p2[0] = a2;
        p2[1] = b2;
        p2[2] = c2;
        p2[3] = d2;
        return pairing(p1, p2);
    }
}

contract Verifier {
    using Pairing for *;
    struct VerifyingKey {
        Pairing.G1Point a;
        Pairing.G2Point b;
        Pairing.G2Point gamma;
        Pairing.G2Point delta;
        Pairing.G1Point[] gamma_abc;
    }
    struct Proof {
        Pairing.G1Point a;
        Pairing.G2Point b;
        Pairing.G1Point c;
    }
    function verifyingKey() pure internal returns (VerifyingKey memory vk) {
        vk.a = Pairing.G1Point(uint256(0x00b98e074c4132dd355d10981a4c8a76fc18c7a3227a3a2f04944c6d089ac65f), uint256(0x07d3a00bdd364cb3f14ba2aae3193fcea361b8d94fd1354b0374407ee91067be));
        vk.b = Pairing.G2Point([uint256(0x0467ef878336029dcb1841d5ee7719397c3ca7d413f6b07bc567d55803aa344e), uint256(0x211efa9a2e8fb503781b7877a8153d99e6b0a68a9ce8b8f60f0c546462ebefa9)], [uint256(0x2abcf347dbd734e052b45c249b363426a7aa4c41eb20dd4c7c2bbd91974db1d0), uint256(0x0261c2bbab220a80e28a983b7f2f3e485afed3f4bc4a56e1e360bef2e301a8ae)]);
        vk.gamma = Pairing.G2Point([uint256(0x00f17431de8ff2fa172cb5112b3759b529fa76b50901589aa8a3698aeaf005cd), uint256(0x026bf748346871c350c2b0eb7d9e269f6731bcaed19c9c26801be0a2b155460e)], [uint256(0x07e67e27bc436f8506a58da2529316e73b7a19b4a4bfce508f7b796ddfd4b822), uint256(0x0bd422ff085906d25b3f6aa362d2223b3448aa37ad9d13eacb1692b8a8212ff7)]);
        vk.delta = Pairing.G2Point([uint256(0x0a350ee3a469d089acab1cbd42b16643b7b6fa1ce6f0ba46809bf42e9d676898), uint256(0x004e29ffc5a7595b97a49e4b6284662467afca687add1a6fe54086bf601a0418)], [uint256(0x0ecb1950e3aac67fe8c17f7d92881d56a70d5704ea36fd5ae20e4a08d400ebc5), uint256(0x1ffeaa82a983411ff6e3ac43441204604d71f4d582815ba92af34f810cf217a2)]);
        vk.gamma_abc = new Pairing.G1Point[](39);
        vk.gamma_abc[0] = Pairing.G1Point(uint256(0x15789245093097411c8e29d70fd36bbaf444845283ddbdcb4ea6c9bfaebb222d), uint256(0x10b1115e5ea4140bb8f0eb58f1e695e3bb6a4eebeaa304dc10d12a79b4e86843));
        vk.gamma_abc[1] = Pairing.G1Point(uint256(0x1cc382a10ca7074f5b20031efaedb5f45e69cadd64b8e88f33bd89e6778d7001), uint256(0x0e8741bf40f625d2960fb5b7ff8d2e2b86298a44cc4002da74de9bb38a5ea7a9));
        vk.gamma_abc[2] = Pairing.G1Point(uint256(0x2e460d84a78bb74c7318bef4a88b0854d07ac780e3c23e8d4aa47f1809879fbe), uint256(0x1a9109baebc64edbff114e3119db092508918d0b40ad47b8c15bffc363751e86));
        vk.gamma_abc[3] = Pairing.G1Point(uint256(0x1e98fb09312faefdcda91fedb634e37b884a9e096e6042cc051b6c424f8cdf62), uint256(0x2f785f84908ea37976f9e47d9ebbb65336a85f8f13b779e6e7a8f40f970a5c3a));
        vk.gamma_abc[4] = Pairing.G1Point(uint256(0x20ed814892eae12c58a59f5b6f068b790d95e19c78d286595d5724b7967b6860), uint256(0x1f7ea7874a89b5ba7b333ac20ee0e0b29dce498821213c2f621bf2f0d4a5b048));
        vk.gamma_abc[5] = Pairing.G1Point(uint256(0x22928dd17c2cb9170915830034d6352acc51c26e2933f1580a8ee02f5e9adb77), uint256(0x12dfb1438f84004cb0a014cd7ec7658395187a70fd74c100ae8e39cdb142bd91));
        vk.gamma_abc[6] = Pairing.G1Point(uint256(0x088601aff3ae2fe53510df7f7a9efb2dccd51047b84344d19d4df8df83a4aaa9), uint256(0x176c8a55c59437950aee567fef64e1caa67665a38c215ba9b45eda6d7e4dc554));
        vk.gamma_abc[7] = Pairing.G1Point(uint256(0x0d17b79c216233c9ab9df8c9464cd38430c06654042472c9d3b32667db4251db), uint256(0x25a136a2c46424d2feae91853847a8fce35c6648e86dfece5ceb2ee6a5cdffb5));
        vk.gamma_abc[8] = Pairing.G1Point(uint256(0x0485395f18640fadedac3cfaa6cc1d746f5f2c39cf63b2f7118c45434c31ff85), uint256(0x2522a84370254b1fd208de89397b091eb4256bd4076b2e5aed23bf34d6662b27));
        vk.gamma_abc[9] = Pairing.G1Point(uint256(0x2cd05a915a14ca4efc1d071350f94e27c58414f43e7709dc67cc79734606887e), uint256(0x07764bdf266b65e1a106ad8a2a88967fb33e3fa51a3d6d43b115ec8dde2674a1));
        vk.gamma_abc[10] = Pairing.G1Point(uint256(0x019edf500a6225fd77aac7fedd817370ba6b6a752c2f994587adfbb1c6b83c4d), uint256(0x165c1f9eb96886c487b133bf43cae38aac58d662ef1b63e165a5dc3e6ba5de0c));
        vk.gamma_abc[11] = Pairing.G1Point(uint256(0x170330ba0ab81eec241743e40f99ebe2fc7e99df61731a71d681d5ce015616de), uint256(0x2c7b38d88fb257fedad1407ff94d73f6aef9719891969da2e7904e660812ca5b));
        vk.gamma_abc[12] = Pairing.G1Point(uint256(0x0b9ff51048bb60d23aa19ecd09e3089e1205d5c3d532f2945b43c99756d8a8f8), uint256(0x23c16e3bf553314d7541c4e827878689a0256e909178c834e417d6c9ee2d3b68));
        vk.gamma_abc[13] = Pairing.G1Point(uint256(0x209c33384f58c94aafcf9bd0d2236d120dda7f8e1f8a34e78a1973e61168348b), uint256(0x2c01f1f701bc7d4e1dd516f0d921d440c2e9e8b3ae892e51c087e0c44451a57c));
        vk.gamma_abc[14] = Pairing.G1Point(uint256(0x11dd1fe1295de392c0d6e93bce20faec2794807fcb4e2ab7d913d25eab7afe00), uint256(0x0164b910acc76cd265604ba1b91045f16d5a8ce574d7e12908db20d25063117a));
        vk.gamma_abc[15] = Pairing.G1Point(uint256(0x20f05dde7c0c1d64f7ef6fd4ccce405ba595ba4351845af4e3b1401218e6eaf3), uint256(0x2e6e962312cb37603eb21b3a552ad90fd31abfdab3a2bf18254dead7c1166191));
        vk.gamma_abc[16] = Pairing.G1Point(uint256(0x28a0177fd2c01063ccb0e73a2bc387e28f380108011a50ec426bd3a973b828e9), uint256(0x0b97f1c0240be2e76840018453962d1384820e4b9617b0a3478732eb8d8c7076));
        vk.gamma_abc[17] = Pairing.G1Point(uint256(0x0775882bfa4b47b35d9041d8567bd57bef8c2b1286d7d68ecf894c3e236b43b2), uint256(0x10c529ef1d1e55e966815616e097c1e68a2ae9237303651533ea05c7b836e735));
        vk.gamma_abc[18] = Pairing.G1Point(uint256(0x2d34357be30baa46966d01e201190f376c4220ee6a2f96511ec8d67f6e99f02d), uint256(0x1f74b388629060de3797f3284770dbc41f79cc47a43450b5dbbcd1546b683a4a));
        vk.gamma_abc[19] = Pairing.G1Point(uint256(0x2dfe561a15196b0cdbd9ce05fff675ebba801c50c29e6ed009cbb9be62b44a12), uint256(0x2b1acfb1b34b465028a5f8e9b2672572ecd2c3690871c6e93bd202144a8d25d0));
        vk.gamma_abc[20] = Pairing.G1Point(uint256(0x034ea4ff19a965ee5eca40dbcc8f8a25bdf34bf8203377e0a1d828a5db4e564d), uint256(0x0ddc61873273dee5410ecfaacc99445383d9633b3bae96e3761c5f40a20595d3));
        vk.gamma_abc[21] = Pairing.G1Point(uint256(0x07510808e69b9676cfb4e6eb150c3110d1362603bcc06a958ae0cad993110471), uint256(0x22b9c3ae68ba9a01f4012c91ac54f058e51fa39edee21b5505e7057016c2adaf));
        vk.gamma_abc[22] = Pairing.G1Point(uint256(0x1d5a0854b15399058a304404d07186398700eccacfd5e4f5f2e3ffa4dc68926d), uint256(0x10fd44ecd539f2e358b4b15527e7b26d26cbc339949917e5d336982476a84d1a));
        vk.gamma_abc[23] = Pairing.G1Point(uint256(0x078c9a647b8f4be86fadd4de07b2562f7d28d86fb9a1feed6f4a1b69a1db1deb), uint256(0x2abe5c86fec864870768494e2c38f820cf241f7f9d49b9015c7ebdb685bf7327));
        vk.gamma_abc[24] = Pairing.G1Point(uint256(0x0d5da8edb4f59eecfaa854c7420b798efdbdd6da10fb4edf9e3508cc36828001), uint256(0x0b291d7ecf34ca6d6574e08b5e47c1b5582e37a983b7f4a206d6a7271a63b82d));
        vk.gamma_abc[25] = Pairing.G1Point(uint256(0x2576c1aade932071c660e4d011d4d5df041a25205f7e83eee8afedc13bef100b), uint256(0x14c8ccf6eb1325207702104663dab6bea7ab5eec425353089cebc9faed1d9147));
        vk.gamma_abc[26] = Pairing.G1Point(uint256(0x27dd615a1ecd421ef81312ac12190f3dd1074e4fea47156f077b52467faf10d3), uint256(0x0ff0e31db9ba5ecd5912a69fb5b08096de7579f46c50cfe40c59a7539e36a4f2));
        vk.gamma_abc[27] = Pairing.G1Point(uint256(0x2c926364ae3f6a56c36b04903949e2e5a8bc77838e8bdfa2f89e632535cdcf94), uint256(0x12341a5dd75a3c604ca0eccb5fd8099c9df7818e8b7b39afda320fc152d28d73));
        vk.gamma_abc[28] = Pairing.G1Point(uint256(0x2a0f890407c909591c9296d0e4d5e4ee10a76e42f4a0ca53ec05628e84bc2be1), uint256(0x0881a2c90eb9a116143c03536461f09591fa0b499ea88ce305ad21992fe3870b));
        vk.gamma_abc[29] = Pairing.G1Point(uint256(0x2b0c4fcefdb40d0fcad693354a9f16f11ad4801e3d1e6ce6872e3c54bc5a5bfb), uint256(0x252b031435e006967d2cc62ee68b7223411ca0f359279af4a8d1498e5e842a52));
        vk.gamma_abc[30] = Pairing.G1Point(uint256(0x2c1a22207520e690d037bca6ad76a91d66e282e9da228fc0e40e1def52ffb2ab), uint256(0x1886485a62617e9327ffab70955d8886ff04319345f475c4af144654ba561430));
        vk.gamma_abc[31] = Pairing.G1Point(uint256(0x2657f10f8314fc3cbf09f3ad0bc8ff29bae412b2fcb1406b5b7fa274c65e8ed0), uint256(0x0c04959a6434e40b9d4cbe3932d953fe618382e424cf0b0861f3e74430546e4f));
        vk.gamma_abc[32] = Pairing.G1Point(uint256(0x28c47a069a8d0f256d927cdb601d7efbbe4fac6ca0ced0ff110156be9bea2480), uint256(0x0fba7b6de8da72cb7d99041c544f6b22e8a0edca11d6cd70318b80921d142582));
        vk.gamma_abc[33] = Pairing.G1Point(uint256(0x10220ab3ba94e33b7a3f22454e78fce88eed774132748c695be2c9d6544a10a5), uint256(0x0b4ae895c4bfbf81b75c35a71a6b5081a5c2c7ca043cfcada82d1c9dcf72edcb));
        vk.gamma_abc[34] = Pairing.G1Point(uint256(0x2b47e397128c51e97a5ec4db7f477bca209d6cf87d5085e8940aaa67faa5c44d), uint256(0x1e77014c89fa6794f69e6d97ea5a7e0994f2516c017f2a53b7ef3982ae6ce64d));
        vk.gamma_abc[35] = Pairing.G1Point(uint256(0x272946117e4f800b661dffd7ac3ab8887105f2ee74b007790916f4935428893c), uint256(0x15cb0509b62df58fbee6c63d85735dfab99f29df88671a3ebcf83b8c22be61fb));
        vk.gamma_abc[36] = Pairing.G1Point(uint256(0x2c7ea77e03cbe23edd521cf2a45e77c989331b4a87ddb08438c907f862f53690), uint256(0x1bffbec58daa15ab946352bb1052a6103885a5b201a36143504544c7dcd3e466));
        vk.gamma_abc[37] = Pairing.G1Point(uint256(0x2dde096dd25afbcdb370acd3c9b57ed0c50759c83558708ddf385324b951698c), uint256(0x024bef449cab992635dea5a4ff70368e059d91aadef2da6dc2140e48430575af));
        vk.gamma_abc[38] = Pairing.G1Point(uint256(0x25a52a09d8a3636d8bc103ee02750c4fbfc349700287033169e86f1ddd921217), uint256(0x04923d1da09d3cdae2c1817d9dd92b56af34395a9adcfed16f2e6ee95bf39fc4));
    }
    function verify(uint[] memory input, Proof memory proof) internal returns (uint) {
        uint256 snark_scalar_field = 21888242871839275222246405745257275088548364400416034343698204186575808495617;
        VerifyingKey memory vk = verifyingKey();
        require(input.length + 1 == vk.gamma_abc.length);
        // Compute the linear combination vk_x
        Pairing.G1Point memory vk_x = Pairing.G1Point(0, 0);
        for (uint i = 0; i < input.length; i++) {
            require(input[i] < snark_scalar_field);
            vk_x = Pairing.addition(vk_x, Pairing.scalar_mul(vk.gamma_abc[i + 1], input[i]));
        }
        vk_x = Pairing.addition(vk_x, vk.gamma_abc[0]);
        if(!Pairing.pairingProd4(
             proof.a, proof.b,
             Pairing.negate(vk_x), vk.gamma,
             Pairing.negate(proof.c), vk.delta,
             Pairing.negate(vk.a), vk.b)) return 1;
        return 0;
    }
    event Verified(string s);
    function verifyTx(
            uint[2] memory a,
            uint[2][2] memory b,
            uint[2] memory c,
            uint[38] memory input
        ) public returns (bool r) {
        Proof memory proof;
        proof.a = Pairing.G1Point(a[0], a[1]);
        proof.b = Pairing.G2Point([b[0][0], b[0][1]], [b[1][0], b[1][1]]);
        proof.c = Pairing.G1Point(c[0], c[1]);
        uint[] memory inputValues = new uint[](input.length);
        for(uint i = 0; i < input.length; i++){
            inputValues[i] = input[i];
        }
        if (verify(inputValues, proof) == 0) {
            emit Verified("Transaction successfully verified.");
            return true;
        } else {
            return false;
        }
    }
}
