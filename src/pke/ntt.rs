use super::{N, ZETA};
use crate::pke::field::{barrett_reduce, mont_reduce};
use crate::pke::polynomial::Poly;

const INV_NTT_REDUCTIONS: [i16; 79] = [
    -1, // after layer 1
    -1, // after layer 2
    16, 17, 48, 49, 80, 81, 112, 113, 144, 145, 176, 177, 208, 209, 240, 241,
    -1, // after layer 3
    0, 1, 32, 33, 34, 35, 64, 65, 96, 97, 98, 99, 128, 129, 160, 161, 162, 163, 192, 193, 224, 225,
    226, 227, -1, // after layer 4
    2, 3, 66, 67, 68, 69, 70, 71, 130, 131, 194, 195, 196, 197, 198, 199, -1, // after layer 5
    4, 5, 6, 7, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143,
    -1, // after layer 6
    -1, //  after layer 7
];

pub fn ntt_generic(p: &mut Poly) {
    // Note that ℤ_q does not have a primitive 512ᵗʰ root of unity (as 512
    // does not divide into q-1) and so we cannot do a regular NTT.  ℤ_q
    // does have a primitive 256ᵗʰ root of unity, the smallest of which
    // is ζ := 17.
    //
    // Recall that our base ring R := ℤ_q[x] / (x²⁵⁶ + 1).  The polynomial
    // x²⁵⁶+1 will not split completely (as its roots would be 512ᵗʰ roots
    // of unity.)  However, it does split almost (using ζ¹²⁸ = -1):
    //
    // x²⁵⁶ + 1 = (x²)¹²⁸ - ζ¹²⁸
    //          = ((x²)⁶⁴ - ζ⁶⁴)((x²)⁶⁴ + ζ⁶⁴)
    //          = ((x²)³² - ζ³²)((x²)³² + ζ³²)((x²)³² - ζ⁹⁶)((x²)³² + ζ⁹⁶)
    //          ⋮
    //          = (x² - ζ)(x² + ζ)(x² - ζ⁶⁵)(x² + ζ⁶⁵) … (x² + ζ¹²⁷)
    //
    // Note that the powers of ζ that appear (from the second line down) are
    // in binary
    //
    // 0100000 1100000
    // 0010000 1010000 0110000 1110000
    // 0001000 1001000 0101000 1101000 0011000 1011000 0111000 1111000
    //         …
    //
    // That is: brv(2), brv(3), brv(4), …, where brv(x) denotes the 7-bit
    // bit-reversal of x.  These powers of ζ are given by the Zetas array.
    //
    // The polynomials x² ± ζⁱ are irreducible and coprime, hence by
    // the Chinese Remainder Theorem we know
    //
    //  ℤ_q[x]/(x²⁵⁶+1) → ℤ_q[x]/(x²-ζ) x … x  ℤ_q[x]/(x²+ζ¹²⁷)
    //
    // given by a ↦ ( a mod x²-ζ, …, a mod x²+ζ¹²⁷ )
    // is an isomorphism, which is the "NTT".  It can be efficiently computed by
    //
    //
    //  a ↦ ( a mod (x²)⁶⁴ - ζ⁶⁴, a mod (x²)⁶⁴ + ζ⁶⁴ )
    //    ↦ ( a mod (x²)³² - ζ³², a mod (x²)³² + ζ³²,
    //        a mod (x²)⁹⁶ - ζ⁹⁶, a mod (x²)⁹⁶ + ζ⁹⁶ )
    //
    //	    etc.
    //
    // If N was 8 then this can be pictured in the following diagram:
    //
    //  https://cnx.org/resources/17ee4dfe517a6adda05377b25a00bf6e6c93c334/File0026.png
    //
    // Each cross is a Cooley-Tukey butterfly: it's the map
    //
    //  (a, b) ↦ (a + ζb, a - ζb)
    //
    // for the appropriate power ζ for that column and row group.

    let mut k = 0;

    let mut l = N / 2;
    while l > 1 {
        let mut offset = 0;
        while offset < N - l {
            k += 1;
            let zeta = ZETA[k] as i32;
            let mut j = offset;
            while j < offset + l {
                let t = mont_reduce(zeta + p.coeffs[j + l] as i32);
                p.coeffs[j + l] = p.coeffs[j] - t;
                p.coeffs[j] = p.coeffs[j] + t;
                j += 1;
            }
            offset += 2 * l;
        }
        l >>= 1;
    }
}

pub fn inv_ntt_generic(p: &mut Poly) {
    let mut k = 127;
    let mut r = -1;

    let mut l = 2;
    while l < N {
        let mut offset = 0;
        while offset < N - l {
            k -= 1;
            let min_zeta = ZETA[k] as i32;
            let mut j = offset;
            while j < offset + l {
                let mut t = p.coeffs[j + l] - p.coeffs[j];
                p.coeffs[j] += p.coeffs[j + l];
                p.coeffs[j + l] = mont_reduce(min_zeta * t as i32);
                j += 1;
            }
            offset += 2 * l;
        }
        loop {
            r += 1;
            let i = INV_NTT_REDUCTIONS[r];
            if i < 0 {
                break;
            }
            p.coeffs[i] = barrett_reduce(p.coeffs[i]);
        }

        l <<= 1;
    }

    let mut j = 0;
    while j < N {
        // Note 1441 = (128)⁻¹ R².  The coefficients are bounded by 9q, so
        // as 1441 * 9 ≈ 2¹⁴ < 2¹⁵, we're within the required bounds
        // for montReduce().
        let red_factor = 1441;
        p.coeffs[j] = mont_reduce(red_factor * p.coeffs[j] as i32);
        j += 1;
    }
}
