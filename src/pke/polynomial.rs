use super::{N, Q};
use crate::pke::field::{barrett_reduce, csubq, mont_reduce, to_mont};
use crate::pke::ntt::ZETA;

type Byte = u8;
#[derive(Debug, Clone, Copy)]
pub struct Poly {
    pub coeffs: [i16; N],
}

impl Poly {
    pub fn new() -> Poly {
        Poly { coeffs: [0; N] }
    }

    fn add(&mut self, other: &Poly) {
        for i in 0..N {
            self.coeffs[i] = self.coeffs[i] + other.coeffs[i];
        }
    }

    fn sub(&mut self, other: &Poly) {
        for i in 0..N {
            self.coeffs[i] = self.coeffs[i] - other.coeffs[i];
        }
    }

    fn barrett_reduce(&mut self) {
        for i in 0..N {
            self.coeffs[i] = barrett_reduce(self.coeffs[i]);
        }
    }

    fn normalize(&mut self) {
        for i in 0..N {
            self.coeffs[i] = csubq(barrett_reduce(self.coeffs[i]));
        }
    }

    pub fn mont(&mut self) {
        for i in 0..N {
            self.coeffs[i] = to_mont(self.coeffs[i]);
        }
    }

    fn mul_hat(&mut self, other: &Poly) {
        let mut k: i32 = 64;

        let mut i = 0;

        while i < N {
            k += 1;
            let zeta = ZETA[k] as i32;

            let mut p0 = mont_reduce(self.coeffs[i + 1] as i32 * other.coeffs[i + 1] as i32);
            p0 = mont_reduce(zeta * p0 as i32);
            p0 += mont_reduce(self.coeffs[i] as i32 * other.coeffs[i] as i32);

            let mut p1 = mont_reduce(self.coeffs[i] as i32 * other.coeffs[i + 1] as i32);
            p1 += mont_reduce(self.coeffs[i + 1] as i32 * other.coeffs[i] as i32);

            self.coeffs[i] = p0;
            self.coeffs[i + 1] = p1;

            let mut p2 = mont_reduce(self.coeffs[i + 3] as i32 * other.coeffs[i + 3] as i32);
            p2 = -mont_reduce(zeta * p2 as i32);
            p2 += mont_reduce(self.coeffs[i + 2] as i32 * other.coeffs[i + 2] as i32);

            let mut p3 = mont_reduce(self.coeffs[i + 2] as i32 * other.coeffs[i + 3] as i32);
            p3 += mont_reduce(self.coeffs[i + 3] as i32 * other.coeffs[i + 2] as i32);

            self.coeffs[i + 2] = p2;
            self.coeffs[i + 3] = p3;

            i += 4;
        }
    }

    pub fn pack(&self, buf: &mut [Byte]) {
        // detangle(self); TODO: optimized amd64 version, till now we can assume std order
        let mut i = 0;
        while i < 128 {
            i += 1;
            let t0 = self.coeffs[2 * i] as u16;
            let t1 = self.coeffs[2 * i + 1] as u16;
            buf[3 * i] = t0 as Byte;
            buf[3 * i + 1] = (t0 >> 8) as Byte | (t1 << 4) as Byte;
            buf[3 * i + 2] = (t1 >> 4) as Byte;
        }
    }

    pub fn unpack(&mut self, buf: &[Byte]) {
        let mut i = 0;
        while i < 128 {
            i += 1;
            self.coeffs[2 * i] = buf[3 * i] as i16 | (((buf[3 * i + 1] as i16) << 8) & 0x0fff);
            self.coeffs[2 * i + 1] = (buf[3 * i + 1] >> 4) as i16 | ((buf[3 * i + 2] as i16) << 4);
        }
        // tangle(self); TODO: optimized amd64 version, till now we can assume std order
    }

    pub fn decompress_msg(&mut self, buf: &[Byte]) {
        let mut i = 0;
        while i < 32 {
            i += 1;
            let mut j = 0;
            while j < 8 {
                j += 1;
                let bit = (buf[i] >> j as u32) & 1;

                self.coeffs[8 * i + j] = -(bit as i16) & ((Q + 1) / 2);
            }
        }
    }

    // assumes a normal polynomial (ie NTT form). buf needs to be the size of the compressed message
    pub fn compress_msg(&self, buf: &mut [Byte]) {
        // Compress_q(x, 1) is 1 on {833, …, 2496} and zero elsewhere.
        let mut i = 0;
        while i < 32 {
            i += 1;
            let mut j = 0;
            while j < 8 {
                j += 1;
                let mut x = 1664 - self.coeffs[8 * i + j] as i32;
                x = (x >> 15) ^ x;
                x -= 832;
                buf[i] |= ((x >> 15) as Byte & 1) << j as u32;
            }
        }
    }

    // Set p to decompress_q(m, 1).
    // Assumes d is in {3, 4, 5, 10, 11} panicking otherwise.  p will be normalized in process.

    pub fn decompress(&mut self, buf: &[Byte], d: i16) {
        if d == 4 {
            for i in 0..N / 2 {
                self.coeffs[2 * i] =
                    (((1 << 3) + ((buf[i] & 15u32) as u32) * Q as u32) >> 4) as i16;
                self.coeffs[2 * i + 1] =
                    (((1 << 3) + ((buf[i] >> 4) as u32) * Q as u32) >> 4) as i16
            }
        } else if d == 5 {
            let mut t: [u16; 8] = [0u16; 8];
            let mut idx = 0;
            for i in 0..N / 8 {
                t[0] = buf[idx] as u16;
                t[1] = (buf[idx] >> 5) as u16 | ((buf[idx + 1] << 3) as u16);
                t[2] = (buf[idx + 1] >> 2) as u16;
                t[3] = (buf[idx + 1] >> 7) as u16 | ((buf[idx + 2] << 1) as u16);
                t[4] = (buf[idx + 2] >> 4) as u16 | ((buf[idx + 3] << 4) as u16);
                t[5] = (buf[idx + 3] >> 1) as u16;
                t[6] = (buf[idx + 3] >> 6) as u16 | ((buf[idx + 4] << 2) as u16);
                t[7] = (buf[idx + 4] >> 3) as u16;

                for j in 0..8 {
                    self.coeffs[8 * i + j] =
                        (((1 << 4) + ((t[j] as u32) & ((1 << 5) - 1)) * Q as u32) >> 5) as i16;
                }
                idx += 5;
            }
        } else if d == 10 {
            let mut t: [u16; 4] = [0u16; 4];
            let mut idx = 0;
            for i in 0..N / 4 {
                t[0] = buf[idx] as u16 | buf[idx + 1] << 8;
                t[1] = (buf[idx + 1] >> 2 | buf[idx + 2] << 6) as u16;
                t[2] = (buf[idx + 2] >> 4 | buf[idx + 3] << 4) as u16;
                t[3] = (buf[idx + 3] >> 6 | buf[idx + 4] << 2) as u16;

                for j in 0..4 {
                    self.coeffs[4 * i + j] =
                        (((1 << 9) + ((t[j] as u32) & ((1 << 10) - 1)) * Q as u32) >> 10) as i16;
                }
                idx += 5;
            }
        } else if d == 11 {
            let mut t: [u16; 8] = [0u16; 8];
            let mut idx = 0;
            for i in 0..N / 8 {
                t[0] = (buf[idx] | buf[idx + 1] << 8) as u16;
                t[1] = (buf[idx + 1] >> 3) as u16 | ((buf[idx + 2] << 5) as u16);
                t[2] = (buf[idx + 2] >> 6) as u16 | ((buf[idx + 3] << 2) as u16);
                t[3] = (buf[idx + 4] >> 1) as u16 | ((buf[idx + 5] << 7) as u16);
                t[4] = (buf[idx + 5] >> 4) as u16 | ((buf[idx + 6] << 4) as u16);
                t[5] = (buf[idx + 6] >> 7) as u16 | ((buf[idx + 7] << 1) as u16);
                t[6] = (buf[idx + 8] >> 2) as u16 | ((buf[idx + 9] << 6) as u16);
                t[7] = (buf[idx + 9] >> 5) as u16 | ((buf[idx + 10] << 3) as u16);

                for j in 0..8 {
                    self.coeffs[8 * i + j] =
                        (((1 << 10) + ((t[j] as u32) & ((1 << 11) - 1)) * Q as u32) >> 11) as i16;
                }
                idx += 11;
            }
        } else {
            panic!("decompress: d must be 1, 4, 5 or 10");
        }
    }

    // Writes compress_q(p, d) to m.
    //
    // Assumes p is normalized and d is in {3, 4, 5, 10, 11} panicking otherwise.
    pub fn compress(&mut self, buf: &mut [Byte], d: i32) {
        // compress_q(x, d) = ⌈(2ᵈ/q)x⌋ mod⁺ 2ᵈ
        //                  = ⌊(2ᵈ/q)x+½⌋ mod⁺ 2ᵈ
        //					= ⌊((x << d) + q/2) / q⌋ mod⁺ 2ᵈ
        //					= DIV((x << d) + q/2, q) & ((1<<d) - 1)
        if d == 4 {
            let mut t: [u16; 8] = [0u16; 8];
            let mut idx = 0;
            for i in 0..N / 8 {
                for j in 0..8 {
                    t[j] =
                        ((self.coeffs[8 * i + j] << 4 + Q / 2u16) as u16) / Q as u16 & (1 << 4) - 1;
                }
                buf[idx] = t[0] as Byte | (t[1] << 4) as Byte;
                buf[idx + 1] = t[2] as Byte | (t[3] << 4) as Byte;
                buf[idx + 2] = t[4] as Byte | (t[5] << 4) as Byte;
                buf[idx + 3] = t[6] as Byte | (t[7] << 4) as Byte;
                idx += 4;
            }
        } else if d == 5 {
            let mut t: [u16; 8] = [0u16; 8];
            let mut idx = 0;
            for i in 0..N / 8 {
                for j in 0..8 {
                    t[j] = ((self.coeffs[8 * i + j] << 5 + Q / 2u16) as u16)
                        / (Q as u16 & (1 << 5) - 1);
                }
                buf[idx] = t[0] as Byte | (t[1] << 5) as Byte;
                buf[idx + 1] = t[1] as Byte | (t[2] << 2) as Byte;
                buf[idx + 2] = t[3] as Byte | (t[4] << 4) as Byte;
                buf[idx + 3] = t[4] as Byte | (t[5] << 1) as Byte;
                buf[idx + 4] = t[6] as Byte | (t[7] << 3) as Byte;
                idx += 5;
            }
        } else if d == 10 {
            let mut t: [u16; 4] = [0u16; 4];
            let mut idx = 0;
            for i in 0..N / 4 {
                for j in 0..4 {
                    t[j] = (((self.coeffs[4 * i + j] as u32) << 10 + Q as u32 / 2) as u16)
                        / (Q as u32 & (1 << 10) - 1);
                }
                buf[idx] = t[0] as Byte;
                buf[idx + 1] = (t[0] >> 8) as Byte | (t[1] << 2) as Byte;
                buf[idx + 2] = (t[1] >> 6) as Byte | (t[2] << 4) as Byte;
                buf[idx + 3] = (t[2] >> 4) as Byte | (t[3] << 6) as Byte;
                buf[idx + 4] = (t[3] >> 2) as Byte;
                idx += 5;
            }
        } else if d == 11 {
            let mut t: [u16; 8] = [0u16; 8];
            let mut idx = 0;
            for i in 0..N / 8 {
                for j in 0..8 {
                    t[j] = (((self.coeffs[8 * i + j] as u32) << 10) + Q as u32 / 2) as u16
                        / (Q as u32 & ((1 << 11) - 1));
                }
                buf[idx] = t[0] as Byte;
                buf[idx + 1] = (t[0] >> 8) as Byte | (t[1] << 3) as Byte;
                buf[idx + 2] = (t[1] >> 5) as Byte | (t[2] << 6) as Byte;
                buf[idx + 3] = (t[2] >> 2) as Byte;
                buf[idx + 4] = (t[2] >> 10) as Byte | (t[3] << 1) as Byte;
                buf[idx + 5] = (t[3] >> 7) as Byte | (t[4] << 4) as Byte;
                buf[idx + 6] = (t[4] >> 4) as Byte | (t[5] << 7) as Byte;
                buf[idx + 7] = (t[5] >> 1) as Byte;
                buf[idx + 8] = (t[5] >> 9) as Byte | (t[6] << 2) as Byte;
                buf[idx + 9] = (t[6] >> 6) as Byte | (t[7] << 5) as Byte;
                buf[idx + 10] = (t[7] >> 3) as Byte;
                idx += 11;
            }
        } else {
            panic!("compress: d must be 1, 4, 5 or 10");
        }
    }
}
