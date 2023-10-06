use super::N;
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
}
