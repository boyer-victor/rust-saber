use super::Q;
pub fn mont_reduce(x: i32) -> i16 {
    let m = (x * 62209) as i16;
    return ((x - m as i32 * Q as i32) as u32) as i16;
}

pub fn to_mont(x: i16) -> i16 {
    // |1353 x| ≤ 1353 2¹⁵ ≤ 13318 q ≤ 2¹⁵ q
    return mont_reduce(x as i32 * 1353); // 1353 = R^2 mod q
}

pub fn barrett_reduce(x: i16) -> i16 {
    // for any x we have x mod q = x - ⌊x/q⌋ q
    // Use 20159/2^36 as an approximation of 1/q
    // but noting that ⌊x 20156/2²⁶⌋ = (20159 x) >> 26

    return x - (((x as i32 * 20159) >> 26) * Q) as i16;
}

pub fn csubq(mut x: i16) -> i16 {
    x -= Q;
    x += (x >> 15) & Q;
    return x;
}
