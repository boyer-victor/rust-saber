use std::arch::asm;

#[cfg(target_feature = "avx2")]
pub unsafe fn tangle(p: &Poly) {
    // TODO: Implement
    // asm!{
    //
    // }
}

#[cfg(target_feature = "avx2")]
pub unsafe fn untangle(p: &Poly) {
    // TODO: Implement
    // asm!{
    //
    // }
}
