use std::arch::asm;

//var a [8]VecVirtual
// 	var b [8]VecVirtual
// 	for i := 0; i < 8; i++ {
// 		a[i] = YMM()
// 		b[i] = YMM() // VEX.256.NP.0F38.W0 can use this opcode
// 	}
//
// 	for j := 0; j < 2; j++ {
// 		for i := 0; i < 8; i++ {
// 			VMOVDQU(Mem{Base: aPtr, Disp: 32 * (8*j + i)}, a[i])
// 		}
// 		for i := 0; i < 8; i++ {
// 			VMOVDQU(Mem{Base: bPtr, Disp: 32 * (8*j + i)}, b[i])
// 		}
// 		for i := 0; i < 8; i++ {
// 			VPADDW(a[i], b[i], b[i])
// 		}
// 		for i := 0; i < 8; i++ {
// 			VMOVDQU(b[i], Mem{Base: pPtr, Disp: 32 * (8*j + i)})
// 		}
// 	}

// #[cfg(target_feature = "avx2")]
pub unsafe fn tangle(p: &Poly) {
    // init registers
    let mut a0: i256;

     asm!{

     }
}

#[cfg(target_feature = "avx2")]
pub unsafe fn untangle(p: &Poly) {
    // TODO: Implement
    // asm!{
    //
    // }
}
