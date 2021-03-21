def xorGate(y, t, fpra, rpra, frna, rrna, kma, kpa, kpic, dma, dmic, da, dic, dia, fprb, rprb, frnb, rrnb, kmb, kpb, dmb, db, dib, fspl1, rspl1, kspl2, fecf, recf, decf, fprc, rprc, frnc, rrnc, kmc, kpc, dmc, dc):
	[Ia, pA, pAIa, RNAp, pAIaRNAp, mA, mIc, A, Ic, Ib, pB, pBIb, pBIbRNAp, mB, B, AxBx, AB, w, ABIc, pC, pCIc, pCIcRNAp, mC, C] = y
	Ia_ = -1 * fpra * Ia * pA + rpra * pAIa + -1 * dia * Ia
	pA_ = -1 * fpra * Ia * pA + rpra * pAIa
	pAIa_ = fpra * Ia * pA + -1 * rpra * pAIa + -1 * frna * pAIa * RNAp + rrna * pAIaRNAp + kma * pAIaRNAp
	RNAp_ = -1 * frna * pAIa * RNAp + rrna * pAIaRNAp + kma * pAIaRNAp + -1 * frnb * pBIb * RNAp + rrnb * pBIbRNAp + kmb * pBIbRNAp + -1 * frnc * pCIc * RNAp + rrnc * pCIcRNAp + kmc * pCIcRNAp
	pAIaRNAp_ = frna * pAIa * RNAp + -1 * rrna * pAIaRNAp + -1 * kma * pAIaRNAp
	mA_ = kma * pAIaRNAp + -1 * kpa * mA + kpa * mA + -1 * dma * mA
	mIc_ = kma * pAIaRNAp + -1 * kpic * mIc + kpic * mIc + -1 * dmic * mIc + kmb * pBIbRNAp + -1 * kpic * mIc + kpic * mIc + -1 * dmic * mIc
	A_ = kpa * mA + -1 * da * A + -1 * fspl1 * A * B + rspl1 * AxBx
	Ic_ = kpic * mIc + -1 * dic * Ic + kpic * mIc + -1 * dic * Ic + -1 * fecf * AB * Ic + recf * ABIc + -1 * fprc * Ic * pC + rprc * pCIc
	Ib_ = -1 * fprb * Ib * pB + rprb * pBIb + -1 * dib * Ib
	pB_ = -1 * fprb * Ib * pB + rprb * pBIb
	pBIb_ = fprb * Ib * pB + -1 * rprb * pBIb + -1 * frnb * pBIb * RNAp + rrnb * pBIbRNAp + kmb * pBIbRNAp
	pBIbRNAp_ = frnb * pBIb * RNAp + -1 * rrnb * pBIbRNAp + -1 * kmb * pBIbRNAp
	mB_ = kmb * pBIbRNAp + -1 * kpb * mB + kpb * mB + -1 * dmb * mB
	B_ = kpb * mB + -1 * db * B + -1 * fspl1 * A * B + rspl1 * AxBx
	AxBx_ = fspl1 * A * B + -1 * rspl1 * AxBx + -1 * kspl2 * AxBx
	AB_ = kspl2 * AxBx + -1 * fecf * AB * Ic + recf * ABIc
	w_ = kspl2 * AxBx
	ABIc_ = fecf * AB * Ic + -1 * recf * ABIc + -1 * decf * ABIc
	pC_ = -1 * fprc * Ic * pC + rprc * pCIc
	pCIc_ = fprc * Ic * pC + -1 * rprc * pCIc + -1 * frnc * pCIc * RNAp + rrnc * pCIcRNAp + kmc * pCIcRNAp
	pCIcRNAp_ = frnc * pCIc * RNAp + -1 * rrnc * pCIcRNAp + -1 * kmc * pCIcRNAp
	mC_ = kmc * pCIcRNAp + -1 * kpc * mC + kpc * mC + -1 * dmc * mC
	C_ = kpc * mC + -1 * dc * C
	return [Ia_, pA_, pAIa_, RNAp_, pAIaRNAp_, mA_, mIc_, A_, Ic_, Ib_, pB_, pBIb_, pBIbRNAp_, mB_, B_, AxBx_, AB_, w_, ABIc_, pC_, pCIc_, pCIcRNAp_, mC_, C_]