def andGate(y, t, fpra, rpra, kma, kpa, dma, da, fprb, rprb, kmb, kpb, dmb, db, kspl, fprc, rprc, kmc, kpc, dmc, dc, dic):
	[Ia, pA, pAIa, mA, A, Ib, pB, pBIb, mB, B, Ic, w, pC, pCIc, mC, C] = y
	Ia_ = -1 * fpra * Ia * pA + rpra * pAIa
	pA_ = -1 * fpra * Ia * pA + rpra * pAIa
	pAIa_ = fpra * Ia * pA + -1 * rpra * pAIa + -1 * kma * pAIa + kma * pAIa
	mA_ = kma * pAIa + -1 * kpa * mA + kpa * mA + -1 * dma * mA
	A_ = kpa * mA + -1 * da * A + -1 * kspl * A * B
	Ib_ = -1 * fprb * Ib * pB + rprb * pBIb
	pB_ = -1 * fprb * Ib * pB + rprb * pBIb
	pBIb_ = fprb * Ib * pB + -1 * rprb * pBIb + -1 * kmb * pBIb + kmb * pBIb
	mB_ = kmb * pBIb + -1 * kpb * mB + kpb * mB + -1 * dmb * mB
	B_ = kpb * mB + -1 * db * B + -1 * kspl * A * B
	Ic_ = kspl * A * B + -1 * fprc * Ic * pC + rprc * pCIc + -1 * dic * Ic
	w_ = kspl * A * B
	pC_ = -1 * fprc * Ic * pC + rprc * pCIc
	pCIc_ = fprc * Ic * pC + -1 * rprc * pCIc + -1 * kmc * pCIc + kmc * pCIc
	mC_ = kmc * pCIc + -1 * kpc * mC + kpc * mC + -1 * dmc * mC
	C_ = kpc * mC + -1 * dc * C
	return [Ia_, pA_, pAIa_, mA_, A_, Ib_, pB_, pBIb_, mB_, B_, Ic_, w_, pC_, pCIc_, mC_, C_]