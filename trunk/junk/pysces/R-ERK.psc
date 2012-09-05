R1:
RetStar + Shc = RetStar_Shc
k1*RetStar*Shc - k2*RetStar_Shc

R2:
RetStar_Shc = RetStar_ShcStar
k3*RetStar_Shc - k4*RetStar_ShcStar


R3:
RetStar_ShcStar = RetStar + ShcStar
k5*RetStar_ShcStar - k6*RetStar*ShcStar


R4:
ShcStar + GrbHalf = ShcStar_GrbHalf
k7*ShcStar*GrbHalf - k8*ShcStar_GrbHalf


R5:
ShcStar_GrbHalf + Sos = ShcStar_GrbHalf_Sos
k9*ShcStar_GrbHalf*Sos - k10*ShcStar_GrbHalf_Sos


R6:
ShcStar_GrbHalf_Sos = ShcStar_GrbHalf_SosStar
k11*ShcStar_GrbHalf_Sos - k12*ShcStar_GrbHalf_SosStar


R7:
ShcStar_GrbHalf_Sos + Ras_GDP = ShcStar_GrbHalf_SosStar_Ras_GDP
k13*ShcStar_GrbHalf_SosStar*Ras_GDP - k14*ShcStar_GrbHalf_SosStar_Ras_GDP


R8:
ShcStar_GrbHalf_Sos + Ras_GDP = ShcStar_GrbHalf_SosStar_Ras_GTP
k15*ShcStar_GrbHalf_SosStar*Ras_GDP - k16*ShcStar_GrbHalf_SosStar_Ras_GTP


R9:
ShcStar_GrbHalf_Sos + Ras_GDP = ShcStar_GrbHalf_SosStar + Ras_GTP
k17*ShcStar_GrbHalf_SosStar*Ras_GDP - k18*ShcStar_GrbHalf_SosStar*Ras_GTP


R10:
Ras_GTP + Raf = Ras_GTP_Raf
k19*Ras_GTP*Raf - k20*Ras_GTP_Raf


R11:
Ras_GTP_Raf = Ras_GDP_pRaf
k21*Ras_GTP_Raf - k22*Ras_GDP_pRaf


R12:
Ras_GDP_pRaf = Ras_GDP + pRaf
k23*Ras_GDP_pRaf - k24*Ras_GDP*pRaf


R13:
pRaf + MEKHalf = pRaf_MEKHalf
k25*pRaf*MEKHalf - k26*pRaf_MEKHalf


R14:
pRaf_MEKHalf = Rafpp_MEKHalf
k27*pRaf_MEKHalf - k28*Rafpp_MEKHalf


R15:
Rafpp_MEKHalf = Raf + ppMEKHalf
k29*Rafpp_MEKHalf - k30*Raf*ppMEKHalf


R16:
ppMEKHalf + ERKHalf = ppMEKHalf_ERKHalf
k31*ppMEKHalf*ERKHalf - k32*ppMEKHalf_ERKHalf


R17:
ppMEKHalf_ERKHalf = pMEKHalf_pERKHalf
k33*ppMEKHalf_ERKHalf - k34*pMEKHalf_pERKHalf


R18:
pMEKHalf_pERKHalf = MEKHalf_ppERKHalf
k35*pMEKHalf_pERKHalf - k36*MEKHalf_ppERKHalf


R19:
MEKHalf_ppERKHalf = ppERKHalf + MEKHalf
k37*MEKHalf_ppERKHalf - k38*ppERKHalf*MEKHalf



s1 RetStar
s2 Shc 
s3 RetStar_Shc
s4 RetStar_ShcStar
s5 ShcStar
s6 GrbHalf 
s7 ShcStar_GrbHalf 
s8 Sos
s9 ShcStar_GrbHalf_Sos
s10 ShcStar_GrbHalf_SosStar
s11 Ras_GDP
s12 ShcStar_GrbHalf_SosStar_Ras_GDP
s13 ShcStar_GrbHalf_SosStar_Ras_GTP
s14 Ras_GTP
s15 Raf 
s16 Ras_GTP_Raf
s17 Ras_GDP_pRaf
s18 pRaf
s19 MEKHalf 
s20 pRaf_MEKHalf
s21 Rafpp_MEKHalf
s22 ppMEKHalf
s23 ERKHalf 
s24 ppMEKHalf_ERKHalf
s25 pMEKHalf_pERKHalf
s26 MEKHalf_ppERKHalf
s27 ppERKHalf 


