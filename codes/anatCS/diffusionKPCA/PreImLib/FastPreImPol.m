function PreImFast= FastPreImPol(TraData,Gamma,KerPara)
TraMat = TraData';
TraMat = TraMat+KerPara.c;
TraMat = TraMat.^KerPara.d;
KxeMat = Gamma'*TraMat;
% Now each row of KxeMat is Kxe
Order = KerPara.d;
c     = KerPara.c;
KxeMatInv=PolyKerInv(KxeMat,Order,c);
PreImFast=KxeMatInv';