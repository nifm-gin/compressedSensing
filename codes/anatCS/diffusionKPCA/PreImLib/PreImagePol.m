%%%%% Polynomia PreImage Function%%%%%%
%%%% 11/03/2014%%%%%%
%%%%%%% Ukash Nakarmi%%%%%%%

% 
% function PreIM=PreImagePol(TraData,TpfData,Gamma,KerPara)
% 
% [nn,mm]=size(TraData);
% tt= size(TpfData,2);
% Order=KerPara.d;
% c=KerPara.c;
% ZZ=zeros(nn,tt);
% for jj=1:tt
%     GammaVec=Gamma(:,jj);
%     for ii=1:nn
%         TraVec=TraData(ii,:)';
%         Kxe(ii)=GammaVec'*TraVec;
%     end
%     KxeInv=PolyKerInv(Kxe,Order,c);
%     KxeInvPrm=KxeInv';
%     ZZ(:,jj)=KxeInvPrm;   
% end
% PreIM=ZZ;
% end


function PreIM=PreImagePol(TraData,TpfData,Gamma,KerPara)

[nn,mm]=size(TraData);
tt= size(TpfData,2);
Order=KerPara.d;
c=KerPara.c;
ZZ=zeros(nn,tt);
for ii=1:tt
    GammaVec=Gamma(:,ii);
    for jj=1:nn
        TraVec=TraData(jj,:)';
        TraVec1=(TraVec+KerPara.c).^KerPara.d;
        Kxe(jj)=GammaVec'*TraVec1;
    end
    KxeInv=PolyKerInv(Kxe,Order,c);
    KxeInvPrm=KxeInv';
    ZZ(:,ii)=KxeInvPrm;   
end
PreIM=ZZ;
end





