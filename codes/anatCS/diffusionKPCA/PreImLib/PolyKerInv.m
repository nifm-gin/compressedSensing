function A=PolyKerInv(Data,Order,Const)
Res=Data.^(1/Order);
Res=Res-Const;
A=Res;
end