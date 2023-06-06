function [SclPreIm,Avg]=PreImWeight(MM,RefData,PreIm)
[a,b]=find(MM==1);
[m,n,z]=size(PreIm);
Avg=zeros(m,n);

SclPreIm=zeros(m,n,z);
for ii=1:length(a)
        tempRatio=0;
        for kk=1:z
           tempRatio=tempRatio+RefData(a(ii),b(ii),kk)/PreIm(a(ii),b(ii),kk);
        end
    Avg(a(ii),b(ii))=tempRatio/z;
end

for ii=1:z
    SclPreIm(:,:,ii)=PreIm(:,:,ii).*Avg;
end