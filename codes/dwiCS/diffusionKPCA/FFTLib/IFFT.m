function im = IFFT(kData,mask)
    %im = ifftshift(ifft(ifftshift(ifft(ifftshift(ifft(kData.*mask,[],1),1),[],2),2),[],3),3);
    im = ifft(ifftshift(ifft(ifftshift(ifft(ifftshift(kData.*mask,1),[],1),2),[],2),3),[],3);
end