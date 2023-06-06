function kData = FFT(im,mask)
  %  kData = fft(fftshift(fft(fftshift(fft(fftshift(im,3),[],3),2),[],2),1),[],1).*mask;
  kData = fftshift(fft(fftshift(fft(fftshift(fft(im,[],3),3),[],2),2),[],1),1).*mask;
end