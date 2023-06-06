function [maps] = coilMaps(DATA)


%DATA = fft2c(imgs_6(:,:,:,1));

[sx,sy,Nc] = size(DATA);
ncalib = 24; % use 24 calibration lines to compute compression
ksize = [6,6]; % kernel size


% Threshold for picking singular vercors of the calibration matrix
% (relative to largest singlular value.

eigThresh_1 = 0.02;

% threshold of eigen vector decomposition in image space.
eigThresh_2 = 0.95;

% crop a calibration area
calib = crop(DATA,[ncalib,ncalib,Nc]);

im = ifft2c(DATA);

[k,S] = dat2Kernel(calib,ksize);
idx = max(find(S >= S(1)*eigThresh_1));

[M,W] = kernelEig(k(:,:,:,1:idx),[sx,sy]);

maps = M(:,:,:,end).*repmat(W(:,:,end)>eigThresh_2,[1,1,Nc]);
end