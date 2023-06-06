disp( '1 -- Initializing ......' );

ref_im = ifftshift(ifft(ifftshift(ifft(ifftshift(ifft(fullySampledKSpace,[],1),1),[],2),2),[],3),3);
ref_im = circshift(ref_im,[0 0 0 0 0]);
ref_sos = squeeze(sos(ref_im,4));


ref_im = permute(ref_im,[2,3,1,4,5]); % recon along ky-kz

ImgData = ref_im;

% Modification - Rev 2 - 2021.04.19 ---
%RefData=abs(ImgData(:,:,:,:,3:end)); 
%NDWI = abs(ImgData(:,:,:,:,1:2));

RefData=abs(ImgData(:,:,:,:,4:end)); 
NDWI = abs(ImgData(:,:,:,:,1:3));
% -------------------------------------

ImgSize = size(RefData);        % Data size in image space, row - column - frame
ref_sos = squeeze(sos(ImgData,4));

DD = sum(sum( abs(RefData), 5),4);  % !!!!!!!!!!!!!!!!!!
DM = max( DD(:) );
MM = ( DD > 0.07 * DM );
MM = imfill( MM, 'holes' );

mask = repmat(MM,[1,1,1,ImgSize(4),ImgSize(5)]);
RefData = RefData.*mask;

RefKData=FFT( RefData, ones( ImgSize ));

IfZeroPad = 0;                % If need to zero pad in k-space
IfShowRaw = 0;                % If show zero padded raw data
IfShowInt = 0;                % If show the initialization result
IfShowCov = 0;                % If show the convergence process

%% Algorithm paramter
% Modification - Rev 2 - 2021.04.19 ----------------------------
%MaxItrNum  = 500;              % total maximum iteration number
MaxItrNum  = 200;              % total maximum iteration number
% --------------------------------------------------------------
KPCAItrNum = 100;              % maximum iteration number for pre-imaging
KPCATraNum = 2000;             % KPCA training data numbers
KPCABlkNum = [ 1, 1 ];         % KPCA block number, must be divided by ImgSize
KPCAPctNum = 10;              % KPCA principal component numbers
KPCARegVal = 1;
ALMRegVal  = 0.01;            % ALM variable regularization paramter
LamdaReg   = 0.01;            % Reg para for Updataing KSpace
KernelMode = 'Poly';          % kernel mode 'Gauss' 'Poly'
GaussSigma = 10;
PolyConst  = 500;
PolyPower  = 3;%3
STh=600;                       % Soft Threhold UN 
SThDec=1;                %Decreasing factor for STh UN


for f = 1:ImgSize(5)
   SamplingMask(:,:,f) = rand_sampleMask4k_ori(ImgSize(1),ImgSize(2),AF,1,15);  
end
Mask3D = permute(repmat(SamplingMask,[1,1,1,ImgSize(3),ImgSize(4)]),[1,2,4,5,3]);

KspcData=RefKData.*Mask3D;
sampledKData = ifft(ifftshift(KspcData,3),[],3);

PhasMask = Mask3D;

if IfZeroPad == 1
    [ SpecData, Mask3D ] = KSpaceZeroFilling( ImgSize, PhasMask, KspcData );
else
    SpecData = KspcData;
    Mask3D   = PhasMask;
end

%% Zero filling reconstruction
ImagData0 = IFFT(SpecData,ones(ImgSize));
ImagData0 = abs( ImagData0 );

%% Low resolution reconstruction
LowData = zeros( ImgSize );
MidValD1=ImgSize(1)/2;
MidValD2=ImgSize(2)/2;
hamming_window_full =  cosine_taper_window(ImgSize(1)/2,ImgSize(2)/2,2,CenSize-2,2,4,2);  
LowData((MidValD1-CenSize/2+1):(MidValD1+CenSize/2),(MidValD2-CenSize/2+1):(MidValD2+CenSize/2),:,:,:)=RefKData((MidValD1-CenSize/2+1):(MidValD1+CenSize/2),(MidValD2-CenSize/2+1):(MidValD2+CenSize/2),:,:,:);
LowImag = abs(IFFT( LowData.*repmat(hamming_window_full,[1,1,ImgSize(3),ImgSize(4),ImgSize(5)]), ones( ImgSize ) ));

LowSos = squeeze(sos(LowImag,4));
for coil = 1:ImgSize(4)
   LowSen(:,:,:,:,coil) = squeeze(LowImag(:,:,:,coil,:))./ LowSos;
end
LowSen = permute(LowSen,[1,2,3,5,4]);