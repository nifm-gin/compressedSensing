clear;
close all;  
clc;


% Modification - Rev 2 - 2021.04.19 -----------------------------------------------------------------------------------------------------------------------------------
% Add BART v0.6.00 path
% set Matlab path and TOOLBOX_PATH environment variable
%bartPath = '/opt/bart-0.6.00';
bartPath = '~/data_ssd/Apps/BART/bart-0.6.00';
addpath(fullfile(bartPath, 'matlab'));
setenv('TOOLBOX_PATH', bartPath);
bart('version')

%load W4_rawdata
%load ('/data_local/data_hdd/alvesd/Pictures/2021.03.10_Test_AF-2x_pics-S/2021.01.20_Mouse_F286_3D-SE-DTI_CS-2x-unif-dens_Matlab-workspace.mat','csUndersampledKSpace')
load ('/data_local/data_hdd/alvesd/Pictures/2021.03.10_Test_AF-2x_pics-S/2021.01.20_Mouse_F286_3D-SE-DTI_CS-2x-unif-dens_Matlab-workspace.mat','fullySampledKSpace')

%d = csUndersampledKSpace;
%d = fullySampledKSpace;

nVolumes = size(fullySampledKSpace,5);
nDWImages = 30;
nB0Images = 3;

globalAF = 2;

AF = nDWImages/(nVolumes/globalAF-nB0Images);
% ---------------------------------------------------------------------------------------------------------------------------------------------------------------------

%AF = 2;
CenSize = 10;
snr = 20;
PC = 10;
TH = 6;

addpath(pwd,'KernelLib' );
addpath(pwd, 'ReconLib' );
addpath(pwd, 'PreImLib' );
addpath(pwd, 'FFTLib' );
addpath('/cbil1/czhang46/Data');


% Slice values to view results ---
whichSliceX = 90;
whichSliceY = 45;
whichSliceZ = 57;
whichChannel = 1;
whichVolume = 1;
% --------------------------------

% Parameters for phase correction step ---
%phaseCorrectionValue = 1.75;
% ----------------------------------------

% Performing phase correction -------------------------------------------------------------------------------------------------------------------------------------------------
%[fullySampledKSpace, recoImageCorrected] =  PhaseCorrection (fullySampledKSpace, 2, phaseCorrectionValue, 1, whichSliceX, whichSliceY, whichSliceZ, whichChannel, whichVolume);
% -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------

% Making k-space for visualization ------------------------
%fullySampledKSpaceView = squeeze(bart('rss 8', fullySampledKSpace));
% ---------------------------------------------------------

% Creating image ----------------------------------------
fullySampledImage = bart('fft -i 7', fullySampledKSpace);
fullySampledImageView = bart('rss 8', fullySampledImage);
fullySampledImageView = squeeze(fullySampledImageView);
% -------------------------------------------------------

% Correcting image with FFTshift -------------------------
fullySampledImageView = fftshift(fullySampledImageView,1);
% --------------------------------------------------------

% ---------------------------------------------------------------------------------------------------------------------------------------------------------------------
%%
initialize;

Timer3 = tic;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
itrNum = zeros(1,ImgSize(3));
disp( '1 -- Initialization is DONE!' );

disp( '==============================================================' );
for xdim = 1:ImgSize(3)
    clear KerPara
disp( strcat('2 -- Kernel Eigen Decomposition on Training ......',num2str(xdim)));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%TrainingDataKPCADTI%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp( KernelMode, 'Gauss' )
    KerPara.Mode   = 'Gauss';
    KerPara.c      = GaussSigma;
    KerPara.d      = NaN;
elseif strcmp( KernelMode, 'Poly' )
    KerPara.Mode   = 'Poly';
    KerPara.c      = PolyConst;
    KerPara.d      = PolyPower;
    KerPara.STh     = STh;
    KerPara.SThDec  = SThDec;
else
end

KerPara.ItrNum  = KPCAItrNum;
KerPara.RegPara = KPCARegVal;
KerPara.ALMPara = ALMRegVal;
KerPara.LamdaReg= LamdaReg;

PreData =  squeeze(abs(LowSos(:,:,xdim,:)));
ReconSize = size(PreData);
MSEVec  = zeros( 1, MaxItrNum );

B = PreData;%zeros( ImgSize );  %%% Comment by UN This is LowResImg
U = zeros(ReconSize);
TraData = PreData;
TraData = permute( TraData, [ 3 1 2 ] );
TraData = reshape( TraData, ReconSize(3), ReconSize(1)*ReconSize(2));
RandIdx = randperm(ReconSize(1)*ReconSize(2), KPCATraNum );
TraData = TraData(:,RandIdx );

Timer1=tic;
% ------ KPCA for training data ------
[ ~, Alpha, D, Kxx, Kxc, ~, ~ ] = KPCA( TraData, TraData, KerPara );
KEigTime=toc(Timer1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp( '2 -- Kernel Eigen Decomposition on Training is DONE!' );
disp( '==============================================================' );
MSEMat=zeros(length(PC),length(TH));
TimerMat=zeros(length(PC),length(TH));
KPCAPctNum=PC;
SThDec=TH;
KerPara.SThDec  = SThDec;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%ReconDTI%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Alpha1= Alpha( :, 1 : KPCAPctNum );
D1= D( 1 : KPCAPctNum );
% ------ projectino for test data, compute Gamma ------
KerPara.Alpha = Alpha1;
KerPara.D     = D1;      
KerPara.Kxx   = Kxx;
KerPara.Kxc   = Kxc;  
MSEVec=zeros(1,500);

clear CurDataPrev;
Timer2 = tic;
Mask = squeeze(Mask3D(:,:,xdim,:,:));
CurData=abs(IFFT2_4D(squeeze(sampledKData(:,:,xdim,:,:)), Mask, 1 ));
ReconSize = size(squeeze(sos(CurData,3)));
RefChanData = squeeze(sos(abs(RefData(:,:,xdim,:,:)),4));
KData = squeeze(sampledKData(:,:,xdim,:,:));
for ItrIndex = 1 : 500 % MaxItrNum   
   %  fprintf( 'Main iteration: %d ......\n', ItrIndex );

%% ==========================  Kernel subproblem ( Gamma and x-subproblem ) ==========================
    BlkData = squeeze(sos(abs(CurData),3));
    TpfData = permute( BlkData, [ 3 1 2 ] );
    TpfData = reshape( TpfData, ReconSize(3), (ReconSize(1)*ReconSize(2)) );   
    [ PrjCoff, ~, ~, ~, ~, Kxy, Kyc ] = KPCA( TpfData, TraData, KerPara );
    % Using Sth like in ktSense
    MaxPjC=max(abs(PrjCoff(:)));
    minPrjCoff=min(abs(PrjCoff(:)));
    %PrjCoff=(abs(PrjCoff)-KerPara.SThDec).*PrjCoff./abs(PrjCoff).*(abs(PrjCoff)>KerPara.SThDec);
    PrjCoff= wthresh(PrjCoff,'h',KerPara.SThDec);
    

    Gamma = KerPara.Alpha * PrjCoff;
    Gamma = Gamma - repmat( mean( Gamma, 1 ), KPCATraNum, 1 ) + 1 / KPCATraNum;
    CurDataPrev=sos(CurData,3);
    %% PreImage Problem. %% Subthresholding is used there within 

    PreIm=FastPreImPol(TraData,Gamma,KerPara);  
    PreIm=reshape(PreIm, ReconSize(3), ReconSize(1), ReconSize(2));
    PreIm=permute(PreIm,[2,3,1]);
    PreIm=abs(PreIm);
     MaxPreIm=max(PreIm(:));
     PreIm=PreIm./max(PreIm(:))*max(abs(CurDataPrev(:)));


    %  CurDataPrev=sos(CurData,3); % This is used to calculate the MSE 
    UpdateImg = zeros(ImgSize(1),ImgSize(2),ImgSize(4),ImgSize(5));
      for frame = 1:ImgSize(5)
          for coil = 1:ImgSize(4)
              UpdateImg(:,:,coil,frame) = PreIm(:,:,frame).*LowSen(:,:,xdim,coil,frame);
          end
      end
      
      UpData = FFT2_4D( UpdateImg, ones( ImgSize(1), ImgSize(2),ImgSize(4),ImgSize(5)), 1 );
      Y=UpData;
      Y = Y.*(Mask==0) + KData; 
      CurData=IFFT2_4D( Y, ones( ImgSize(1), ImgSize(2),ImgSize(4),ImgSize(5)), 1 );
      CurData(isnan(CurData))=0;
 
  
   
   ErrData = abs(RefChanData)-abs(CurDataPrev);
   MSEVec(ItrIndex) = mean( abs( ErrData(:) ).^2 );
   aa=MSEVec(ItrIndex);

    temp = abs(squeeze(sos(abs(CurData),3))-abs(CurDataPrev));
    DiffCur=norm(temp(:));
    if(ItrIndex >2)
        if(MSEVec(ItrIndex-1)< MSEVec(ItrIndex)|| DiffCur<2.8e-15)
            aa=MSEVec(ItrIndex-1);           
            break;
        end
    end
end%%%  End MaxItrNum


PrjTime=toc(Timer2);
ItrTime=PrjTime;
reconImg(:,:,:,:,xdim) = abs(CurData);
end

reconImg = permute(reconImg,[5,1,2,3,4]);
recon = permute(squeeze(sos(NDWI,4)),[3,1,2,4,5]);

% Modification - Rev 2 - 2021.04.19 --------------------------
%recon(:,:,:,3:(ImgSize(5)+2)) = squeeze(sos(abs(reconImg),4));
recon(:,:,:,4:(ImgSize(5)+3)) = squeeze(sos(abs(reconImg),4));
% ------------------------------------------------------------

totalRecoTime = toc(Timer3)

recon = fftshift(recon,1);

% ------------------------------------------------------------

%% Save as NifTi
% Modification - Rev 2 - 2021.04.19 -----------------------------------------------------------------------------------------------------------
Mask3DAllVolumes = ones(size(Mask3D,1),size(Mask3D,2),size(Mask3D,3),size(Mask3D,4),nVolumes);

Mask3DAllVolumes(:,:,:,:,(nB0Images+1):end) = Mask3D;

filename = strcat('W4__Recon3D_overall-AF',num2str(globalAF));
save(filename,'recon','AF','globalAF','Mask3D','Mask3DAllVolumes');
disp('Workspace saved.')

niftiwrite(single(abs(recon)),strcat('KLR-CS-reco_lambReg-',num2str(LamdaReg),'_ALMReg-',num2str(ALMRegVal),'_itmax-',num2str(MaxItrNum),'_overall-AF-',num2str(globalAF),'.nii'));
disp('Image saved.')

% ---------------------------------------------------------------------------------------------------------------------------------------------



%% Correct offset

%pixelOffset = 11;
%reconCorrected = zeros(size(recon,1),size(recon,2),size(recon,3),size(recon,4));
%reconCorrected(:,(pixelOffset:end),:,:) = recon(:,1:(size(reconCorrected,2)-pixelOffset+1),:,:);
%reconCorrected(:,1:(pixelOffset-1),:,:) = recon(:,(size(reconCorrected,2)-pixelOffset+2):end,:,:);

%filename = strcat('W4__Recon3DCorrected_AF',num2str(AF));
%save(filename,'reconCorrected','AF');
%niftiwrite(single(abs(reconCorrected)),strcat('KLR-CS-reco_lamb1-',num2str(LamdaReg),'_lamb2-',num2str(ALMRegVal),'_itmax-',num2str(MaxItrNum),'_Corr-offset.nii'));


%% Show results
whichVolume = 30;


figure('units', 'normalized', 'outerposition', [0 0 1 1]);
subplot(3,2,1)
imagesc(squeeze(fullySampledImageView(:,:,whichSliceZ,whichVolume)))
title('FS')
subplot(3,2,2)
imagesc(squeeze(recon(:,:,whichSliceZ,whichVolume)))
title('KLR CS')
%subplot(3,3,3)
%imagesc(squeeze(reconCorrected(:,:,whichSliceZ,whichVolume)))
%title('KLR-CS after offset correction')

subplot(3,2,3)
imagesc(squeeze(fullySampledImageView(:,whichSliceY,:,whichVolume)))
subplot(3,2,4)
imagesc(squeeze(recon(:,whichSliceY,:,whichVolume)))
%subplot(3,3,6)
%imagesc(squeeze(reconCorrected(:,whichSliceY,:,whichVolume)))

subplot(3,2,5)
imagesc(squeeze(fullySampledImageView(whichSliceX,:,:,whichVolume)))
subplot(3,2,6)
imagesc(squeeze(recon(whichSliceX,:,:,whichVolume)))
%subplot(3,3,9)
%imagesc(squeeze(reconCorrected(whichSliceX,:,:,whichVolume)))

colormap gray
sgtitle(['Comparison FS vs KLR CS - Volume: ',num2str(whichVolume)])
% -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
