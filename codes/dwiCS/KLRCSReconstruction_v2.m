function [klrCSImageReconstructed] = KLRCSReconstruction_v2 (fullySampledKSpace, Mask3DAllVolumes, lambdaRegArray, ALMRegArray, maxItrArray, correctImage, correctionBeforeReco, fftShiftCorrection, offsetCorrection, B0SamplingStrategy, nB0s, centerNormX, centerNormY, centerNormZ, cubeSize, savePathAndFileName, pathTxtDoc, showResults, whichSliceX, whichSliceY, whichSliceZ)

    fprintf("\n-------------------------------------------------------------\n Starting KLR reconstruction... \n-------------------------------------------------------------\n")    
    
    close all

    if (B0SamplingStrategy=="UndersampledB0")    
        fprintf("\nUndersampled B0 strategy chosen...\n")
        
        for lambdaReg=lambdaRegArray
            for ALMReg=ALMRegArray
                for maxItr=maxItrArray             
                    
                    % --------------------------------------------------------------------------------------------
                    % Apply 3D KLR Compressed Sensing reconstruction
                    % --------------------------------------------------------------------------------------------
                    close all;  
                                  
                           
                    undersampledKSpace = fullySampledKSpace.*Mask3DAllVolumes;
                    fprintf("\n Size undersampledKSpace:")
                    size(undersampledKSpace)
                    
                    fprintf("\n Difference between fully sampled and undersampled k-space (0 if an undersampled k-space is provided as input): %f\n",sum(sum(sum(sum(sum(abs(fullySampledKSpace-undersampledKSpace)))))));
                                                                                
                    
                    %AF = 2;
                    CenSize = 10;
                    %snr = 20;
                    PC = 10;
                    TH = 6;

                    disp( '1 -- Initializing ......' );

                    ref_im = ifftshift(ifft(ifftshift(ifft(ifftshift(ifft(fullySampledKSpace,[],1),1),[],2),2),[],3),3);
                    ref_im = circshift(ref_im,[0 0 0 0 0]);
                    ref_sos = squeeze(sos(ref_im,4));


                    ref_im = permute(ref_im,[2,3,1,4,5]); % recon along ky-kz

                    ImgData = ref_im;



                    RefData=abs(ImgData(:,:,:,:,(nB0s+1):end)); 
                    NDWI = abs(ImgData(:,:,:,:,1:nB0s));
                    % -------------------------------------

                    ImgSize = size(RefData);        % Data size in image space, row - column - frame
                    ref_sos = squeeze(sos(ImgData,4));

                    DD = sum(sum( abs(RefData), 5),4);  % !!!!!!!!!!!!!!!!!!
                    DM = max( DD(:) );
                    MM = ( DD > 0.07 * DM );
                    MM = imfill( MM, 'holes' );

                    mask = repmat(MM,[1,1,1,ImgSize(4),ImgSize(5)]);
                    
                    fprintf("\n Size mask:")
                    size(mask)
                    if (showResults)
                        figure
                        subplot(1,3,1)
                        imshow(squeeze(mask(:,:,whichSliceX,1,1)),[])
                        
                        subplot(1,3,2)
                        imshow(squeeze(mask(:,whichSliceZ,:,1,1)),[]) 
                        
                        subplot(1,3,3)
                        imshow(squeeze(mask(whichSliceY,:,:,1,1)),[]) 
                    end
                    
                    RefData = RefData.*mask;

                    RefKData=FFT( RefData, ones( ImgSize )); 
                    fprintf("\n Size RefKData:")
                    size(RefKData)
                    

                    IfZeroPad = 0;                % If need to zero pad in k-space
                    IfShowRaw = 0;                % If show zero padded raw data
                    IfShowInt = 0;                % If show the initialization result
                    IfShowCov = 0;                % If show the convergence process

                    %% Algorithm paramter
                    % Modification - Rev 2 - 2021.04.19 ----------------------------
                    %MaxItrNum  = 500;              % total maximum iteration number
                    %maxItrNum  = 200;              % total maximum iteration number
                    % --------------------------------------------------------------
                    KPCAItrNum = 100;              % maximum iteration number for pre-imaging
                    KPCATraNum = 2000;             % KPCA training data numbers
                    KPCABlkNum = [ 1, 1 ];         % KPCA block number, must be divided by ImgSize
                    KPCAPctNum = 10;              % KPCA principal component numbers
                    KPCARegVal = 1;
                    %ALMRegVal  = 0.01;            % ALM variable regularization paramter
                    %LamdaReg   = 0.01;            % Reg para for Updataing KSpace
                    KernelMode = 'Poly';          % kernel mode 'Gauss' 'Poly'
                    GaussSigma = 10;
                    PolyConst  = 500;
                    PolyPower  = 3;%3
                    STh=600;                       % Soft Threhold UN 
                    SThDec=1;                %Decreasing factor for STh UN


                    % Modification - Rev 3 - 2021.04.27 --------------------------------------------------------

                    %for f = 1:ImgSize(5)
                    %{
                    for f = 1:nVolumes
                       SamplingMask(:,:,f) = rand_sampleMask4k_ori(ImgSize(1),ImgSize(2),AF,1,15);  
                    end
                    %Mask3D = permute(repmat(SamplingMask,[1,1,1,ImgSize(3),ImgSize(4)]),[1,2,4,5,3]);

                    Mask3DAllVolumes = permute(repmat(SamplingMask,[1,1,1,ImgSize(3),ImgSize(4)]),[1,2,4,5,3]);
                    Mask3D = Mask3DAllVolumes(:,:,:,:,(nB0s+1):end);
                    %}
                    % -----------------------------------------------------------------------------------------
                    
                    Mask3DAllVolumes = permute(Mask3DAllVolumes,[2 3 1 4 5]);
                    fprintf("\n Size Mask3DAllVolumes:")
                    size(Mask3DAllVolumes)
                    
                    Mask3D = Mask3DAllVolumes(:,:,:,:,(nB0s+1):end);
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
                    KerPara.ALMPara = ALMReg;
                    KerPara.LamdaReg= lambdaReg;

                    PreData =  squeeze(abs(LowSos(:,:,xdim,:)));
                    ReconSize = size(PreData);
                    MSEVec  = zeros( 1, maxItr );

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
                    MSEVec=zeros(1,maxItr);

                    clear CurDataPrev;
                    Timer2 = tic;
                    Mask = squeeze(Mask3D(:,:,xdim,:,:));
                    CurData=abs(IFFT2_4D(squeeze(sampledKData(:,:,xdim,:,:)), Mask, 1 ));
                    ReconSize = size(squeeze(sos(CurData,3)));
                    RefChanData = squeeze(sos(abs(RefData(:,:,xdim,:,:)),4));
                    KData = squeeze(sampledKData(:,:,xdim,:,:));
                    for ItrIndex = 1 : maxItr % MaxItrNum   
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

                    disp('  ')
                    disp('KLR CS reconstruction of diff. weighted volumes done!!!')
                    disp('  ')

                    reconImg = permute(reconImg,[5,1,2,3,4]);


                    % Modification - Rev 3 - 2021.04.27 --------------------------
                    recon = permute(squeeze(sos(NDWI,4)),[3,1,2,4,5]);
                    reconB0Norm = permute(squeeze(sos(NDWI,4)),[3,1,2,4,5]);

                    
                    
                    
                    
                    

                    %
                    % CS reconstruction of B0 volumes ---

                    disp('  ')
                    disp('Preparing CS reconstruction of B0 volumes...')
                    disp('  ')

                    % Correct mask orientation ------------------------------
                    %mask3DB0Volumes = permute(Mask3DAllVolumes,[3 1 2 4 5]);
                    %mask3DB0Volumes = mask3DB0Volumes(:,:,:,:,1:nB0s);
                    % -------------------------------------------------------


                    csUndersampledKSpace = undersampledKSpace(:,:,:,:,1:nB0s);


                    % Making k-space for visualization ----------------------------
                    csUndersampledKSpaceView = squeeze(bart('rss 8', csUndersampledKSpace));
                    % -------------------------------------------------------------

                    % Creating image -------------------------------------------
                    csZeroPaddingImage = bart('fft -i 7', csUndersampledKSpace);
                    csZeroPaddingImageView = bart('rss 8', csZeroPaddingImage);
                    csZeroPaddingImageView = squeeze(csZeroPaddingImageView);
                    % ----------------------------------------------------------

                    % Correcting image with FFTshift ---------------------------
                    csZeroPaddingImageView = fftshift(csZeroPaddingImageView,1);
                    % ----------------------------------------------------------
                    
                    whichVolume = 1;
                    if(showResults)
                        close all
                        figure;
                        colormap gray
                        subplot(2,3,1)
                        imagesc(abs(squeeze(csUndersampledKSpaceView(:,:,whichSliceZ,whichVolume))).^0.25)

                        subplot(2,3,2)
                        imagesc(abs(squeeze(csUndersampledKSpaceView(:,whichSliceY,:,whichVolume))).^0.25)

                        subplot(2,3,3)
                        imagesc(abs(squeeze(csUndersampledKSpaceView(whichSliceX,:,:,whichVolume))).^0.25)

                        subplot(2,3,4)
                        imagesc(abs(squeeze(csZeroPaddingImageView(:,:,whichSliceZ,whichVolume))))

                        subplot(2,3,5)
                        imagesc(abs(squeeze(csZeroPaddingImageView(:,whichSliceY,:,whichVolume))))

                        subplot(2,3,6)
                        imagesc(abs(squeeze(csZeroPaddingImageView(whichSliceX,:,:,whichVolume))))
                    end

                    close all
                    pause(2)
                    % --------------------------------------------------------------------------------------------
                    % Coil sensitivity maps
                    % --------------------------------------------------------------------------------------------
                    sensitivityMapsAllVolumes = zeros(size(csUndersampledKSpace,1),size(csUndersampledKSpace,2),size(csUndersampledKSpace,3),size(csUndersampledKSpace,4),size(csUndersampledKSpace,5));

                    fprintf("\n Size csUndersampledKSpace:")
                    size(csUndersampledKSpace)
                    
                    for volume=1:size(csUndersampledKSpace,5)
                        disp(['Calculating coil sensitivity maps... - Volume: ', num2str(volume)])
                        tic
                        KOrigBart_maps = bart('ecalib -c0. -m1', squeeze(csUndersampledKSpace(:,:,:,:,volume)));
                        coilSensMapTime = toc;
                        sensitivityMapsAllVolumes(:,:,:,:,volume) = KOrigBart_maps;

                        sensitivityMapsAllVolumesView = fftshift(sensitivityMapsAllVolumes,1);

                        if(showResults && ((volume==1 || volume==(size(undersampledKSpace,5)))))
                            figure('WindowState','maximized');
                            %suptitle({'Coil sensitivity maps',['Volume: ', num2str(volume)],'',''});
                            sgtitle({'Coil sensitivity maps',['Volume: ', num2str(volume)]});
                            subplot(3,4,1)
                            imagesc(abs(squeeze(sensitivityMapsAllVolumesView(whichSliceX,:,:,1,volume))));
                            title(['YZ plane - Ch: 1 - Slice: ', num2str(whichSliceX)]);
                            %colorbar
                            caxis([0 1])
                            subplot(3,4,2)
                            imagesc(abs(squeeze(sensitivityMapsAllVolumesView(whichSliceX,:,:,2,volume))));
                            title(['YZ plane - Ch: 2 - Slice: ', num2str(whichSliceX)]);
                            %colorbar
                            caxis([0 1])
                            subplot(3,4,3)
                            imagesc(abs(squeeze(sensitivityMapsAllVolumesView(whichSliceX,:,:,3,volume))));
                            title(['YZ plane - Ch: 3 - Slice: ', num2str(whichSliceX)]);
                            %colorbar
                            caxis([0 1])
                            subplot(3,4,4)
                            imagesc(abs(squeeze(sensitivityMapsAllVolumesView(whichSliceX,:,:,4,volume))));
                            title(['YZ plane - Ch: 4 - Slice: ', num2str(whichSliceX)]);
                            colorbar
                            caxis([0 1])

                            subplot(3,4,5)
                            imagesc(abs(squeeze(sensitivityMapsAllVolumesView(:,whichSliceY,:,1,volume))));
                            title(['XZ plane - Ch: 1 - Slice: ', num2str(whichSliceY)]);
                            %colorbar
                            caxis([0 1])
                            subplot(3,4,6)
                            imagesc(abs(squeeze(sensitivityMapsAllVolumesView(:,whichSliceY,:,2,volume))));
                            title(['XZ plane - Ch: 2 - Slice: ', num2str(whichSliceY)]);
                            %colorbar
                            caxis([0 1])
                            subplot(3,4,7)
                            imagesc(abs(squeeze(sensitivityMapsAllVolumesView(:,whichSliceY,:,3,volume))));
                            title(['XZ plane - Ch: 3 - Slice: ', num2str(whichSliceY)]);
                            %colorbar
                            caxis([0 1])
                            subplot(3,4,8)
                            imagesc(abs(squeeze(sensitivityMapsAllVolumesView(:,whichSliceY,:,4,volume))));
                            title(['XZ plane - Ch: 4 - Slice: ', num2str(whichSliceY)]);
                            colorbar
                            caxis([0 1])

                            subplot(3,4,9)
                            imagesc(abs(squeeze(sensitivityMapsAllVolumesView(:,:,whichSliceZ,1,volume))));
                            title(['XY plane - Ch: 1 - Slice: ', num2str(whichSliceZ)]);
                            %colorbar
                            caxis([0 1])
                            subplot(3,4,10)
                            imagesc(abs(squeeze(sensitivityMapsAllVolumesView(:,:,whichSliceZ,2,volume))));
                            title(['XY plane - Ch: 2 - Slice: ', num2str(whichSliceZ)]);
                            %colorbar
                            caxis([0 1])
                            subplot(3,4,11)
                            imagesc(abs(squeeze(sensitivityMapsAllVolumesView(:,:,whichSliceZ,3,volume))));
                            title(['XY plane - Ch: 3 - Slice: ', num2str(whichSliceZ)]);
                            %colorbar
                            caxis([0 1])
                            subplot(3,4,12)
                            imagesc(abs(squeeze(sensitivityMapsAllVolumesView(:,:,whichSliceZ,4,volume))));
                            title(['XY plane - Ch: 4 - Slice: ', num2str(whichSliceZ)]);
                            colorbar
                            caxis([0 1])
                            pause(1)
                            terminusName = '_Cor-Hor-Sag-view';
                            imageExtension = '.jpg';
                            saveas(gcf,strcat(savePathAndFileName, '_0_Coil-sensitivity-maps_Volume-',num2str(volume),  terminusName, imageExtension));
                        end

                    end
                    
                    
                    coilSensitivityMapsCorrected = sensitivityMapsAllVolumes;
                    if (correctImage)
                        fprintf("\nCorrectig coil sensitivity maps using fftshift...\n")
                        if (fftShiftCorrection(1))
                            coilSensitivityMapsCorrected = fftshift(coilSensitivityMapsCorrected,1);
                        end
                        if (fftShiftCorrection(2))
                            coilSensitivityMapsCorrected = fftshift(coilSensitivityMapsCorrected,2);
                        end
                        if (fftShiftCorrection(3))
                            coilSensitivityMapsCorrected = fftshift(coilSensitivityMapsCorrected,3);                                                                
                        end


                        coilSensitivityMapsCorrected = circshift(coilSensitivityMapsCorrected, offsetCorrection);

                    end
                    % Save .nii image in a specific folder ------------------------------------------------------------------------------------
                    niftiwrite((single(abs(coilSensitivityMapsCorrected))),strcat(savePathAndFileName,'_coilSensitivityMaps_lamb1-',num2str(0.005),'_lamb2-',num2str(0.002),'_maxItr-',num2str(maxItr),'.nii'))
                    % -------------------------------------------------------------------------------------------------------------------------


                    % Save .nii image in a specific folder ------------------------------------------------------------------------------------
                    %niftiwrite((single(abs(sensitivityMapsAllVolumesView))),strcat(savePathAndFileNme,'_coilSensitivityMapsAllVolumesView.nii'))
                    % -------------------------------------------------------------------------------------------------------------------------


                    % --------------------------------------------------------------------------------------------
                    % Apply 3D Compressed Sensing reconstruction
                    % --------------------------------------------------------------------------------------------
                    csReconstructedImageAllVolumes = zeros(size(csUndersampledKSpace,1),size(csUndersampledKSpace,2),size(csUndersampledKSpace,3),size(csUndersampledKSpace,5));


                    disp('  ')
                    disp('Starting CS reconstruction of B0 volumes...')
                    disp('  ')

                    for volume=1:size(csUndersampledKSpace,5)
                        disp(['CS image reconstruction... - Volume: ', num2str(volume)])
                        tic
                        csReconstructedImage = bart('pics -l1 -r0.005 -R T:7:0:0.002 -i200', squeeze(csUndersampledKSpace(:,:,:,:,volume)), squeeze(sensitivityMapsAllVolumes(:,:,:,:,volume)));
                        csRecoTime = toc;
                        csReconstructedImageAllVolumes(:,:,:,volume) = csReconstructedImage;

                        if (showResults)
                            figure;
                            subplot(1,3,1)
                            imshow(abs(squeeze(csReconstructedImageAllVolumes(whichSliceX,:,:,volume))),[]);
                            title({'3D Compressed Sensing reconstruction',['Volume: ',num2str(volume)],['Y-Z plan - slice ', num2str(whichSliceX)]});

                            subplot(1,3,2)
                            imshow(abs(squeeze(csReconstructedImageAllVolumes(:,whichSliceY,:,volume))),[]);
                            title({'3D Compressed Sensing reconstruction',['Volume: ',num2str(volume)],['X-Z plan - slice ', num2str(whichSliceY)]});

                            subplot(1,3,3)               
                            imshow(abs(squeeze(csReconstructedImageAllVolumes(:,:,whichSliceZ,volume))),[]);
                            title({'3D Compressed Sensing reconstruction',['Volume: ',num2str(volume)],['X-Y plan - slice ', num2str(whichSliceZ)]});
                        end
                        % Save .nii image for each volume in a specific folder ------------------------------------------------------------------------------------
                        %niftiwrite((single(abs((csReconstructedImage)))),strcat(savePathAndFileNme,'_csReconstructedImage_Volume-',num2str(volume),'_single-abs.nii'))
                        %niftiwrite((single(abs(fftshift(csReconstructedImage,1)))),strcat(savePathAndFileNme,'_csReconstructedImage_Volume-',num2str(volume),'_fftshift-1','_single-abs.nii'))
                        %niftiwrite((single(abs(fftshift(csReconstructedImage,2)))),strcat(savePathAndFileNme,'_csReconstructedImage_Volume-',num2str(volume),'_fftshift-2','_single-abs.nii'))
                        %niftiwrite((single(abs(fftshift(csReconstructedImage,3)))),strcat(savePathAndFileNme,'_csReconstructedImage_Volume-',num2str(volume),'_fftshift-3','_single-abs.nii'))
                        % -------------------------------------------------------------------------------------------------------------------------

                        pause(1)
                    end

                    % -----------------------------------

                    disp('  ')
                    disp('CS reconstruction of B0 volumes done!!!')
                    disp('  ')





                    %niftiwrite((single(abs(fftshift(csReconstructedImageAllVolumes,1)))),strcat(savePathAndFileNme,'_csReconstructedImageB0Volumes_single-abs.nii'))

                    % Normalize B0 reconstructed ---
                    %csReconstructedImageAllVolumesB0Norm = csReconstructedImageAllVolumes;
                    
                    if (correctionBeforeReco)
                        [csReconstructedImageAllVolumesB0Norm, B0RegularizationVector] =  imageNormalizationWrt3DReference (permute(squeeze(sos(NDWI,4)),[3,1,2,4,5]), csReconstructedImageAllVolumes, centerNormX, centerNormY, centerNormZ, cubeSize);
                        %csReconstructedImageAllVolumesB0Norm(:,:,:,1:nB0s) = csReconstructedImageAllVolumesOnlyNormalizedB0;
                        %csReconstructedImageAllVolumesB0Norm = fftshift(csReconstructedImageAllVolumesB0Norm,1);
                    else
                        [csReconstructedImageAllVolumesB0Norm, B0RegularizationVector] =  imageNormalizationWrt3DReference (fftshift(permute(squeeze(sos(NDWI,4)),[3,1,2,4,5]),1), fftshift(csReconstructedImageAllVolumes,1), centerNormX, centerNormY, centerNormZ, cubeSize);
                        %csReconstructedImageAllVolumesB0Norm(:,:,:,1:nB0s) = csReconstructedImageAllVolumesOnlyNormalizedB0;
                        csReconstructedImageAllVolumesB0Norm = fftshift(csReconstructedImageAllVolumesB0Norm,1);
                    end
                    %niftiwrite((single(abs(fftshift(csReconstructedImageAllVolumesB0Norm,1)))),strcat(savePathAndFileNme,'_csReconstructedImageB0VolumesNorm_single-abs.nii'))
                    % ------------------------------



                    % ------------------------------------------------------------

                    % Modification - Rev 2 - 2021.04.19 --------------------------
                    %recon(:,:,:,3:(ImgSize(5)+2)) = squeeze(sos(abs(reconImg),4));
                    recon(:,:,:,1:nB0s) = csReconstructedImageAllVolumes;
                    recon(:,:,:,(nB0s+1):(ImgSize(5)+nB0s)) = squeeze(sos(abs(reconImg),4));

                    reconB0Norm(:,:,:,1:nB0s) = csReconstructedImageAllVolumesB0Norm;
                    reconB0Norm(:,:,:,(nB0s+1):(ImgSize(5)+nB0s)) = squeeze(sos(abs(reconImg),4));
                    % ------------------------------------------------------------

                    totalRecoTime = toc(Timer3)
                    
                    

                    %recon = fftshift(recon,1);
                    %reconB0Norm = fftshift(reconB0Norm,1);

                    % ------------------------------------------------------------

                    % Save as NifTi
                    % Modification - Rev 2 - 2021.04.19 -----------------------------------------------------------------------------------------------------------
                    %Mask3DAllVolumes = ones(size(Mask3D,1),size(Mask3D,2),size(Mask3D,3),size(Mask3D,4),nVolumes);

                    %Mask3DAllVolumes(:,:,:,:,(nB0s+1):end) = Mask3D;

                    %filename = strcat(savePathAndFileNme,'W4__Recon3D_overall-AF',num2str(globalAF),'_ConvCS-B0-norm.mat');
                    %save(filename,'recon','reconB0Norm','AF','globalAF','Mask3D','Mask3DAllVolumes');
                    %disp('Workspace saved.')

                    %niftiwrite(single(abs(recon)),strcat(savePathAndFileNme,'KLR-CS-reco_lambReg-',num2str(lamdaReg),'_ALMReg-',num2str(ALMReg),'_itmax-',num2str(maxItr),'_overall-AF-',num2str(globalAF),'_ConvCS-B0.nii'));
                    %niftiwrite(single(abs(reconB0Norm)),strcat(savePathAndFileNme,'KLR-CS-reco_lambReg-',num2str(lamdaReg),'_ALMReg-',num2str(ALMReg),'_itmax-',num2str(maxItr),'_overall-AF-',num2str(globalAF),'_ConvCS-B0-norm.nii'));
                    %niftiwrite(single(abs(fullySampledImageView)),strcat(savePathAndFileNameFullySampled,'.nii'))

                    klrCSImageReconstructed = reconB0Norm;

                    %[csReconstructedkSpace, csReconstructedImage] =  MultiChannelCSReconstruction (undersampledKSpace, '3D-CS', 1, 1, 1, showResults, whichSliceX, whichSliceY, whichSliceZ);
                    if (correctImage)
                        fprintf("\nCorrecting offset of reconstructed image...\n")
                        % Correcting image with FFTshift -------------------------------------------
                        if (fftShiftCorrection(1))
                            klrCSImageReconstructed = fftshift(klrCSImageReconstructed,1);
                        end
                        if (fftShiftCorrection(2))
                            klrCSImageReconstructed = fftshift(klrCSImageReconstructed,2);
                        end
                        if (fftShiftCorrection(3))
                            klrCSImageReconstructed = fftshift(klrCSImageReconstructed,3);
                        end
                        % --------------------------------------------------------------------------

                        % Correct offset --------------------------------------------------------------------------
                        klrCSImageReconstructed = circshift(klrCSImageReconstructed, offsetCorrection);
                        % -----------------------------------------------------------------------------------------

                    end
                    
                    if (showResults)    
                        figure;
                        whichVolume = 1;
                        subplot(4,3,1)
                        imshow(abs(squeeze(klrCSImageReconstructed(whichSliceX,:,:,whichVolume))),[]);
                        title({'KLR CS (for DWIs) and conv. CS (for B0s) reconstruction',['Volume: ',num2str(whichVolume)],['Y-Z plan - slice ', num2str(whichSliceX)]});

                        subplot(4,3,2)
                        imshow(abs(squeeze(klrCSImageReconstructed(:,whichSliceY,:,whichVolume))),[]);
                        title({'KLR CS (for DWIs) and conv. CS (for B0s) reconstruction',['Volume: ',num2str(whichVolume)],['X-Z plan - slice ', num2str(whichSliceY)]});

                        subplot(4,3,3)               
                        imshow(abs(squeeze(klrCSImageReconstructed(:,:,whichSliceZ,whichVolume))),[]);
                        title({'KLR CS (for DWIs) and conv. CS (for B0s) reconstruction',['Volume: ',num2str(whichVolume)],['X-Y plan - slice ', num2str(whichSliceZ)]});
                        
                        figure;
                        whichVolume = nB0s;
                        subplot(4,3,4)
                        imshow(abs(squeeze(klrCSImageReconstructed(whichSliceX,:,:,whichVolume))),[]);
                        title({'KLR CS (for DWIs) and conv. CS (for B0s) reconstruction',['Volume: ',num2str(whichVolume)],['Y-Z plan - slice ', num2str(whichSliceX)]});

                        subplot(4,3,5)
                        imshow(abs(squeeze(klrCSImageReconstructed(:,whichSliceY,:,whichVolume))),[]);
                        title({'KLR CS (for DWIs) and conv. CS (for B0s) reconstruction',['Volume: ',num2str(whichVolume)],['X-Z plan - slice ', num2str(whichSliceY)]});

                        subplot(4,3,6)               
                        imshow(abs(squeeze(klrCSImageReconstructed(:,:,whichSliceZ,whichVolume))),[]);
                        title({'KLR CS (for DWIs) and conv. CS (for B0s) reconstruction',['Volume: ',num2str(whichVolume)],['X-Y plan - slice ', num2str(whichSliceZ)]});
                        
                        figure;
                        whichVolume = nB0s + 1;
                        subplot(4,3,7)
                        imshow(abs(squeeze(klrCSImageReconstructed(whichSliceX,:,:,whichVolume))),[]);
                        title({'KLR CS (for DWIs) and conv. CS (for B0s) reconstruction',['Volume: ',num2str(whichVolume)],['Y-Z plan - slice ', num2str(whichSliceX)]});

                        subplot(4,3,8)
                        imshow(abs(squeeze(klrCSImageReconstructed(:,whichSliceY,:,whichVolume))),[]);
                        title({'KLR CS (for DWIs) and conv. CS (for B0s) reconstruction',['Volume: ',num2str(whichVolume)],['X-Z plan - slice ', num2str(whichSliceY)]});

                        subplot(4,3,9)               
                        imshow(abs(squeeze(klrCSImageReconstructed(:,:,whichSliceZ,whichVolume))),[]);
                        title({'KLR CS (for DWIs) and conv. CS (for B0s) reconstruction',['Volume: ',num2str(whichVolume)],['X-Y plan - slice ', num2str(whichSliceZ)]});
                        
                        figure;
                        whichVolume = size(klrCSImageReconstructed,4);
                        subplot(4,3,10)
                        imshow(abs(squeeze(klrCSImageReconstructed(whichSliceX,:,:,whichVolume))),[]);
                        title({'KLR CS (for DWIs) and conv. CS (for B0s) reconstruction',['Volume: ',num2str(whichVolume)],['Y-Z plan - slice ', num2str(whichSliceX)]});

                        subplot(4,3,11)
                        imshow(abs(squeeze(klrCSImageReconstructed(:,whichSliceY,:,whichVolume))),[]);
                        title({'KLR CS (for DWIs) and conv. CS (for B0s) reconstruction',['Volume: ',num2str(whichVolume)],['X-Z plan - slice ', num2str(whichSliceY)]});

                        subplot(4,3,12)               
                        imshow(abs(squeeze(klrCSImageReconstructed(:,:,whichSliceZ,whichVolume))),[]);
                        title({'KLR CS (for DWIs) and conv. CS (for B0s) reconstruction',['Volume: ',num2str(whichVolume)],['X-Y plan - slice ', num2str(whichSliceZ)]});
                    end


                    % Save .nii image in a specific folder ------------------------------------------------------------------------------------                
                    pathAndfilename = strcat(savePathAndFileName,'_KLRCS_lambReg-',num2str(lambdaReg),'_ALMReg-',num2str(ALMReg),'_maxItr-',num2str(maxItr),'.nii');
                    niftiwrite((single(abs(klrCSImageReconstructed))),pathAndfilename);
                    fprintf("\nReconstructed image saved as a NIfTI file.\n")
                    % -------------------------------------------------------------------------------------------------------------------------


                    % Pause step 3 ---
                    pause(3)
                    % ----------------

                    fileID = fopen(pathTxtDoc,'a');
                    fprintf(fileID,'%s\r\n',pathAndfilename);
                    fclose(fileID);
                    fprintf("\nPath of reconstructed image added in a txt file.\n")
                end
            end
        end
        
        
        
    elseif (B0SamplingStrategy=="FullySampledB0")
        fprintf("\nFully sampled B0 strategy chosen...\n")
        % --------------------------------------------------------------------------------------------
        % Coil sensitivity maps
        % --------------------------------------------------------------------------------------------
        for lambdaReg=lambdaRegArray
            for ALMReg=ALMRegArray
                for maxItr=maxItrArray
                    lambdaReg
                    ALMReg
                    maxItr

                     % --------------------------------------------------------------------------------------------
                    % Apply 3D KLR Compressed Sensing reconstruction
                    % --------------------------------------------------------------------------------------------
                    close all;  
                                  
                           
                    undersampledKSpace = fullySampledKSpace.*Mask3DAllVolumes;
                    
                    fprintf("\n Size undersampledKSpace:")
                    size(undersampledKSpace)
                    
                    fprintf("\n Difference between fully sampled and undersampled k-space (0 if an undersampled k-space is provided as input): %f\n",sum(sum(sum(sum(sum(abs(fullySampledKSpace-undersampledKSpace)))))));
                                        
                    
                    %AF = 2;
                    CenSize = 10;
                    %snr = 20;
                    PC = 10;
                    TH = 6;

                    disp( '1 -- Initializing ......' );

                    ref_im = ifftshift(ifft(ifftshift(ifft(ifftshift(ifft(fullySampledKSpace,[],1),1),[],2),2),[],3),3);
                    ref_im = circshift(ref_im,[0 0 0 0 0]);
                    ref_sos = squeeze(sos(ref_im,4));


                    ref_im = permute(ref_im,[2,3,1,4,5]); % recon along ky-kz
                    
                    ImgData = ref_im;



                    RefData=abs(ImgData(:,:,:,:,(nB0s+1):end)); 
                    NDWI = abs(ImgData(:,:,:,:,1:nB0s));
                    % -------------------------------------

                    ImgSize = size(RefData);        % Data size in image space, row - column - frame
                    ref_sos = squeeze(sos(ImgData,4));

                    DD = sum(sum( abs(RefData), 5),4);  % !!!!!!!!!!!!!!!!!!
                    DM = max( DD(:) );
                    MM = ( DD > 0.07 * DM );
                    MM = imfill( MM, 'holes' );

                    mask = repmat(MM,[1,1,1,ImgSize(4),ImgSize(5)]);
                    
                    fprintf("\n Size mask:")
                    size(mask)
                    if (showResults)
                        figure
                        subplot(1,3,1)
                        imshow(squeeze(mask(:,:,whichSliceX,1,1)),[])
                        
                        subplot(1,3,2)
                        imshow(squeeze(mask(:,whichSliceZ,:,1,1)),[]) 
                        
                        subplot(1,3,3)
                        imshow(squeeze(mask(whichSliceY,:,:,1,1)),[]) 
                    end
                    
                    %RefData = RefData.*mask;

                    RefKData=FFT( RefData, ones( ImgSize )); 
                    fprintf("\n Size RefKData:")
                    size(RefKData)
                    

                    IfZeroPad = 0;                % If need to zero pad in k-space
                    IfShowRaw = 0;                % If show zero padded raw data
                    IfShowInt = 0;                % If show the initialization result
                    IfShowCov = 0;                % If show the convergence process

                    %% Algorithm paramter
                    % Modification - Rev 2 - 2021.04.19 ----------------------------
                    %MaxItrNum  = 500;              % total maximum iteration number
                    %maxItrNum  = 200;              % total maximum iteration number
                    % --------------------------------------------------------------
                    KPCAItrNum = 100;              % maximum iteration number for pre-imaging
                    KPCATraNum = 2000;             % KPCA training data numbers
                    KPCABlkNum = [ 1, 1 ];         % KPCA block number, must be divided by ImgSize
                    KPCAPctNum = 10;              % KPCA principal component numbers
                    KPCARegVal = 1;
                    %ALMRegVal  = 0.01;            % ALM variable regularization paramter
                    %LamdaReg   = 0.01;            % Reg para for Updataing KSpace
                    KernelMode = 'Poly';          % kernel mode 'Gauss' 'Poly'
                    GaussSigma = 10;
                    PolyConst  = 500;
                    PolyPower  = 3;%3
                    STh=600;                       % Soft Threhold UN 
                    SThDec=1;                %Decreasing factor for STh UN


                    % Modification - Rev 3 - 2021.04.27 --------------------------------------------------------

                    %for f = 1:ImgSize(5)
                    %{
                    for f = 1:nVolumes
                       SamplingMask(:,:,f) = rand_sampleMask4k_ori(ImgSize(1),ImgSize(2),AF,1,15);  
                    end
                    %Mask3D = permute(repmat(SamplingMask,[1,1,1,ImgSize(3),ImgSize(4)]),[1,2,4,5,3]);

                    Mask3DAllVolumes = permute(repmat(SamplingMask,[1,1,1,ImgSize(3),ImgSize(4)]),[1,2,4,5,3]);
                    Mask3D = Mask3DAllVolumes(:,:,:,:,(nB0s+1):end);
                    %}
                    % -----------------------------------------------------------------------------------------
                    
                    Mask3DAllVolumes = permute(Mask3DAllVolumes,[2 3 1 4 5]);
                    fprintf("\n Size Mask3DAllVolumes:")
                    size(Mask3DAllVolumes)
                    
                    undersampledKSpace = permute(undersampledKSpace,[2 3 1 4 5]); % Changed on 2021.08.30
                    size(undersampledKSpace(:,:,:,:,(nB0s+1):end))
                    
                    Mask3D = Mask3DAllVolumes(:,:,:,:,(nB0s+1):end);
                    KspcData=RefKData.*Mask3D; % Changed on 2021.08.30
                    %KspcData=RefKData; % Changed on 2021.08.30
                    
                    
                    
                    KspcData2 = undersampledKSpace(:,:,:,:,(nB0s+1):end); % Changed on 2021.08.30
                    
                    sampledKData = ifft(ifftshift(KspcData,3),[],3);
                    sampledKData2 = ifft(ifftshift(KspcData2,3),[],3);

                    undersampledKSpace = permute(undersampledKSpace,[3 1 2 4 5]); % Changed on 2021.08.30
                    size(undersampledKSpace)
                    
                    PhasMask = Mask3D;

                    if IfZeroPad == 1
                        [ SpecData, Mask3D ] = KSpaceZeroFilling( ImgSize, PhasMask, KspcData );
                    else
                        SpecData = KspcData;
                        SpecData2 = KspcData2;
                        
                        Mask3D   = PhasMask;
                    end

                    %% Zero filling reconstruction
                    ImagData0 = IFFT(SpecData,ones(ImgSize));
                    ImagData0 = abs( ImagData0 );
                    
                    ImagData02 = IFFT(SpecData2,ones(ImgSize));
                    ImagData02 = abs( ImagData02 );

                    %% Low resolution reconstruction
                    LowData = zeros( ImgSize );
                    MidValD1=ImgSize(1)/2;
                    MidValD2=ImgSize(2)/2;
                    hamming_window_full =  cosine_taper_window(ImgSize(1)/2,ImgSize(2)/2,2,CenSize-2,2,4,2);  
                    % Original ---
                    %LowData((MidValD1-CenSize/2+1):(MidValD1+CenSize/2),(MidValD2-CenSize/2+1):(MidValD2+CenSize/2),:,:,:)=RefKData((MidValD1-CenSize/2+1):(MidValD1+CenSize/2),(MidValD2-CenSize/2+1):(MidValD2+CenSize/2),:,:,:);
                    % ---
                    LowData((MidValD1-CenSize/2+1):(MidValD1+CenSize/2),(MidValD2-CenSize/2+1):(MidValD2+CenSize/2),:,:,:)=KspcData((MidValD1-CenSize/2+1):(MidValD1+CenSize/2),(MidValD2-CenSize/2+1):(MidValD2+CenSize/2),:,:,:);
                    LowImag = abs(IFFT( LowData.*repmat(hamming_window_full,[1,1,ImgSize(3),ImgSize(4),ImgSize(5)]), ones( ImgSize ) ));

                    LowSos = squeeze(sos(LowImag,4));
                    for coil = 1:ImgSize(4)
                       LowSen(:,:,:,:,coil) = squeeze(LowImag(:,:,:,coil,:))./ LowSos;
                    end
                    LowSen = permute(LowSen,[1,2,3,5,4]);

                    
                    % Code added on 2021.11.28 to study the evolution of MeanSquareError and kernelDiff ---
                    meanSquareErrorArray = NaN(ImgSize(3), maxItr);
                    diffCurArray = NaN(ImgSize(3), maxItr);
                    nItrArray = zeros(ImgSize(3), 1);
                    % -------------------------------------------------------------------------------------
                    
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
                    KerPara.ALMPara = ALMReg;
                    KerPara.LamdaReg= lambdaReg;

                    PreData =  squeeze(abs(LowSos(:,:,xdim,:)));
                    ReconSize = size(PreData);
                    MSEVec  = zeros( 1, maxItr );

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
                    MSEVec=zeros(1,maxItr);

                    clear CurDataPrev;
                    Timer2 = tic;
                    Mask = squeeze(Mask3D(:,:,xdim,:,:));
                    CurData=abs(IFFT2_4D(squeeze(sampledKData(:,:,xdim,:,:)), Mask, 1 ));
                    ReconSize = size(squeeze(sos(CurData,3)));
                    RefChanData = squeeze(sos(abs(RefData(:,:,xdim,:,:)),4));
                    KData = squeeze(sampledKData(:,:,xdim,:,:));
                    
                    for ItrIndex = 1 : maxItr % MaxItrNum   
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
                        
                        
                        % Code added on 2021.11.28 to study the evolution of MeanSquareError and kernelDiff ---
                        meanSquareErrorArray(xdim, ItrIndex) = mean( abs( ErrData(:) ).^2 );
                        diffCurArray(xdim, ItrIndex) = norm(temp(:));
                        nItrArray(xdim) = ItrIndex;
                        % -------------------------------------------------------------------------------------
                    
                        if(ItrIndex >2)
                            if(MSEVec(ItrIndex-1)< MSEVec(ItrIndex)|| DiffCur<2.8e-15)
                                aa=MSEVec(ItrIndex-1);           
                                break;
                            end
                        end
                    end%%%  End MaxItrNum

                    ItrIndex
                    PrjTime=toc(Timer2);
                    ItrTime=PrjTime;
                    reconImg(:,:,:,:,xdim) = abs(CurData);
                    end

                    totalRecoTime = toc(Timer3)
                    
                    disp('  ')
                    disp('KLR CS reconstruction of diff. weighted volumes done!!!')
                    disp('  ')

                    reconImg = permute(reconImg,[5,1,2,3,4]);
                    %recon = permute(squeeze(sos(NDWI,4)),[3,1,2,4,5]);
                    %recon(:,:,:,4:(ImgSize(5)+3)) = squeeze(sos(abs(reconImg),4));

                    fprintf("\n Size NDWI:")
                    size(NDWI)
                    fprintf("\n Size reconImg:")
                    size(reconImg)
                    
                    klrCSImageReconstructed = zeros(size(undersampledKSpace,1),size(undersampledKSpace,2),size(undersampledKSpace,3),size(undersampledKSpace,5));
                    
                    klrCSImageReconstructed(:,:,:,1:nB0s) = permute(squeeze(sos(NDWI,4)),[3,1,2,4,5]);
                    klrCSImageReconstructed(:,:,:,(nB0s+1):end) = squeeze(sos(abs(reconImg),4));

                    %[csReconstructedkSpace, csReconstructedImage] =  MultiChannelCSReconstruction (undersampledKSpace, '3D-CS', 1, 1, 1, showResults, whichSliceX, whichSliceY, whichSliceZ);
                    if (correctImage)
                        fprintf("\nCorrecting offset of reconstructed image...\n")
                        % Correcting image with FFTshift -------------------------------------------
                        if (fftShiftCorrection(1))
                            klrCSImageReconstructed = fftshift(klrCSImageReconstructed,1);
                        end
                        if (fftShiftCorrection(2))
                            klrCSImageReconstructed = fftshift(klrCSImageReconstructed,2);
                        end
                        if (fftShiftCorrection(3))
                            klrCSImageReconstructed = fftshift(klrCSImageReconstructed,3);
                        end
                        % --------------------------------------------------------------------------

                        % Correct offset --------------------------------------------------------------------------
                        klrCSImageReconstructed = circshift(klrCSImageReconstructed, offsetCorrection);
                        % -----------------------------------------------------------------------------------------

                    end
                    
                    if (showResults)    
                        figure;
                        whichVolume = 1;
                        subplot(4,3,1)
                        imshow(abs(squeeze(klrCSImageReconstructed(whichSliceX,:,:,whichVolume))),[]);
                        title({'KLR Compressed Sensing reconstruction',['Volume: ',num2str(whichVolume)],['Y-Z plan - slice ', num2str(whichSliceX)]});

                        subplot(4,3,2)
                        imshow(abs(squeeze(klrCSImageReconstructed(:,whichSliceY,:,whichVolume))),[]);
                        title({'KLR Compressed Sensing reconstruction',['Volume: ',num2str(whichVolume)],['X-Z plan - slice ', num2str(whichSliceY)]});

                        subplot(4,3,3)               
                        imshow(abs(squeeze(klrCSImageReconstructed(:,:,whichSliceZ,whichVolume))),[]);
                        title({'KLR Compressed Sensing reconstruction',['Volume: ',num2str(whichVolume)],['X-Y plan - slice ', num2str(whichSliceZ)]});
                        
                        figure;
                        whichVolume = nB0s;
                        subplot(4,3,4)
                        imshow(abs(squeeze(klrCSImageReconstructed(whichSliceX,:,:,whichVolume))),[]);
                        title({'KLR Compressed Sensing reconstruction',['Volume: ',num2str(whichVolume)],['Y-Z plan - slice ', num2str(whichSliceX)]});

                        subplot(4,3,5)
                        imshow(abs(squeeze(klrCSImageReconstructed(:,whichSliceY,:,whichVolume))),[]);
                        title({'KLR Compressed Sensing reconstruction',['Volume: ',num2str(whichVolume)],['X-Z plan - slice ', num2str(whichSliceY)]});

                        subplot(4,3,6)               
                        imshow(abs(squeeze(klrCSImageReconstructed(:,:,whichSliceZ,whichVolume))),[]);
                        title({'KLR Compressed Sensing reconstruction',['Volume: ',num2str(whichVolume)],['X-Y plan - slice ', num2str(whichSliceZ)]});
                        
                        figure;
                        whichVolume = nB0s + 1;
                        subplot(4,3,7)
                        imshow(abs(squeeze(klrCSImageReconstructed(whichSliceX,:,:,whichVolume))),[]);
                        title({'KLR Compressed Sensing reconstruction',['Volume: ',num2str(whichVolume)],['Y-Z plan - slice ', num2str(whichSliceX)]});

                        subplot(4,3,8)
                        imshow(abs(squeeze(klrCSImageReconstructed(:,whichSliceY,:,whichVolume))),[]);
                        title({'KLR Compressed Sensing reconstruction',['Volume: ',num2str(whichVolume)],['X-Z plan - slice ', num2str(whichSliceY)]});

                        subplot(4,3,9)               
                        imshow(abs(squeeze(klrCSImageReconstructed(:,:,whichSliceZ,whichVolume))),[]);
                        title({'KLR Compressed Sensing reconstruction',['Volume: ',num2str(whichVolume)],['X-Y plan - slice ', num2str(whichSliceZ)]});
                        
                        figure;
                        whichVolume = size(klrCSImageReconstructed,4);
                        subplot(4,3,10)
                        imshow(abs(squeeze(klrCSImageReconstructed(whichSliceX,:,:,whichVolume))),[]);
                        title({'KLR Compressed Sensing reconstruction',['Volume: ',num2str(whichVolume)],['Y-Z plan - slice ', num2str(whichSliceX)]});

                        subplot(4,3,11)
                        imshow(abs(squeeze(klrCSImageReconstructed(:,whichSliceY,:,whichVolume))),[]);
                        title({'KLR Compressed Sensing reconstruction',['Volume: ',num2str(whichVolume)],['X-Z plan - slice ', num2str(whichSliceY)]});

                        subplot(4,3,12)               
                        imshow(abs(squeeze(klrCSImageReconstructed(:,:,whichSliceZ,whichVolume))),[]);
                        title({'KLR Compressed Sensing reconstruction',['Volume: ',num2str(whichVolume)],['X-Y plan - slice ', num2str(whichSliceZ)]});
                    end


                    % Save .nii image in a specific folder ------------------------------------------------------------------------------------                
                    pathAndfilename = strcat(savePathAndFileName,'_KLRCS_lambReg-',num2str(lambdaReg),'_ALMReg-',num2str(ALMReg),'_maxItr-',num2str(maxItr),'.nii');
                    niftiwrite((single(abs(klrCSImageReconstructed))),pathAndfilename);
                    fprintf("\nReconstructed image saved as a NIfTI file.\n")
                    % -------------------------------------------------------------------------------------------------------------------------


                    % Pause step 3 ---
                    pause(3)
                    % ----------------

                    fileID = fopen(pathTxtDoc,'a');
                    fprintf(fileID,'%s\r\n',pathAndfilename);
                    fclose(fileID);
                    fprintf("\nPath of reconstructed image added in a txt file.\n")
                    

                    
                end
            end
        end   
        
    end
end