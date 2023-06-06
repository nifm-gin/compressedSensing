%% Compare Bruker acquisitions - Fully sampled k-space vs CS




% -------------------------------------------------------------------------------------------------------------------------------
% -------------------------------------------------------------------------------------------------------------------------------
% -------------------------------------------------------------------------------------------------------------------------------

%% Add BART v0.8.00 path
% set Matlab path and TOOLBOX_PATH environment variable
%bartPath = '/opt/bart-0.6.00';
bartPath = '~/Apps/Bart/bart-0.8.00';
addpath(fullfile(bartPath, 'matlab'));
setenv('TOOLBOX_PATH', bartPath);
bart('version')

addpath(pwd,'diffusionKPCA/KernelLib' );
addpath(pwd, 'diffusionKPCA/ReconLib' );
addpath(pwd, 'diffusionKPCA/PreImLib' );
addpath(pwd, 'diffusionKPCA/FFTLib' );
addpath(pwd, 'diffusionKPCA' );


%%
% ====================================================================================================
% ====================================================================================================
% ====================================================================================================
% MAP6 M10, M17, M15, M16, M19, M18 - Undersampled reconstruction - AF=2
% ====================================================================================================
close all
clear
clc

dataInfo.type = "ParavisionRawFID";
dataInfo.processingType = "UndersampledDataReconstruction";


dataInfo.path = [...
"/data_network/summer_projects/alvesrod/Current/2021_CompresSense_DA_HM_JCD_EB/Datasets/Acquisitions/2021.08.27_MAP6-WT-Het_StudyForPaper/20210827_150342_Mouse_ex_vivo_MAP6_M10_2021_08_27_Mouse_ex__1_1/17/csFid_ch-",...
"/data_network/summer_projects/alvesrod/Current/2021_CompresSense_DA_HM_JCD_EB/Datasets/Acquisitions/2021.08.27_MAP6-WT-Het_StudyForPaper/20210829_082311_Mouse_ex_vivo_MAP6_M17_2021_08_29_Mouse_ex__1_1/10/csFid_ch-",...
"/data_network/summer_projects/alvesrod/Current/2021_CompresSense_DA_HM_JCD_EB/Datasets/Acquisitions/2021.08.27_MAP6-WT-Het_StudyForPaper/20210830_090019_Mouse_ex_vivo_MAP6_M15_2021_08_30_Mouse_ex__1_1/10/csFid_ch-",...
"/data_network/summer_projects/alvesrod/Current/2021_CompresSense_DA_HM_JCD_EB/Datasets/Acquisitions/2021.08.27_MAP6-WT-Het_StudyForPaper/20210904_141307_Mouse_ex_vivo_MAP6_M16_2021_09_04_Mouse_ex__1_1/11/csFid_ch-",...
"/data_network/summer_projects/alvesrod/Current/2021_CompresSense_DA_HM_JCD_EB/Datasets/Acquisitions/2021.08.27_MAP6-WT-Het_StudyForPaper/20210905_151247_Mouse_ex_vivo_MAP6_M19_2021_09_06_Mouse_ex__1_1/10/csFid_ch-",...
"/data_network/summer_projects/alvesrod/Current/2021_CompresSense_DA_HM_JCD_EB/Datasets/Acquisitions/2021.08.27_MAP6-WT-Het_StudyForPaper/20210907_162349_Mouse_ex_vivo_MAP6_M18_2021_09_07_Mouse_ex__1_2/11/csFid_ch-"];

dataInfo.acquisitionParameters = [...
"/data_network/summer_projects/alvesrod/Current/2021_CompresSense_DA_HM_JCD_EB/Datasets/Acquisitions/2021.08.27_MAP6-WT-Het_StudyForPaper/20210827_150342_Mouse_ex_vivo_MAP6_M10_2021_08_27_Mouse_ex__1_1/17/Acquisition-parameters.txt",...
"/data_network/summer_projects/alvesrod/Current/2021_CompresSense_DA_HM_JCD_EB/Datasets/Acquisitions/2021.08.27_MAP6-WT-Het_StudyForPaper/20210829_082311_Mouse_ex_vivo_MAP6_M17_2021_08_29_Mouse_ex__1_1/10/Acquisition-parameters.txt",...
"/data_network/summer_projects/alvesrod/Current/2021_CompresSense_DA_HM_JCD_EB/Datasets/Acquisitions/2021.08.27_MAP6-WT-Het_StudyForPaper/20210830_090019_Mouse_ex_vivo_MAP6_M15_2021_08_30_Mouse_ex__1_1/10/Acquisition-parameters.txt",...
"/data_network/summer_projects/alvesrod/Current/2021_CompresSense_DA_HM_JCD_EB/Datasets/Acquisitions/2021.08.27_MAP6-WT-Het_StudyForPaper/20210904_141307_Mouse_ex_vivo_MAP6_M16_2021_09_04_Mouse_ex__1_1/11/Acquisition-parameters.txt",...
"/data_network/summer_projects/alvesrod/Current/2021_CompresSense_DA_HM_JCD_EB/Datasets/Acquisitions/2021.08.27_MAP6-WT-Het_StudyForPaper/20210905_151247_Mouse_ex_vivo_MAP6_M19_2021_09_06_Mouse_ex__1_1/10/Acquisition-parameters.txt",...
"/data_network/summer_projects/alvesrod/Current/2021_CompresSense_DA_HM_JCD_EB/Datasets/Acquisitions/2021.08.27_MAP6-WT-Het_StudyForPaper/20210907_162349_Mouse_ex_vivo_MAP6_M18_2021_09_07_Mouse_ex__1_2/11/Acquisition-parameters.txt"];

dataInfo.undersamplingPatternParameters = [...
"/data_network/summer_projects/alvesrod/Current/2021_CompresSense_DA_HM_JCD_EB/Datasets/Acquisitions/2021.08.27_MAP6-WT-Het_StudyForPaper/20210827_150342_Mouse_ex_vivo_MAP6_M10_2021_08_27_Mouse_ex__1_1/17/20210814_225850_UndersamplingPattern_180x90x110_AF-2_overallAF-1.9995_DWIAF-2.2222_nVol-33_nB0-3_fsCenterArea-706.8583_MonteCarlo-VariableDensity_Multi3D_FullySampledB0_BinAcqGuideArray.txt",...
"/data_network/summer_projects/alvesrod/Current/2021_CompresSense_DA_HM_JCD_EB/Datasets/Acquisitions/2021.08.27_MAP6-WT-Het_StudyForPaper/20210829_082311_Mouse_ex_vivo_MAP6_M17_2021_08_29_Mouse_ex__1_1/10/20210814_225850_UndersamplingPattern_180x90x110_AF-2_overallAF-1.9995_DWIAF-2.2222_nVol-33_nB0-3_fsCenterArea-706.8583_MonteCarlo-VariableDensity_Multi3D_FullySampledB0_BinAcqGuideArray.txt",...
"/data_network/summer_projects/alvesrod/Current/2021_CompresSense_DA_HM_JCD_EB/Datasets/Acquisitions/2021.08.27_MAP6-WT-Het_StudyForPaper/20210830_090019_Mouse_ex_vivo_MAP6_M15_2021_08_30_Mouse_ex__1_1/10/20210814_225850_UndersamplingPattern_180x90x110_AF-2_overallAF-1.9995_DWIAF-2.2222_nVol-33_nB0-3_fsCenterArea-706.8583_MonteCarlo-VariableDensity_Multi3D_FullySampledB0_BinAcqGuideArray.txt",...
"/data_network/summer_projects/alvesrod/Current/2021_CompresSense_DA_HM_JCD_EB/Datasets/Acquisitions/2021.08.27_MAP6-WT-Het_StudyForPaper/20210904_141307_Mouse_ex_vivo_MAP6_M16_2021_09_04_Mouse_ex__1_1/11/20210814_225850_UndersamplingPattern_180x90x110_AF-2_overallAF-1.9995_DWIAF-2.2222_nVol-33_nB0-3_fsCenterArea-706.8583_MonteCarlo-VariableDensity_Multi3D_FullySampledB0_BinAcqGuideArray.txt",...
"/data_network/summer_projects/alvesrod/Current/2021_CompresSense_DA_HM_JCD_EB/Datasets/Acquisitions/2021.08.27_MAP6-WT-Het_StudyForPaper/20210905_151247_Mouse_ex_vivo_MAP6_M19_2021_09_06_Mouse_ex__1_1/10/20210814_225850_UndersamplingPattern_180x90x110_AF-2_overallAF-1.9995_DWIAF-2.2222_nVol-33_nB0-3_fsCenterArea-706.8583_MonteCarlo-VariableDensity_Multi3D_FullySampledB0_BinAcqGuideArray.txt",...
"/data_network/summer_projects/alvesrod/Current/2021_CompresSense_DA_HM_JCD_EB/Datasets/Acquisitions/2021.08.27_MAP6-WT-Het_StudyForPaper/20210907_162349_Mouse_ex_vivo_MAP6_M18_2021_09_07_Mouse_ex__1_2/11/20210814_225850_UndersamplingPattern_180x90x110_AF-2_overallAF-1.9995_DWIAF-2.2222_nVol-33_nB0-3_fsCenterArea-706.8583_MonteCarlo-VariableDensity_Multi3D_FullySampledB0_BinAcqGuideArray.txt"];

dataInfo.saveDataPath = "/data_local/data_hdd/alvesrod/2022.01.14_MAP6-WT-Het_StudyForCSPaper/Recos_RealAcquisitions/AF-2/";

dataInfo.saveDataName = [...
"Acq-2021.08.27_M10-WT_3D-SE-DTI_Mouse-brain_Type-1-seq",...
"Acq-2021.08.29_M17-MAP6-Het_3D-SE-DTI_Mouse-brain_Type-1-seq",...
"Acq-2021.08.30_M15-WT_3D-SE-DTI_Mouse-brain_Type-1-seq",...
"Acq-2021.09.04_M16-WT_3D-SE-DTI_Mouse-brain_Type-1-seq",...
"Acq-2021.09.05_M19-MAP6-Het_3D-SE-DTI_Mouse-brain_Type-1-seq",...
"Acq-2021.09.07_M18-MAP6-Het_3D-SE-DTI_Mouse-brain_Type-1-seq"];


% Parameters of the undersampling pattern must match the ones used for the acquisition ---
dataInfo.B0SamplingStrategy = ["FullySampledB0"];
dataInfo.AF = [2];
dataInfo.undersamplingPatternType = ["MonteCarlo-VariableDensity"];
dataInfo.undersamplingPatternStrategy = ["Multi3D"];
dataInfo.undersamplingPatternFullySampledRegionArea = [(pi*15^2)];

dataInfo.actual4DAF = 1.9995;
dataInfo.correctedAF = 2.2222;
% ----------------------------------------------------------------------------------------


dataInfo.reconstructionMethod = ["ZeroPadding","ConventionalCS","KLRCS"];

dataInfo.CSParamLamba1 = [0.005];
dataInfo.CSParamLamba2 = [0.002];
dataInfo.CSParamNItrMax = [200];

dataInfo.KLRCSParamALMReg = [0.1];
dataInfo.KLRCSParamLambdaReg = [0.1];
dataInfo.KLRCSParamNItrMax = [500];

dataInfo.xCenterCubeForB0Norm = [63,70,69,71,58,64];
dataInfo.yCenterCubeForB0Norm = [74,75,75,73,71,70];
dataInfo.zCenterCubeForB0Norm = [69,69,73,71,69,66];
dataInfo.cubeSizeForB0Norm = [5,5,5,5,5,5];


dataInfo.imageCorrection = 1;
dataInfo.offsetCorrection = [0,-25,7,0,0];

dataInfo.fftshiftCorrection = [1,0,0,0,0];
dataInfo.correctionBeforeReconstruction = 0;
dataInfo.showResults = 1;
dataInfo.slicesToShow = [90,45,55,1,1];
dataInfo.sliceRotation = [180,0,90];
dataInfo.saveFileType = ["NifTI"];

% ====================================================================================================
% ====================================================================================================
% ====================================================================================================



% ####################################################################################################
% ####################################################################################################
% ####################################################################################################





%%


close all
clc

offsetCorrectionOriginal = dataInfo.offsetCorrection;
fftshiftCorrectionOriginal = dataInfo.fftshiftCorrection;



fprintf("\n======================================================================\nVerification if all paths exist:\n")

for nSubject=1:size(dataInfo.path,2)
    
    if (exist(strcat(dataInfo.path(:,nSubject),"0")) && exist((dataInfo.acquisitionParameters(:,nSubject))) && exist((dataInfo.undersamplingPatternParameters(:,nSubject))) && exist(dataInfo.saveDataPath))
        fprintf("All paths are right for subject %d:\n",nSubject)
        
    else
        fprintf("\nERROR: Either the paths for the subject %d or the one for saving data do not exist!!!\nPlease verify the following paths:\n",nSubject)
        dataInfo.path(:,nSubject)
        dataInfo.acquisitionParameters(:,nSubject)
        dataInfo.undersamplingPatternParameters(:,nSubject)
        dataInfo.saveDataPath
        return
    end
    
    
end
fprintf("\nAll paths actually exist\n======================================================================\n")
    

fprintf("\nAttention: The matrix size of all subjects must be the same, given that the same undersampling patterns will be used at all!!!\n")


if dataInfo.processingType=="SimulationCSFromFullySampledData"
                
    
    fprintf('\nSimulation of CS acquisition by undersampling a fully sampled k-space.\n')
    
    
    if (dataInfo.loadReadyUndersamplinPatterns)
        fileID = fopen(dataInfo.pathTxtDocContainingUndersamplinPatternPaths,'r');
        pathTxtDoc = textscan(fileID,'%s','delimiter','\n');
        fclose(fileID);
        
        %nUndersamplingPatterns = str2num(pathTxtDoc{1,1}{1});
        nUndersamplingPatterns = size(pathTxtDoc{1,1},1);
        
                
            

                
        for nSubject=1:size(dataInfo.path,2)
            
            dataInfo.offsetCorrection = offsetCorrectionOriginal;
            dataInfo.fftshiftCorrection = fftshiftCorrectionOriginal;
                                                                    
            disp(' ')
            dataPath = dataInfo.path(nSubject);
            acquisitionFilePath = dataInfo.acquisitionParameters(nSubject);

            disp('============================================================')
            % Defining path to save fully sampled data -------------------------------------------       
            timeNow = datevec(now);
            savePath = (strcat(dataInfo.saveDataPath,num2str(timeNow(1)),num2str(timeNow(2),'%.2d'),num2str(timeNow(3),'%.2d'),'_',num2str(timeNow(4),'%.2d'),num2str(timeNow(5),'%.2d'),num2str(timeNow(6),'%02.0f'),'_',dataInfo.saveDataName(nSubject)));
            mkdir (savePath);

            savePathAndFileName = strcat(savePath,'/',dataInfo.saveDataName(nSubject));
            % ------------------------------------------------------------------------------------                       

            disp(strcat('Subject No: ',num2str(nSubject)))
            disp(strcat('Subject: ',dataPath))


            % Read acquisition file in order to reconstruct k-space correctly -------------------------------------------------------------------------------------------------------------------------
            [acquisitionTypeID, origReadoutDim, origPhase1Dim, origPhase2Dim, nEchos, finalReadoutDim, finalPhase1Dim, finalPhase2Dim, nChannels, nB0s, nDWIs, nVolumes] = readAcquisitionFile(acquisitionFilePath);
            % -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------



            % Identify acquisition type --------------------------
            if (acquisitionTypeID==1)
                fsAcquisitionType = "3D-FLASH";

            elseif (acquisitionTypeID==2)
                fsAcquisitionType = "3D-SE-DTI_DWLoopFirst";
                
            elseif (acquisitionTypeID==3)
                fsAcquisitionType = "3D-SE-DTI_DWLoopLast";

            else
                disp(' ')
                disp('Acquisition type out of allowed cases!!!');
                disp(' ')
                %return

            end
            % ---------------------------------------------------
            
            
            % Create txt document containing all the reconstruction for future data processing ---
            pathListOfRecoTxtDoc = strcat(savePathAndFileName,"_listOfReconstructions.txt");
            fileID = fopen(pathListOfRecoTxtDoc,'w');            
            fclose(fileID);
            % ------------------------------------------------------------------------------------


            % Retrieving an organized k-space -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
            if (dataInfo.type == "ParavisionRawFID")
                disp(' ')
                disp('Reconstruction Paravision raw FID files...')

                [fullySampledKSpace, ~, ~] = fullySampledDataReconstructionFromMultiChannelRawFID(dataPath,origReadoutDim, origPhase1Dim, origPhase2Dim, nEchos, finalReadoutDim, finalPhase1Dim, finalPhase2Dim, nChannels, nVolumes, fsAcquisitionType, dataInfo.imageCorrection, dataInfo.fftshiftCorrection, dataInfo.offsetCorrection, dataInfo.showResults, dataInfo.slicesToShow(1), dataInfo.slicesToShow(2), dataInfo.slicesToShow(3), dataInfo.slicesToShow(5), dataInfo.sliceRotation(1), dataInfo.sliceRotation(2), dataInfo.sliceRotation(3), savePathAndFileName, pathListOfRecoTxtDoc);
                
                
                if(dataInfo.correctionBeforeReconstruction)
                    close all
                    fprintf("\n####################################################\nCorrection before reconstruction\n##################################################");
                    %ref_im = ifftshift(ifft(ifftshift(ifft(ifftshift(ifft(fullySampledKSpace,[],1),1),[],2),2),[],3),3);                                                
                    ref_im = (ifft((ifft((ifft(fullySampledKSpace,[],1)),[],2)),[],3));                                                

                    ref_im = circshift(ref_im,dataInfo.offsetCorrection);

                    if (dataInfo.fftshiftCorrection(1))
                        ref_im = fftshift(ref_im,1);
                    end
                    if (dataInfo.fftshiftCorrection(2))
                        ref_im = fftshift(ref_im,2);
                    end
                    if (dataInfo.fftshiftCorrection(3))
                        ref_im = fftshift(ref_im,3);
                    end
                    if (dataInfo.fftshiftCorrection(4))
                        ref_im = fftshift(ref_im,4);
                    end
                    if (dataInfo.fftshiftCorrection(5))
                        ref_im = fftshift(ref_im,5);
                    end


                    fprintf("\n Showing images with k-space corrected before reconstruction...\n");
                    showImage(ref_im, dataInfo.showResults, dataInfo.slicesToShow(1), dataInfo.slicesToShow(2), dataInfo.slicesToShow(3), dataInfo.slicesToShow(5), dataInfo.sliceRotation(1), dataInfo.sliceRotation(2), dataInfo.sliceRotation(3));                                
                    
                    pause(5)
                    
                    ImgSize = size(fullySampledKSpace);
                    fullySampledKSpaceCorrected = FFT( ref_im, ones( ImgSize )); 
                    showImage(fullySampledKSpaceCorrected, dataInfo.showResults, dataInfo.slicesToShow(1), dataInfo.slicesToShow(2), dataInfo.slicesToShow(3), dataInfo.slicesToShow(5), 0, 0, 0);
                    
                    dataInfo.offsetCorrection = [0,0,0,0,0];
                    dataInfo.fftshiftCorrection = [0,0,0,0,0];
                    
                    fullySampledKSpace = fullySampledKSpaceCorrected;
                    
                    ref_im = ifftshift(ifft(ifftshift(ifft(ifftshift(ifft(fullySampledKSpace,[],1),1),[],2),2),[],3),3);                                                                                        
                    showImage(ref_im, dataInfo.showResults, dataInfo.slicesToShow(1), dataInfo.slicesToShow(2), dataInfo.slicesToShow(3), dataInfo.slicesToShow(5), dataInfo.sliceRotation(1), dataInfo.sliceRotation(2), dataInfo.sliceRotation(3));                                
                    
                    %return
                    
                    % Enlever 
                    % -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
                    %fullySampledKSpaceCorrected = fftshift(fft(fftshift(fft(fftshift(fft(ref_im,[],1),1),[],2),2),[],3),3);                
                    %{
                    fullySampledKSpaceCorrected1 = (fft((fft((fft(ref_im,[],1)),[],2)),[],3));                
                    fullySampledKSpaceCorrected2 = fftshift(fft(fftshift(fft(fftshift(fft(ref_im,[],1),1),[],2),2),[],3),3);                

                    fullySampledKSpaceCorrectedView1 = squeeze(bart('rss 8', fullySampledKSpaceCorrected1));
                    fullySampledKSpaceCorrectedView2 = squeeze(bart('rss 8', fullySampledKSpaceCorrected2));

                    showImage(fullySampledKSpaceCorrectedView1, dataInfo.showResults, dataInfo.slicesToShow(1), dataInfo.slicesToShow(2), dataInfo.slicesToShow(3), dataInfo.slicesToShow(5), 0, 0, 0);                                
                    showImage(fullySampledKSpaceCorrectedView2, dataInfo.showResults, dataInfo.slicesToShow(1), dataInfo.slicesToShow(2), dataInfo.slicesToShow(3), dataInfo.slicesToShow(5), 0, 0, 0);                                

                    %fullySampledImageCorrected = bart('fft -i 7', fullySampledKSpaceCorrected);
                    fullySampledImageCorrected1 = ifftshift(ifft(ifftshift(ifft(ifftshift(ifft(fullySampledKSpaceCorrected1,[],1),1),[],2),2),[],3),3);

                    fullySampledImageCorrected2 = (ifft((ifft((ifft(fullySampledKSpaceCorrected2,[],1)),[],2)),[],3));
                    %fullySampledImageCorrected = (ifft((ifft((ifft(fullySampledKSpaceCorrected,[],1)),[],2)),[],3));                
                    fullySampledImageCorrectedView1 = bart('rss 8', fullySampledImageCorrected1);
                    fullySampledImageCorrectedView2 = bart('rss 8', fullySampledImageCorrected2);

                    fullySampledImageCorrectedView1 = squeeze(fullySampledImageCorrectedView1);
                    fullySampledImageCorrectedView2 = squeeze(fullySampledImageCorrectedView2);
                    
                    showImage(fullySampledImageCorrectedView1, dataInfo.showResults, dataInfo.slicesToShow(1), dataInfo.slicesToShow(2), dataInfo.slicesToShow(3), dataInfo.slicesToShow(5), dataInfo.sliceRotation(1), dataInfo.sliceRotation(2), dataInfo.sliceRotation(3));                                
                    showImage(fullySampledImageCorrectedView2, dataInfo.showResults, dataInfo.slicesToShow(1), dataInfo.slicesToShow(2), dataInfo.slicesToShow(3), dataInfo.slicesToShow(5), dataInfo.sliceRotation(1), dataInfo.sliceRotation(2), dataInfo.sliceRotation(3));                                
                    %}
                    % -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
                
                end
                
                disp('Fully sampled data reconstruction and saving process done!!!')
                disp('============================================================')


            elseif dataInfo.type == "k-spaceFromMatlabWorkspace"
                % To be coded !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                disp(' ')
                disp('Reading k-space from .mat file...')
                disp('K-space reading done!!!')
                disp(' ')
                % To be coded !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            else
                disp(' ')
                disp('Data type out of allowed cases!!!');
                disp(' ')
                %return                            
            end
            % -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------                                    
            
            
            
            %for nMask=2:(nUndersamplingPatterns+1)
            %parfor nMask=1:nUndersamplingPatterns
            for nMask=1:nUndersamplingPatterns
                fprintf("\n#############################################################\nProcessing undersampling pattern N: %d\n#############################################################\n",nMask)
                undersamplingPatternData=load(pathTxtDoc{1,1}{nMask});
                AF = undersamplingPatternData.AF;
                B0SamplingStrategy = undersamplingPatternData.B0SamplingStrategy;
                actual4DAF = undersamplingPatternData.actual4DAF;
                correctedAF = undersamplingPatternData.correctedAF;
                maskDiff = undersamplingPatternData.maskDiff;
                nB0s = undersamplingPatternData.nB0s;
                undersamplingPattern = undersamplingPatternData.undersamplingPattern5D; % Temp code - The correct is undersamplingPattern = undersamplingPatternData.undersamplingPattern
                undersamplingPatternFullySampledRegionArea = undersamplingPatternData.undersamplingPatternFullySampledRegionArea;
                undersamplingPatternStrategy = undersamplingPatternData.undersamplingPatternStrategy;
                undersamplingPatternType = undersamplingPatternData.undersamplingPatternType;
                
                % Undersampling k-space ----------------------------------------------------------------------------------
                fprintf(strcat('Undersampling k-space using ',undersamplingPatternType,' undersampling pattern and ',undersamplingPatternStrategy,' strategy...\n\n'))

            

                undersampledKSpace = fullySampledKSpace.*undersamplingPattern;
                
                % Enlever
                % -----------------------------------------------------------------------
                %undersampledKSpace1 = fullySampledKSpaceCorrected1.*undersamplingPattern;
                %undersampledKSpace2 = fullySampledKSpaceCorrected2.*undersamplingPattern;
                % -----------------------------------------------------------------------


                fprintf('Undersampling done!!!\n\n')
                % --------------------------------------------------------------------------------------------------------
                    
                
                for recoMeth=1:size(dataInfo.reconstructionMethod,2)
                    
                    reconstructionMethod = dataInfo.reconstructionMethod(recoMeth);
                    

                    % Reconstruction image from undersampled k-space ----------------------------------------------------------------------------------
                    if (reconstructionMethod=="ConventionalCS")
                        fprintf(strcat('Reconsctruction of undersampled k-space via ',reconstructionMethod,', using a',undersamplingPatternType,' undersampling pattern and ',undersamplingPatternStrategy,' strategy...\n\n'))
                    
                        savePathAndFileNameAfterReco = strcat(savePathAndFileName,"_AF-",num2str(AF),"_OverallAF-",num2str(actual4DAF),"_DWIAF-",num2str(correctedAF),"_",undersamplingPatternType,"_",undersamplingPatternStrategy,"_FSArea-",num2str(undersamplingPatternFullySampledRegionArea),"_",B0SamplingStrategy);

                        [conventionalCSImageReconstructed, coilSensitivityMapsViewCorrected] = conventionalCSReconstruction_v2 (undersampledKSpace, dataInfo.CSParamLamba1, dataInfo.CSParamLamba2, dataInfo.CSParamNItrMax, dataInfo.imageCorrection, dataInfo.fftshiftCorrection, dataInfo.offsetCorrection, B0SamplingStrategy, nB0s, dataInfo.xCenterCubeForB0Norm(nSubject), dataInfo.yCenterCubeForB0Norm(nSubject), dataInfo.zCenterCubeForB0Norm(nSubject), dataInfo.cubeSizeForB0Norm(nSubject), savePathAndFileNameAfterReco, pathListOfRecoTxtDoc, dataInfo.showResults, dataInfo.slicesToShow(1), dataInfo.slicesToShow(2), dataInfo.slicesToShow(3));

                        % Enlever
                        % --------------------------------------------------------------------------------------------------------
                        %[conventionalCSImageReconstructed1, coilSensitivityMapsViewCorrected1] = conventionalCSReconstruction_v2 (undersampledKSpace1, dataInfo.CSParamLamba1(1), dataInfo.CSParamLamba2(1), dataInfo.CSParamNItrMax(1), 0, dataInfo.fftshiftCorrection, dataInfo.offsetCorrection, B0SamplingStrategy, nB0s, dataInfo.xCenterCubeForB0Norm(nSubject), dataInfo.yCenterCubeForB0Norm(nSubject), dataInfo.zCenterCubeForB0Norm(nSubject), dataInfo.cubeSizeForB0Norm(nSubject), strcat(savePathAndFileNameAfterReco,"_1"), pathListOfRecoTxtDoc, dataInfo.showResults, dataInfo.slicesToShow(1), dataInfo.slicesToShow(2), dataInfo.slicesToShow(3));
                        %[conventionalCSImageReconstructed2, coilSensitivityMapsViewCorrected2] = conventionalCSReconstruction_v2 (undersampledKSpace2, dataInfo.CSParamLamba1(1), dataInfo.CSParamLamba2(1), dataInfo.CSParamNItrMax(1), 0, dataInfo.fftshiftCorrection, dataInfo.offsetCorrection, B0SamplingStrategy, nB0s, dataInfo.xCenterCubeForB0Norm(nSubject), dataInfo.yCenterCubeForB0Norm(nSubject), dataInfo.zCenterCubeForB0Norm(nSubject), dataInfo.cubeSizeForB0Norm(nSubject), strcat(savePathAndFileNameAfterReco,"_2"), pathListOfRecoTxtDoc, dataInfo.showResults, dataInfo.slicesToShow(1), dataInfo.slicesToShow(2), dataInfo.slicesToShow(3));
                        
                        %save('conventionalCSImageReconstructed1.mat','conventionalCSImageReconstructed1');
                        %save('conventionalCSImageReconstructed2.mat','conventionalCSImageReconstructed2');
                        % ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
                        fprintf('Conventional CS reconstruction done!!!\n\n')
                        close all
                        
                    elseif (reconstructionMethod=="KLRCS")
                        fprintf(strcat('Reconsctruction of undersampled k-space via ',reconstructionMethod,', using a',undersamplingPatternType,' undersampling pattern and ',undersamplingPatternStrategy,' strategy...\n\n'))                    
                        savePathAndFileNameAfterReco = strcat(savePathAndFileName,"_AF-",num2str(AF),"_OverallAF-",num2str(actual4DAF),"_DWIAF-",num2str(correctedAF),"_",undersamplingPatternType,"_",undersamplingPatternStrategy,"_FSArea-",num2str(undersamplingPatternFullySampledRegionArea),"_",B0SamplingStrategy);                                       
                        
                        %undersampledKSpace = fullySampledKSpace.*undersamplingPattern;
                        
                        [klrCSImageReconstructed] = KLRCSReconstruction_v2 (fullySampledKSpace, undersamplingPattern, dataInfo.KLRCSParamLambdaReg, dataInfo.KLRCSParamALMReg, dataInfo.KLRCSParamNItrMax, dataInfo.imageCorrection, dataInfo.correctionBeforeReconstruction, dataInfo.fftshiftCorrection, dataInfo.offsetCorrection, B0SamplingStrategy, nB0s, dataInfo.xCenterCubeForB0Norm(nSubject), dataInfo.yCenterCubeForB0Norm(nSubject), dataInfo.zCenterCubeForB0Norm(nSubject), dataInfo.cubeSizeForB0Norm(nSubject), savePathAndFileNameAfterReco, pathListOfRecoTxtDoc, dataInfo.showResults, dataInfo.slicesToShow(1), dataInfo.slicesToShow(2), dataInfo.slicesToShow(3));
                        
                        fprintf('KLR CS reconstruction done!!!\n\n')
                        
                    elseif (reconstructionMethod=="ZeroPadding")
                        fprintf(strcat('Reconsctruction of undersampled k-space via ',reconstructionMethod,', using a',undersamplingPatternType,' undersampling pattern and ',undersamplingPatternStrategy,' strategy...\n\n'))
                    
                        savePathAndFileNameAfterReco = strcat(savePathAndFileName,"_AF-",num2str(AF),"_OverallAF-",num2str(actual4DAF),"_DWIAF-",num2str(correctedAF),"_",undersamplingPatternType,"_",undersamplingPatternStrategy,"_FSArea-",num2str(undersamplingPatternFullySampledRegionArea),"_",B0SamplingStrategy);
                        
                        [zpImageReconstructed] = zeroPaddingReconstruction (undersampledKSpace, dataInfo.imageCorrection, dataInfo.fftshiftCorrection, dataInfo.offsetCorrection, dataInfo.showResults, dataInfo.slicesToShow(1), dataInfo.slicesToShow(2), dataInfo.slicesToShow(3), dataInfo.slicesToShow(5), dataInfo.sliceRotation(1), dataInfo.sliceRotation(2), dataInfo.sliceRotation(3), savePathAndFileNameAfterReco, pathListOfRecoTxtDoc);
                        
                        fprintf('Zero padding reconstruction done!!!\n\n')
                    end
                    
                    
                    % --------------------------------------------------------------------------------------------------------
                    
                end
            end
        end
    end
    
    
    
    
    
    
    

elseif dataInfo.processingType=="FullySampledDataReconstruction"    
    fprintf('\nDirect reconstruction of fully sampled data via FT...\n')    
    
    
    for nSubject=1:size(dataInfo.path,2)
            
        %dataInfo.offsetCorrection = offsetCorrectionOriginal;
        %dataInfo.fftshiftCorrection = fftshiftCorrectionOriginal;

        
        dataPath = dataInfo.path(nSubject);
        acquisitionFilePath = dataInfo.acquisitionParameters(nSubject);

        fprintf('============================================================')
        % Defining path to save fully sampled data -------------------------------------------       
        timeNow = datevec(now);
        savePath = (strcat(dataInfo.saveDataPath,num2str(timeNow(1)),num2str(timeNow(2),'%.2d'),num2str(timeNow(3),'%.2d'),'_',num2str(timeNow(4),'%.2d'),num2str(timeNow(5),'%.2d'),num2str(timeNow(6),'%02.0f'),'_',dataInfo.saveDataName(nSubject)));
        mkdir (savePath);

        savePathAndFileName = strcat(savePath,'/',dataInfo.saveDataName(nSubject));
        % ------------------------------------------------------------------------------------                       

        disp(strcat('Subject No: ',num2str(nSubject)))
        disp(strcat('Subject: ',dataPath))


        % Read acquisition file in order to reconstruct k-space correctly -------------------------------------------------------------------------------------------------------------------------
        [acquisitionTypeID, origReadoutDim, origPhase1Dim, origPhase2Dim, nEchos, finalReadoutDim, finalPhase1Dim, finalPhase2Dim, nChannels, nB0s, nDWIs, nVolumes] = readAcquisitionFile(acquisitionFilePath);
        % -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------



        % Identify acquisition type --------------------------
        if (acquisitionTypeID==1)
            fsAcquisitionType = "3D-FLASH";

        elseif (acquisitionTypeID==2)
            fsAcquisitionType = "3D-SE-DTI_DWLoopFirst";

        elseif (acquisitionTypeID==3)
            fsAcquisitionType = "3D-SE-DTI_DWLoopLast";

        else
            disp(' ')
            disp('Acquisition type out of allowed cases!!!');
            disp(' ')
            %return

        end
        % ---------------------------------------------------


        % Create txt document containing all the reconstruction for future data processing ---
        pathListOfRecoTxtDoc = strcat(savePathAndFileName,"_listOfReconstructions.txt");
        fileID = fopen(pathListOfRecoTxtDoc,'w');            
        fclose(fileID);
        % ------------------------------------------------------------------------------------


        % Retrieving an organized k-space -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
        if (dataInfo.type == "ParavisionRawFID")
            disp(' ')
            disp('Reconstruction Paravision raw FID files...')

            [fullySampledKSpace, ~, ~] = fullySampledDataReconstructionFromMultiChannelRawFID(dataPath,origReadoutDim, origPhase1Dim, origPhase2Dim, nEchos, finalReadoutDim, finalPhase1Dim, finalPhase2Dim, nChannels, nVolumes, fsAcquisitionType, dataInfo.imageCorrection, dataInfo.fftshiftCorrection, dataInfo.offsetCorrection, dataInfo.showResults, dataInfo.slicesToShow(1), dataInfo.slicesToShow(2), dataInfo.slicesToShow(3), dataInfo.slicesToShow(5), dataInfo.sliceRotation(1), dataInfo.sliceRotation(2), dataInfo.sliceRotation(3), savePathAndFileName, pathListOfRecoTxtDoc);


            if(dataInfo.correctionBeforeReconstruction)
                close all
                fprintf("\n####################################################\nCorrection before reconstruction\n##################################################");
                %ref_im = ifftshift(ifft(ifftshift(ifft(ifftshift(ifft(fullySampledKSpace,[],1),1),[],2),2),[],3),3);                                                
                ref_im = (ifft((ifft((ifft(fullySampledKSpace,[],1)),[],2)),[],3));                                                

                ref_im = circshift(ref_im,dataInfo.offsetCorrection);

                if (dataInfo.fftshiftCorrection(1))
                    ref_im = fftshift(ref_im,1);
                end
                if (dataInfo.fftshiftCorrection(2))
                    ref_im = fftshift(ref_im,2);
                end
                if (dataInfo.fftshiftCorrection(3))
                    ref_im = fftshift(ref_im,3);
                end
                if (dataInfo.fftshiftCorrection(4))
                    ref_im = fftshift(ref_im,4);
                end
                if (dataInfo.fftshiftCorrection(5))
                    ref_im = fftshift(ref_im,5);
                end


                fprintf("\n Showing images with k-space corrected before reconstruction...\n");
                showImage(ref_im, dataInfo.showResults, dataInfo.slicesToShow(1), dataInfo.slicesToShow(2), dataInfo.slicesToShow(3), dataInfo.slicesToShow(5), dataInfo.sliceRotation(1), dataInfo.sliceRotation(2), dataInfo.sliceRotation(3));                                

                pause(5)

                ImgSize = size(fullySampledKSpace);
                fullySampledKSpaceCorrected = FFT( ref_im, ones( ImgSize )); 
                showImage(fullySampledKSpaceCorrected, dataInfo.showResults, dataInfo.slicesToShow(1), dataInfo.slicesToShow(2), dataInfo.slicesToShow(3), dataInfo.slicesToShow(5), 0, 0, 0);

                dataInfo.offsetCorrection = [0,0,0,0,0];
                dataInfo.fftshiftCorrection = [0,0,0,0,0];

                fullySampledKSpace = fullySampledKSpaceCorrected;

                ref_im = ifftshift(ifft(ifftshift(ifft(ifftshift(ifft(fullySampledKSpace,[],1),1),[],2),2),[],3),3);                                                                                        
                showImage(ref_im, dataInfo.showResults, dataInfo.slicesToShow(1), dataInfo.slicesToShow(2), dataInfo.slicesToShow(3), dataInfo.slicesToShow(5), dataInfo.sliceRotation(1), dataInfo.sliceRotation(2), dataInfo.sliceRotation(3));                                

                
            end

            disp('Fully sampled data reconstruction and saving process done!!!')
            disp('============================================================')
        end
    end
    
elseif dataInfo.processingType=="UndersampledDataReconstruction"
    disp(' ')
    disp('Reconstruction of undersampled data by ZP, conventional CS or Kernel-Low Rank CS')
    disp(' ')
    
    % Attention - Correct it (the right one is the one below with some modifications) ---
    dataInfo.offsetCorrection = offsetCorrectionOriginal;
    dataInfo.fftshiftCorrection = fftshiftCorrectionOriginal;
    % -----------------------------------------------------------------------------------
    %parfor nSubject=1:size(dataInfo.path,2)
    for nSubject=1:size(dataInfo.path,2)
            
        %dataInfo.offsetCorrection = offsetCorrectionOriginal;
        %dataInfo.fftshiftCorrection = fftshiftCorrectionOriginal;

        disp(' ')
        dataPath = dataInfo.path(nSubject);
        acquisitionFilePath = dataInfo.acquisitionParameters(nSubject);
        undersamplingPatternFilePath = dataInfo.undersamplingPatternParameters(nSubject);

        disp('============================================================\n')
        % Defining path to save fully sampled data -------------------------------------------       
        timeNow = datevec(now);
        savePath = (strcat(dataInfo.saveDataPath,num2str(timeNow(1)),num2str(timeNow(2),'%.2d'),num2str(timeNow(3),'%.2d'),'_',num2str(timeNow(4),'%.2d'),num2str(timeNow(5),'%.2d'),num2str(timeNow(6),'%02.0f'),'_',dataInfo.saveDataName(nSubject)));
        mkdir (savePath);

        savePathAndFileName = strcat(savePath,'/',dataInfo.saveDataName(nSubject));
        % ------------------------------------------------------------------------------------                       

        disp(strcat('Subject No: ',num2str(nSubject)))
        disp(strcat('Subject: ',dataPath))


        % Read acquisition file in order to reconstruct k-space correctly -------------------------------------------------------------------------------------------------------------------------
        [acquisitionTypeID, origReadoutDim, origPhase1Dim, origPhase2Dim, nEchos, finalReadoutDim, finalPhase1Dim, finalPhase2Dim, nChannels, nB0s, nDWIs, nVolumes] = readAcquisitionFile(acquisitionFilePath);
        % -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------



        % Identify acquisition type --------------------------
        if (acquisitionTypeID==1)
            acquisitionType = "3D-FLASH";

        elseif (acquisitionTypeID==2)
            acquisitionType = "CS-3D-SE-DTI_DWLoopFirst";

        elseif (acquisitionTypeID==3)
            acquisitionType = "CS-3D-SE-DTI_DWLoopLast";
            
        elseif (acquisitionTypeID==4)
            acquisitionType = "CS-3D-SE-DTI_DWLoopFirst_BinaryGuideArray";
            

        else
            disp(' ')
            disp('Acquisition type out of allowed cases!!!');
            disp(' ')
            %return

        end
        % ---------------------------------------------------


        % Create txt document containing all the reconstruction for future data processing ---
        pathListOfRecoTxtDoc = strcat(savePathAndFileName,"_listOfReconstructions.txt");
        fileID = fopen(pathListOfRecoTxtDoc,'w');            
        fclose(fileID);
        % ------------------------------------------------------------------------------------


        % Retrieving an organized k-space -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
        if (dataInfo.type == "ParavisionRawFID")
            disp(' ')
            disp('Reconstruction Paravision raw FID files...')

            %[undersampledKSpace, ~, ~] = undersampledDataReconstructionFromMultiChannelRawFID(dataPath,origReadoutDim, origPhase1Dim, origPhase2Dim, nEchos, finalReadoutDim, finalPhase1Dim, finalPhase2Dim, nChannels, nVolumes, fsAcquisitionType, dataInfo.imageCorrection, dataInfo.fftshiftCorrection, dataInfo.offsetCorrection, dataInfo.showResults, dataInfo.slicesToShow(1), dataInfo.slicesToShow(2), dataInfo.slicesToShow(3), dataInfo.slicesToShow(5), dataInfo.sliceRotation(1), dataInfo.sliceRotation(2), dataInfo.sliceRotation(3), savePathAndFileName, pathListOfRecoTxtDoc);
            [undersampledKSpace, undersampledKSpaceView, undersampledZPImage, undersampledZPImageView, undersamplingPattern] = undersampledDataReconstructionFromMultiChannelRawFID(dataPath, undersamplingPatternFilePath, origReadoutDim, origPhase1Dim, origPhase2Dim, nEchos, finalReadoutDim, finalPhase1Dim, finalPhase2Dim, nChannels, nVolumes, acquisitionType, dataInfo.imageCorrection, dataInfo.fftshiftCorrection, dataInfo.offsetCorrection, dataInfo.showResults, dataInfo.slicesToShow(1), dataInfo.slicesToShow(2), dataInfo.slicesToShow(3), dataInfo.slicesToShow(5), dataInfo.sliceRotation(1), dataInfo.sliceRotation(2), dataInfo.sliceRotation(3), savePathAndFileName, pathListOfRecoTxtDoc);

            if(dataInfo.correctionBeforeReconstruction)
                % To code for undersampled k-spaces
                
                %{
                close all
                fprintf("\n####################################################\nCorrection before reconstruction\n##################################################");
                %ref_im = ifftshift(ifft(ifftshift(ifft(ifftshift(ifft(fullySampledKSpace,[],1),1),[],2),2),[],3),3);                                                
                ref_im = (ifft((ifft((ifft(fullySampledKSpace,[],1)),[],2)),[],3));                                                

                ref_im = circshift(ref_im,dataInfo.offsetCorrection);

                if (dataInfo.fftshiftCorrection(1))
                    ref_im = fftshift(ref_im,1);
                end
                if (dataInfo.fftshiftCorrection(2))
                    ref_im = fftshift(ref_im,2);
                end
                if (dataInfo.fftshiftCorrection(3))
                    ref_im = fftshift(ref_im,3);
                end
                if (dataInfo.fftshiftCorrection(4))
                    ref_im = fftshift(ref_im,4);
                end
                if (dataInfo.fftshiftCorrection(5))
                    ref_im = fftshift(ref_im,5);
                end


                fprintf("\n Showing images with k-space corrected before reconstruction...\n");
                showImage(ref_im, dataInfo.showResults, dataInfo.slicesToShow(1), dataInfo.slicesToShow(2), dataInfo.slicesToShow(3), dataInfo.slicesToShow(5), dataInfo.sliceRotation(1), dataInfo.sliceRotation(2), dataInfo.sliceRotation(3));                                

                pause(5)

                ImgSize = size(fullySampledKSpace);
                fullySampledKSpaceCorrected = FFT( ref_im, ones( ImgSize )); 
                showImage(fullySampledKSpaceCorrected, dataInfo.showResults, dataInfo.slicesToShow(1), dataInfo.slicesToShow(2), dataInfo.slicesToShow(3), dataInfo.slicesToShow(5), 0, 0, 0);

                dataInfo.offsetCorrection = [0,0,0,0,0];
                dataInfo.fftshiftCorrection = [0,0,0,0,0];

                fullySampledKSpace = fullySampledKSpaceCorrected;

                ref_im = ifftshift(ifft(ifftshift(ifft(ifftshift(ifft(fullySampledKSpace,[],1),1),[],2),2),[],3),3);                                                                                        
                showImage(ref_im, dataInfo.showResults, dataInfo.slicesToShow(1), dataInfo.slicesToShow(2), dataInfo.slicesToShow(3), dataInfo.slicesToShow(5), dataInfo.sliceRotation(1), dataInfo.sliceRotation(2), dataInfo.sliceRotation(3));                                

                %return

                % Enlever 
                % -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
                %fullySampledKSpaceCorrected = fftshift(fft(fftshift(fft(fftshift(fft(ref_im,[],1),1),[],2),2),[],3),3);                
                %{
                fullySampledKSpaceCorrected1 = (fft((fft((fft(ref_im,[],1)),[],2)),[],3));                
                fullySampledKSpaceCorrected2 = fftshift(fft(fftshift(fft(fftshift(fft(ref_im,[],1),1),[],2),2),[],3),3);                

                fullySampledKSpaceCorrectedView1 = squeeze(bart('rss 8', fullySampledKSpaceCorrected1));
                fullySampledKSpaceCorrectedView2 = squeeze(bart('rss 8', fullySampledKSpaceCorrected2));

                showImage(fullySampledKSpaceCorrectedView1, dataInfo.showResults, dataInfo.slicesToShow(1), dataInfo.slicesToShow(2), dataInfo.slicesToShow(3), dataInfo.slicesToShow(5), 0, 0, 0);                                
                showImage(fullySampledKSpaceCorrectedView2, dataInfo.showResults, dataInfo.slicesToShow(1), dataInfo.slicesToShow(2), dataInfo.slicesToShow(3), dataInfo.slicesToShow(5), 0, 0, 0);                                

                %fullySampledImageCorrected = bart('fft -i 7', fullySampledKSpaceCorrected);
                fullySampledImageCorrected1 = ifftshift(ifft(ifftshift(ifft(ifftshift(ifft(fullySampledKSpaceCorrected1,[],1),1),[],2),2),[],3),3);

                fullySampledImageCorrected2 = (ifft((ifft((ifft(fullySampledKSpaceCorrected2,[],1)),[],2)),[],3));
                %fullySampledImageCorrected = (ifft((ifft((ifft(fullySampledKSpaceCorrected,[],1)),[],2)),[],3));                
                fullySampledImageCorrectedView1 = bart('rss 8', fullySampledImageCorrected1);
                fullySampledImageCorrectedView2 = bart('rss 8', fullySampledImageCorrected2);

                fullySampledImageCorrectedView1 = squeeze(fullySampledImageCorrectedView1);
                fullySampledImageCorrectedView2 = squeeze(fullySampledImageCorrectedView2);

                showImage(fullySampledImageCorrectedView1, dataInfo.showResults, dataInfo.slicesToShow(1), dataInfo.slicesToShow(2), dataInfo.slicesToShow(3), dataInfo.slicesToShow(5), dataInfo.sliceRotation(1), dataInfo.sliceRotation(2), dataInfo.sliceRotation(3));                                
                showImage(fullySampledImageCorrectedView2, dataInfo.showResults, dataInfo.slicesToShow(1), dataInfo.slicesToShow(2), dataInfo.slicesToShow(3), dataInfo.slicesToShow(5), dataInfo.sliceRotation(1), dataInfo.sliceRotation(2), dataInfo.sliceRotation(3));                                
                %}
                % -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
                %}
            end

            disp('Undersampled data reconstruction and ZP saving process done!!!')
            disp('============================================================')


        elseif dataInfo.type == "k-spaceFromMatlabWorkspace"
            % To be coded !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            disp(' ')
            disp('Reading k-space from .mat file...')
            disp('K-space reading done!!!')
            disp(' ')
            % To be coded !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        else
            disp(' ')
            disp('Data type out of allowed cases!!!');
            disp(' ')
            %return                            
        end
        % -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------                                    
       

        %undersampledKSpace = fullySampledKSpace.*undersamplingPattern;
        

        fprintf('Undersampling k-space reconstruction done!!!\n\n')
        % --------------------------------------------------------------------------------------------------------


        for recoMeth=1:size(dataInfo.reconstructionMethod,2)

            reconstructionMethod = dataInfo.reconstructionMethod(recoMeth);


            % Reconstruction image from undersampled k-space ----------------------------------------------------------------------------------
            if (reconstructionMethod=="ConventionalCS")
                fprintf(strcat('Reconsctruction of undersampled k-space via ',reconstructionMethod,', using a',dataInfo.undersamplingPatternType,' undersampling pattern and ',dataInfo.undersamplingPatternStrategy,' strategy...\n\n'))

                savePathAndFileNameAfterReco = strcat(savePathAndFileName,"_AF-",num2str(dataInfo.AF),"_OverallAF-",num2str(dataInfo.actual4DAF),"_DWIAF-",num2str(dataInfo.correctedAF),"_",dataInfo.undersamplingPatternType,"_",dataInfo.undersamplingPatternStrategy,"_FSArea-",num2str(dataInfo.undersamplingPatternFullySampledRegionArea),"_",dataInfo.B0SamplingStrategy);

                [conventionalCSImageReconstructed, coilSensitivityMapsViewCorrected] = conventionalCSReconstruction_v2 (undersampledKSpace, dataInfo.CSParamLamba1, dataInfo.CSParamLamba2, dataInfo.CSParamNItrMax, dataInfo.imageCorrection, dataInfo.fftshiftCorrection, dataInfo.offsetCorrection, dataInfo.B0SamplingStrategy, nB0s, dataInfo.xCenterCubeForB0Norm(nSubject), dataInfo.yCenterCubeForB0Norm(nSubject), dataInfo.zCenterCubeForB0Norm(nSubject), dataInfo.cubeSizeForB0Norm(nSubject), savePathAndFileNameAfterReco, pathListOfRecoTxtDoc, dataInfo.showResults, dataInfo.slicesToShow(1), dataInfo.slicesToShow(2), dataInfo.slicesToShow(3));

                % Enlever
                % --------------------------------------------------------------------------------------------------------
                %[conventionalCSImageReconstructed1, coilSensitivityMapsViewCorrected1] = conventionalCSReconstruction_v2 (undersampledKSpace1, dataInfo.CSParamLamba1(1), dataInfo.CSParamLamba2(1), dataInfo.CSParamNItrMax(1), 0, dataInfo.fftshiftCorrection, dataInfo.offsetCorrection, B0SamplingStrategy, nB0s, dataInfo.xCenterCubeForB0Norm(nSubject), dataInfo.yCenterCubeForB0Norm(nSubject), dataInfo.zCenterCubeForB0Norm(nSubject), dataInfo.cubeSizeForB0Norm(nSubject), strcat(savePathAndFileNameAfterReco,"_1"), pathListOfRecoTxtDoc, dataInfo.showResults, dataInfo.slicesToShow(1), dataInfo.slicesToShow(2), dataInfo.slicesToShow(3));
                %[conventionalCSImageReconstructed2, coilSensitivityMapsViewCorrected2] = conventionalCSReconstruction_v2 (undersampledKSpace2, dataInfo.CSParamLamba1(1), dataInfo.CSParamLamba2(1), dataInfo.CSParamNItrMax(1), 0, dataInfo.fftshiftCorrection, dataInfo.offsetCorrection, B0SamplingStrategy, nB0s, dataInfo.xCenterCubeForB0Norm(nSubject), dataInfo.yCenterCubeForB0Norm(nSubject), dataInfo.zCenterCubeForB0Norm(nSubject), dataInfo.cubeSizeForB0Norm(nSubject), strcat(savePathAndFileNameAfterReco,"_2"), pathListOfRecoTxtDoc, dataInfo.showResults, dataInfo.slicesToShow(1), dataInfo.slicesToShow(2), dataInfo.slicesToShow(3));

                %save('conventionalCSImageReconstructed1.mat','conventionalCSImageReconstructed1');
                %save('conventionalCSImageReconstructed2.mat','conventionalCSImageReconstructed2');
                % ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
                fprintf('Conventional CS reconstruction done!!!\n\n')
                close all

            elseif (reconstructionMethod=="KLRCS")
                fprintf(strcat('Reconsctruction of undersampled k-space via ',reconstructionMethod,', using a',dataInfo.undersamplingPatternType,' undersampling pattern and ',dataInfo.undersamplingPatternStrategy,' strategy...\n\n'))                    
                savePathAndFileNameAfterReco = strcat(savePathAndFileName,"_AF-",num2str(dataInfo.AF),"_OverallAF-",num2str(dataInfo.actual4DAF),"_DWIAF-",num2str(dataInfo.correctedAF),"_",dataInfo.undersamplingPatternType,"_",dataInfo.undersamplingPatternStrategy,"_FSArea-",num2str(dataInfo.undersamplingPatternFullySampledRegionArea),"_",dataInfo.B0SamplingStrategy);                                       
                
                
                [klrCSImageReconstructed] = KLRCSReconstruction_v2 (undersampledKSpace, undersamplingPattern, dataInfo.KLRCSParamLambdaReg, dataInfo.KLRCSParamALMReg, dataInfo.KLRCSParamNItrMax, dataInfo.imageCorrection, dataInfo.correctionBeforeReconstruction, dataInfo.fftshiftCorrection, dataInfo.offsetCorrection, dataInfo.B0SamplingStrategy, nB0s, dataInfo.xCenterCubeForB0Norm(nSubject), dataInfo.yCenterCubeForB0Norm(nSubject), dataInfo.zCenterCubeForB0Norm(nSubject), dataInfo.cubeSizeForB0Norm(nSubject), savePathAndFileNameAfterReco, pathListOfRecoTxtDoc, dataInfo.showResults, dataInfo.slicesToShow(1), dataInfo.slicesToShow(2), dataInfo.slicesToShow(3));

                fprintf('KLR CS reconstruction done!!!\n\n')

            elseif (reconstructionMethod=="ZeroPadding")
                
                fprintf(strcat('Reconsctruction of undersampled k-space via ',reconstructionMethod,', using a',dataInfo.undersamplingPatternType,' undersampling pattern and ',dataInfo.undersamplingPatternStrategy,' strategy...\n\n'))                    
                    
                savePathAndFileNameAfterReco = strcat(savePathAndFileName,"_AF-",num2str(dataInfo.AF),"_OverallAF-",num2str(dataInfo.actual4DAF),"_DWIAF-",num2str(dataInfo.correctedAF),"_",dataInfo.undersamplingPatternType,"_",dataInfo.undersamplingPatternStrategy,"_FSArea-",num2str(dataInfo.undersamplingPatternFullySampledRegionArea),"_",dataInfo.B0SamplingStrategy);                                       

                [zpImageReconstructed] = zeroPaddingReconstruction (undersampledKSpace, dataInfo.imageCorrection, dataInfo.fftshiftCorrection, dataInfo.offsetCorrection, dataInfo.showResults, dataInfo.slicesToShow(1), dataInfo.slicesToShow(2), dataInfo.slicesToShow(3), dataInfo.slicesToShow(5), dataInfo.sliceRotation(1), dataInfo.sliceRotation(2), dataInfo.sliceRotation(3), savePathAndFileNameAfterReco, pathListOfRecoTxtDoc);

                fprintf('Zero padding reconstruction done!!!\n\n')

            end


            % --------------------------------------------------------------------------------------------------------

        end
        
    end
    
    
    
    

else
    disp(' ')
    disp('Option chosen does not match any case!!! Please, choose one among the available options!!!')
    disp(' ')

end


