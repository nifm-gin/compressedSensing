function [undersampledKSpace, undersampledKSpaceView, undersampledZPImage, undersampledZPImageView, undersamplingPatterRetrieved] = undersampledDataReconstructionFromMultiChannelRawFID(dataPath, acquisitionCSMaskTxtPath, origReadoutDim, origPhase1Dim, origPhase2Dim, NEchos, finalReadoutDim, finalPhase1Dim, finalPhase2Dim, nChannels, nVolumes, acquisitionType, correctImage, fftShiftCorrection, offsetCorrection, showResults, whichSliceX, whichSliceY, whichSliceZ, whichVolume, rotAngleSliceX, rotAngleSliceY, rotAngleSliceZ, savePathAndFileNameundersampled, pathTxtDoc)

% -------------------------------------------------------------------------------------------------------------------------
% -------------------------------------------------------------------------------------------------------------------------
% Fully sampled -----------------------------------------------------------------------------------------------------------
% -------------------------------------------------------------------------------------------------------------------------
% -------------------------------------------------------------------------------------------------------------------------

% Redifining whichSliceX for first view (wihtout fftshift correction) ---
whichSliceXTemp = whichSliceX;
whichSliceX = 1;
% -----------------------------------------------------------------------

undersampledKSpace = zeros(finalReadoutDim,finalPhase1Dim,finalPhase2Dim,nChannels,nVolumes);
undersamplingPatterRetrieved = zeros(finalReadoutDim,finalPhase1Dim,finalPhase2Dim,nChannels,nVolumes);

fprintf("Number of channels to reconstruct: %d\n",nChannels);
for channel=0:(nChannels-1)    
    fprintf("Reconstructing channel %d\n", channel);
    fileID = strcat(dataPath,num2str(channel));
    [undersampledKSpacePerChannel, ~, csMaskRetrieved] =  readParavisionUndersampledRawFID (fileID, acquisitionCSMaskTxtPath, origReadoutDim, origPhase1Dim, origPhase2Dim, NEchos, finalReadoutDim, finalPhase1Dim, finalPhase2Dim, nChannels, nVolumes, acquisitionType, showResults, whichSliceX, whichSliceY, whichSliceZ, whichVolume);           
    undersampledKSpace(:,:,:,channel+1,:) = squeeze(undersampledKSpacePerChannel);    
    undersamplingPatterRetrieved(:,:,:,channel+1,:) = squeeze(csMaskRetrieved);
    
end

%[undersampledKSpace0, undersampledImage0] =  readParavisionUndersampledRawFID (fileID, origReadoutDim, origPhase1Dim, origPhase2Dim, NEchos, finalReadoutDim, finalPhase1Dim, finalPhase2Dim, nChannels, nVolumes, dataTypeFS, showResults, whichSliceX, whichSliceY, whichSliceZ, whichVolume);


%fileID = strcat(dataPath,'1');
%[undersampledKSpace1, ~, csMaskRetrieved1] =  readParavisionUndersampledRawFID (fileID, acquisitionCSMaskTxtPath, origReadoutDim, origPhase1Dim, origPhase2Dim, NEchos, finalReadoutDim, finalPhase1Dim, finalPhase2Dim, nChannels, nVolumes, acquisitionType, showResults, whichSliceX, whichSliceY, whichSliceZ, whichVolume);




%undersampledKSpace(:,:,:,1,:) = squeeze(undersampledKSpace0);
%undersampledKSpace(:,:,:,2,:) = squeeze(undersampledKSpace1);
%undersampledKSpace(:,:,:,3,:) = squeeze(undersampledKSpace2);
%undersampledKSpace(:,:,:,4,:) = squeeze(undersampledKSpace3);


% Redifining whichSliceX for correct view ---
whichSliceX = whichSliceXTemp;
% -------------------------------------------

% Performing phase correction -------------------------------------------------------------------------------------------------------------------------------------------------
%[undersampledKSpace, recoImageCorrected] =  PhaseCorrection (undersampledKSpace, 2, phaseCorrectionValue, 1, whichSliceX, whichSliceY, whichSliceZ, whichChannel, whichVolume);
% -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------

% Making k-space for visualization ------------------------
undersampledKSpaceView = squeeze(bart('rss 8', undersampledKSpace));
% ---------------------------------------------------------

% Creating image ----------------------------------------
undersampledZPImage = bart('fft -i 7', undersampledKSpace);
undersampledZPImageView = bart('rss 8', undersampledZPImage);
undersampledZPImageView = squeeze(undersampledZPImageView);
% -------------------------------------------------------


if (correctImage)
    % Correcting image with FFTshift --------------------------------------------
    if (fftShiftCorrection(1))
        undersampledZPImageView = fftshift(undersampledZPImageView,1);
    end
    if (fftShiftCorrection(2))
        undersampledZPImageView = fftshift(undersampledZPImageView,2);
    end
    if (fftShiftCorrection(3))
        undersampledZPImageView = fftshift(undersampledZPImageView,3);
    end
    % ---------------------------------------------------------------------------
    
    
    % Correct image offset ----------------------------------------------------
    undersampledZPImageView = circshift(undersampledZPImageView, offsetCorrection);
    % -------------------------------------------------------------------------

end



if (showResults)
    close all

    figure;
    colormap gray
    subplot(2,3,1)
    imagesc(abs(squeeze(undersampledKSpaceView(:,:,whichSliceZ,whichVolume))).^0.25)

    subplot(2,3,2)
    imagesc(abs(squeeze(undersampledKSpaceView(:,whichSliceY,:,whichVolume))).^0.25)

    subplot(2,3,3)
    imagesc(abs(squeeze(undersampledKSpaceView(whichSliceX,:,:,whichVolume))).^0.25)

    subplot(2,3,4)
    imagesc(abs(squeeze(undersampledZPImageView(:,:,whichSliceZ,whichVolume))))

    subplot(2,3,5)
    imagesc(abs(squeeze(undersampledZPImageView(:,whichSliceY,:,whichVolume))))

    subplot(2,3,6)
    imagesc(abs(squeeze(undersampledZPImageView(whichSliceX,:,:,whichVolume))))



    figure;
    sgtitle('Rotated views')
    colormap gray
    subplot(1,3,1)
    imagesc(imrotate(abs(squeeze(undersampledZPImageView(:,:,whichSliceZ,whichVolume))),rotAngleSliceZ))

    subplot(1,3,2)
    imagesc(imrotate(abs(squeeze(undersampledZPImageView(:,whichSliceY,:,whichVolume))),rotAngleSliceY))

    subplot(1,3,3)
    imagesc(imrotate(abs(squeeze(undersampledZPImageView(whichSliceX,:,:,whichVolume))),rotAngleSliceX))

end

disp(' ')
disp("Saving undersampled image...");
% Save .nii image in a specific folder -------------------------------------------
niftiwrite(single(abs(undersampledZPImageView)),strcat(savePathAndFileNameundersampled,'_ZP-reconstruction.nii'));
niftiwrite(single(abs(fftshift(undersampledZPImage,1))),strcat(savePathAndFileNameundersampled,'_ZP-reconstruction_4-ch.nii'))
disp(strcat("Done!!! Undersampled image saved as NifTI in: ",savePathAndFileNameundersampled));
disp(' ')
% --------------------------------------------------------------------------------

fileID = fopen(pathTxtDoc,'a');
fprintf(fileID,'%s\r\n',strcat(savePathAndFileNameundersampled,'_ZP-reconstruction.nii'));
fclose(fileID);
fprintf("\nPath of reconstructed undersampled image added in a txt file.\n")

% Pause step 1 ---
pause(2)
% ----------------