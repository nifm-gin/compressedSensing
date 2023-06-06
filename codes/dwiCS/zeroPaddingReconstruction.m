function [undersampledImageView] = zeroPaddingReconstruction (undersampledKSpace, correctImage, fftShiftCorrection, offsetCorrection, showResults, whichSliceX, whichSliceY, whichSliceZ, whichVolume, rotAngleSliceX, rotAngleSliceY, rotAngleSliceZ, savePathAndFileNameUndersampled, pathTxtDoc)

% Making k-space for visualization ------------------------
undersampledKSpaceView = squeeze(bart('rss 8', undersampledKSpace));
% ---------------------------------------------------------

% Creating image ----------------------------------------
undersampledImage = bart('fft -i 7', undersampledKSpace);
undersampledImageView = bart('rss 8', undersampledImage);
undersampledImageView = squeeze(undersampledImageView);
% -------------------------------------------------------


if (correctImage)
    % Correcting image with FFTshift --------------------------------------------
    if (fftShiftCorrection(1))
        undersampledImageView = fftshift(undersampledImageView,1);
    end
    if (fftShiftCorrection(2))
        undersampledImageView = fftshift(undersampledImageView,2);
    end
    if (fftShiftCorrection(3))
        undersampledImageView = fftshift(undersampledImageView,3);
    end
    % ---------------------------------------------------------------------------
    
    
    % Correct image offset ----------------------------------------------------
    undersampledImageView = circshift(undersampledImageView, offsetCorrection);
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
    imagesc(abs(squeeze(undersampledImageView(:,:,whichSliceZ,whichVolume))))

    subplot(2,3,5)
    imagesc(abs(squeeze(undersampledImageView(:,whichSliceY,:,whichVolume))))

    subplot(2,3,6)
    imagesc(abs(squeeze(undersampledImageView(whichSliceX,:,:,whichVolume))))



    figure;
    sgtitle('Rotated views')
    colormap gray
    subplot(1,3,1)
    imagesc(imrotate(abs(squeeze(undersampledImageView(:,:,whichSliceZ,whichVolume))),rotAngleSliceZ))

    subplot(1,3,2)
    imagesc(imrotate(abs(squeeze(undersampledImageView(:,whichSliceY,:,whichVolume))),rotAngleSliceY))

    subplot(1,3,3)
    imagesc(imrotate(abs(squeeze(undersampledImageView(whichSliceX,:,:,whichVolume))),rotAngleSliceX))

end

disp(' ')
disp("Saving zero-padding image...");
% Save .nii image in a specific folder -------------------------------------------
niftiwrite(single(abs(undersampledImageView)),strcat(savePathAndFileNameUndersampled,'_ZeroPadding.nii'));
niftiwrite(single(abs(fftshift(undersampledImage,1))),strcat(savePathAndFileNameUndersampled,'_ZeroPadding_4-ch.nii'))
disp(strcat("Done!!! Zero-padding image saved as NifTI in: ",savePathAndFileNameUndersampled));
disp(' ')
% --------------------------------------------------------------------------------

fileID = fopen(pathTxtDoc,'a');
fprintf(fileID,'%s\r\n',strcat(savePathAndFileNameUndersampled,'_ZeroPadding.nii'));
fclose(fileID);
fprintf("\nPath of reconstructed zero-padding image added in a txt file.\n")

end