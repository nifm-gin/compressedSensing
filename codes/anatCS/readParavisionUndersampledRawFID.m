function [kSpaceReconstructed, imageZPReconstructed, csMaskRetrieved] =  readParavisionUndersampledRawFID (fidPath, acquisitionCSMaskTxtPath, origReadoutDim, origPhase1Dim, origPhase2Dim, NEchos, finalReadoutDim, finalPhase1Dim, finalPhase2Dim, nChannels, nVolumes, acquisitionType, showResults, whichSliceX, whichSliceY, whichSliceZ, whichVolume)

    switch acquisitionType
        case 'CS-3D-FLASH'
            disp('Compressed Sensing 3D Flash selected')
        
            
            fileid = fopen(fidPath,'r','ieee-le');                        

            imagData = fread(fileid,[inf],'float64');  % ACQ_word_size vaut 'int16' ou 'int32' en général, sa valeur est écrite dans le fichier ACQP
            fclose(fileid);

            imagData = imagData(1:2:end)+i*imagData(2:2:end); %Séparer les parties réelles et les parties imaginaires

            figure;
            plot(1:size(imagData,1),abs(imagData));
            title('Bruker FID');



           
            %% Read corresponding txt file and organization of data
            originalkSpaceReconstructed = zeros(origReadoutDim,finalPhase1Dim,finalPhase2Dim);

            fileID = fopen(acquisitionCSMaskTxtPath,'r');
            allData = fscanf(fileID,'%f');
            fclose(fileID);

            whichCase = allData(1);
            readoutMatrixSize = allData(2);
            phase1MatrixSize = allData(3);
            phase2MatrixSize = allData(4);

            %numberOfSlices = allData(5);
            numberOfSlices = phase2MatrixSize;
            nVolumes = allData(5);
            accelerationFactor = allData(6);
            centerSquareSize = allData(7);

            variableDensityOption = allData(8);
            ellipseOption = allData(9);
            seedValue = allData(10);

            totalNumberOfLines = allData(11);

            posInAllDataOfFirstSlice = 12;
            
            
            %% Data organization
            disp('ImagData size')
            size(imagData)            
            origReadoutDim
            disp('ImagData size / ReadOut dim')
            size(imagData,1)/origReadoutDim
            totalNumberOfLines
            originalOrganizedData = zeros(origReadoutDim,totalNumberOfLines);
            
            % If fewer lines were acquired -----------------------------------------------------------
            %originalOrganizedDataTemp = reshape(imagData, origReadoutDim, []);
            %originalOrganizedData(:,1:size(originalOrganizedDataTemp,2)) = originalOrganizedDataTemp;
            % ----------------------------------------------------------------------------------------           
            
            %imagData = imagData(1:(origReadoutDim*totalNumberOfLines))
            
            originalOrganizedData = reshape(imagData, origReadoutDim, totalNumberOfLines);
            
            figure;
            imagesc(abs(squeeze(originalOrganizedData)).^0.25);
            title('Original lines acquired');

            %organizedData = originalOrganizedData(1:finalReadoutDim,:);
            
            

            posData = 1;
            globalPos = posInAllDataOfFirstSlice;
            for volume=1:nVolumes
                
                for readPos=1:numberOfSlices
                    whichSlice = allData(globalPos);
                    globalPos = globalPos + 1;

                    numberOfLinesCurrentSlice = allData(globalPos);
                    globalPos = globalPos + 1;

                    for posLine=1:numberOfLinesCurrentSlice
                        whichLine = allData(globalPos);
                        globalPos = globalPos + 1;

                        originalkSpaceReconstructed(:,whichLine,whichSlice,volume) = originalOrganizedData(:,posData);
                        posData = posData + 1;
                        posLine;
                    end


                end
            end
            
            kSpaceReconstructed = originalkSpaceReconstructed(1:finalReadoutDim,:,:);
            

            %kSpaceReconstructed = originalOrganizedData(1:finalReadoutDim,:,round((origPhase2Dim-finalPhase2Dim)/2):(round((origPhase2Dim-finalPhase2Dim)/2)+finalPhase2Dim-1));
            imageZPReconstructed = fftshift(ifft(fftshift(kSpaceReconstructed,1),[],1),1);
            imageZPReconstructed = fftshift(ifft(fftshift(imageZPReconstructed,2),[],2),2);
            imageZPReconstructed = fftshift(ifft(fftshift(imageZPReconstructed,3),[],3),3);


            
            if (showResults)
                disp('Showing results ...');
                
                figure;
                subplot(1,3,1)
                imagesc(abs(squeeze(originalkSpaceReconstructed(:,:,whichSliceZ))).^0.25);
                title({'Original undersampled k-space','(with readout zero-padding by Paravision)','X-Y plan'});

                subplot(1,3,2)
                imshow(abs(squeeze(kSpaceReconstructed(:,:,whichSliceZ))).^0.25,[]);
                title({'Undersampled K-space cropped','(without readout zero-padding by Paravision)','X-Y plan'});

                subplot(1,3,3)
                imshow(abs(squeeze(imageZPReconstructed(:,:,whichSliceZ))),[]);
                title({'Zero-filling image','X-Y plan'});
                
                
                figure;
                subplot(1,3,1)
                imagesc(abs(squeeze(originalkSpaceReconstructed(:,whichSliceY,:))).^0.25);
                title({'Original undersampled k-space','(with readout zero-padding by Paravision)','X-Z plan'});

                subplot(1,3,2)
                imshow(abs(squeeze(kSpaceReconstructed(:,whichSliceY,:))).^0.25,[]);
                title({'Undersampled K-space cropped','(without readout zero-padding by Paravision)','X-Z plan'});

                subplot(1,3,3)
                imshow(abs(squeeze(imageZPReconstructed(:,whichSliceY,:))),[]);
                title({'Zero-filling image','X-Z plan'});
                
                
                figure;
                subplot(1,3,1)
                imagesc(abs(squeeze(originalkSpaceReconstructed(whichSliceX,:,:))).^0.25);
                title({'Original undersampled k-space','(with readout zero-padding by Paravision)','Y-Z plan'});

                subplot(1,3,2)
                imshow(abs(squeeze(kSpaceReconstructed(whichSliceX,:,:))).^0.25,[]);
                title({'Undersampled K-space cropped','(without readout zero-padding by Paravision)','Y-Z plan'});

                subplot(1,3,3)
                imshow(abs(squeeze(imageZPReconstructed(whichSliceX,:,:))),[]);
                title({'Zero-filling image','Y-Z plan'});
            
                csMaskRetrieved = squeeze(abs(kSpaceReconstructed) > 0);
                figure;
                imshow(squeeze(csMaskRetrieved(whichSliceX,:,:)));
                title('Compressed Sensing mask retrieved');
                
            else                
                disp('Show results not selected!');
                
            end
        
        
        
        
        
        case 'CS-2D-FLASH'
            disp('Compressed Sensing 2D-MS Flash selected')
        
        
        case 'CS-3D-SE-DTI_DWLoopFirst'
            disp('Compressed Sensing 3D DTI Spin Echo - Same undersampling pattern for all 3D volumes (Mono3D masks) - DW loop as the first one')
        
            
            fileid = fopen(fidPath,'r','ieee-le');                        

            imagData = fread(fileid,[inf],'float64');  % ACQ_word_size vaut 'int16' ou 'int32' en général, sa valeur est écrite dans le fichier ACQP
            fclose(fileid);

            imagData = imagData(1:2:end)+i*imagData(2:2:end); %Séparer les parties réelles et les parties imaginaires

            figure;
            plot(1:size(imagData,1),abs(imagData));
            title('Bruker FID');



           
            %% Read corresponding txt file and organization of data
            originalkSpaceReconstructed = zeros(origReadoutDim,finalPhase1Dim,finalPhase2Dim,nVolumes);

            fileID = fopen(acquisitionCSMaskTxtPath,'r');
            allData = fscanf(fileID,'%f');
            fclose(fileID);

            whichCase = allData(1);
            readoutMatrixSize = allData(2);
            phase1MatrixSize = allData(3);
            phase2MatrixSize = allData(4);

            numberOfSlices = allData(5);
            accelerationFactor = allData(6);
            centerSquareSize = allData(7);

            variableDensityOption = allData(8);
            ellipseOption = allData(9);
            seedValue = allData(10);

            totalNumberOfLines = allData(11);

            posInAllDataOfFirstSlice = 12;
            
            
            %% Data organization
            size(imagData)
            origReadoutDim
            totalNumberOfLines
            
            
            % To do - Make acquisitions always with right number of lines
            % !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            % or less in order to vanish Temp variable and having always the same sumber of lines than "totalNumberOfLines"          
            %originalOrganizedData = zeros(origReadoutDim,totalNumberOfLines,nVolumes);
            %originalOrganizedDataTemp = reshape(imagData, origReadoutDim, totalNumberOfLines, nVolumes);
            %originalOrganizedData(:,1:size(originalOrganizedDataTemp,2),:) = originalOrganizedDataTemp;
            % !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            
            originalOrganizedData = reshape(imagData, origReadoutDim, nVolumes, totalNumberOfLines);
            originalOrganizedData = permute(originalOrganizedData,[1 3 2]);
            
            figure;
            imagesc(abs(squeeze(originalOrganizedData(:,:,whichVolume))).^0.25);
            title('Original lines acquired');

            %organizedData = originalOrganizedData(1:finalReadoutDim,:);
            
            

            
            for volume=1:nVolumes
                posData = 1;
                globalPos = posInAllDataOfFirstSlice;
                for readPos=1:numberOfSlices
                    whichSlice = allData(globalPos);
                    globalPos = globalPos + 1;

                    numberOfLinesCurrentSlice = allData(globalPos);
                    globalPos = globalPos + 1;

                    for posLine=1:numberOfLinesCurrentSlice
                        whichLine = allData(globalPos);
                        globalPos = globalPos + 1;

                        originalkSpaceReconstructed(:,whichLine,whichSlice,volume) = originalOrganizedData(:,posData,volume);
                        posData = posData + 1;
                        posLine;
                    end


                end
            end
            
            kSpaceReconstructed = originalkSpaceReconstructed(1:finalReadoutDim,:,:,:);
            

            %kSpaceReconstructed = originalOrganizedData(1:finalReadoutDim,:,round((origPhase2Dim-finalPhase2Dim)/2):(round((origPhase2Dim-finalPhase2Dim)/2)+finalPhase2Dim-1));
            imageZPReconstructed = fftshift(ifft(fftshift(kSpaceReconstructed,1),[],1),1);
            imageZPReconstructed = fftshift(ifft(fftshift(imageZPReconstructed,2),[],2),2);
            imageZPReconstructed = fftshift(ifft(fftshift(imageZPReconstructed,3),[],3),3);


            
            if (showResults)
                disp('Showing results ...');
                
                figure;
                subplot(1,3,1)
                imagesc(abs(squeeze(originalkSpaceReconstructed(:,:,whichSliceZ,whichVolume))).^0.25);
                title({'Original undersampled k-space','(with readout zero-padding by Paravision)','X-Y plan'});

                subplot(1,3,2)
                imshow(abs(squeeze(kSpaceReconstructed(:,:,whichSliceZ,whichVolume))).^0.25,[]);
                title({'Undersampled K-space cropped','(without readout zero-padding by Paravision)','X-Y plan'});

                subplot(1,3,3)
                imshow(abs(squeeze(imageZPReconstructed(:,:,whichSliceZ,whichVolume))),[]);
                title({'Zero-filling image','X-Y plan'});
                
                
                figure;
                subplot(1,3,1)
                imagesc(abs(squeeze(originalkSpaceReconstructed(:,whichSliceY,:,whichVolume))).^0.25);
                title({'Original undersampled k-space','(with readout zero-padding by Paravision)','X-Z plan'});

                subplot(1,3,2)
                imshow(abs(squeeze(kSpaceReconstructed(:,whichSliceY,:,whichVolume))).^0.25,[]);
                title({'Undersampled K-space cropped','(without readout zero-padding by Paravision)','X-Z plan'});

                subplot(1,3,3)
                imshow(abs(squeeze(imageZPReconstructed(:,whichSliceY,:,whichVolume))),[]);
                title({'Zero-filling image','X-Z plan'});
                
                
                figure;
                subplot(1,3,1)
                imagesc(abs(squeeze(originalkSpaceReconstructed(whichSliceX,:,:,whichVolume))).^0.25);
                title({'Original undersampled k-space','(with readout zero-padding by Paravision)','Y-Z plan'});

                subplot(1,3,2)
                imshow(abs(squeeze(kSpaceReconstructed(whichSliceX,:,:,whichVolume))).^0.25,[]);
                title({'Undersampled K-space cropped','(without readout zero-padding by Paravision)','Y-Z plan'});

                subplot(1,3,3)
                imshow(abs(squeeze(imageZPReconstructed(whichSliceX,:,:,whichVolume))),[]);
                title({'Zero-filling image','Y-Z plan'});
            
                csMaskRetrieved = squeeze(abs(kSpaceReconstructed) > 0);
                figure;
                imshow(squeeze(csMaskRetrieved(whichSliceX,:,:,whichVolume)),[]);
                title('Compressed Sensing mask retrieved - Y-Z plan - Volume 1');
                
            else                
                disp('Show results not selected!');
                
            end
            
        
        case 'CS-3D-SE-DTI_DWLoopLast'
            disp('Compressed Sensing 3D DTI Spin Echo - One undersampling pattern per 3D volume (Multi3D masks) - DW loop as the last one')
        
            
            fileid = fopen(fidPath,'r','ieee-le');                        

            imagData = fread(fileid,[inf],'float64');  % ACQ_word_size vaut 'int16' ou 'int32' en général, sa valeur est écrite dans le fichier ACQP
            fclose(fileid);

            imagData = imagData(1:2:end)+i*imagData(2:2:end); %Séparer les parties réelles et les parties imaginaires

            figure;
            plot(1:size(imagData,1),abs(imagData));
            title('Bruker FID');
            %pause(10)


           
            %% Read corresponding txt file and organization of data
            originalkSpaceReconstructed = zeros(origReadoutDim,finalPhase1Dim,finalPhase2Dim,nVolumes);

            fileID = fopen(acquisitionCSMaskTxtPath,'r');
            allData = fscanf(fileID,'%f');
            fclose(fileID);

            whichCase = allData(1);
            readoutMatrixSize = allData(2);
            phase1MatrixSize = allData(3);
            phase2MatrixSize = allData(4);

            %numberOfSlices = allData(5);
            numberOfSlices = phase2MatrixSize;
            nVolumes = allData(5);
            accelerationFactor = allData(6);
            centerSquareSize = allData(7);

            variableDensityOption = allData(8);
            ellipseOption = allData(9);
            seedValue = allData(10);

            totalNumberOfLines = allData(11);

            posInAllDataOfFirstSlice = 12;
            
            
            %% Data organization
            size(imagData)
            origReadoutDim
            totalNumberOfLines
            
            
            % To do - Make acquisitions always with right number of lines
            % !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            % or less in order to vanish Temp variable and having always the same sumber of lines than "totalNumberOfLines"          
            %originalOrganizedData = zeros(origReadoutDim,totalNumberOfLines,nVolumes);
            %originalOrganizedDataTemp = reshape(imagData, origReadoutDim, totalNumberOfLines, nVolumes);
            %originalOrganizedData(:,1:size(originalOrganizedDataTemp,2),:) = originalOrganizedDataTemp;
            % !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            
            %originalOrganizedData = reshape(imagData, origReadoutDim, nVolumes, totalNumberOfLines);
            
            originalOrganizedData = reshape(imagData, origReadoutDim, totalNumberOfLines);
            %originalOrganizedData = reshape(imagData(1:(end-origReadoutDim)), origReadoutDim, totalNumberOfLines);
           
            %originalOrganizedData = permute(originalOrganizedData,[1 3 2]);
            
            %figure;
            %imagesc(abs(squeeze(originalOrganizedData(:,:,whichVolume))).^0.25);
            %title('Original lines acquired');

            %organizedData = originalOrganizedData(1:finalReadoutDim,:);
            
            

            posData = 1;
            globalPos = posInAllDataOfFirstSlice;
            for volume=1:nVolumes
                
                for readPos=1:numberOfSlices
                    whichSlice = allData(globalPos);
                    globalPos = globalPos + 1;

                    numberOfLinesCurrentSlice = allData(globalPos);
                    globalPos = globalPos + 1;

                    for posLine=1:numberOfLinesCurrentSlice
                        whichLine = allData(globalPos);
                        globalPos = globalPos + 1;

                        originalkSpaceReconstructed(:,whichLine,whichSlice,volume) = originalOrganizedData(:,posData);
                        posData = posData + 1;
                        posLine;
                    end


                end
            end
            
            kSpaceReconstructed = originalkSpaceReconstructed(1:finalReadoutDim,:,:,:);
            

            %kSpaceReconstructed = originalOrganizedData(1:finalReadoutDim,:,round((origPhase2Dim-finalPhase2Dim)/2):(round((origPhase2Dim-finalPhase2Dim)/2)+finalPhase2Dim-1));
            imageZPReconstructed = fftshift(ifft(fftshift(kSpaceReconstructed,1),[],1),1);
            imageZPReconstructed = fftshift(ifft(fftshift(imageZPReconstructed,2),[],2),2);
            imageZPReconstructed = fftshift(ifft(fftshift(imageZPReconstructed,3),[],3),3);


            
            if (showResults)
                disp('Showing results ...');
                
                figure;
                subplot(1,3,1)
                imagesc(abs(squeeze(originalkSpaceReconstructed(:,:,whichSliceZ,whichVolume))).^0.25);
                title({'Original undersampled k-space','(with readout zero-padding by Paravision)','X-Y plan'});

                subplot(1,3,2)
                imshow(abs(squeeze(kSpaceReconstructed(:,:,whichSliceZ,whichVolume))).^0.25,[]);
                title({'Undersampled K-space cropped','(without readout zero-padding by Paravision)','X-Y plan'});

                subplot(1,3,3)
                imshow(abs(squeeze(imageZPReconstructed(:,:,whichSliceZ,whichVolume))),[]);
                title({'Zero-filling image','X-Y plan'});
                
                
                figure;
                subplot(1,3,1)
                imagesc(abs(squeeze(originalkSpaceReconstructed(:,whichSliceY,:,whichVolume))).^0.25);
                title({'Original undersampled k-space','(with readout zero-padding by Paravision)','X-Z plan'});

                subplot(1,3,2)
                imshow(abs(squeeze(kSpaceReconstructed(:,whichSliceY,:,whichVolume))).^0.25,[]);
                title({'Undersampled K-space cropped','(without readout zero-padding by Paravision)','X-Z plan'});

                subplot(1,3,3)
                imshow(abs(squeeze(imageZPReconstructed(:,whichSliceY,:,whichVolume))),[]);
                title({'Zero-filling image','X-Z plan'});
                
                
                figure;
                subplot(1,3,1)
                imagesc(abs(squeeze(originalkSpaceReconstructed(whichSliceX,:,:,whichVolume))).^0.25);
                title({'Original undersampled k-space','(with readout zero-padding by Paravision)','Y-Z plan'});

                subplot(1,3,2)
                imshow(abs(squeeze(kSpaceReconstructed(whichSliceX,:,:,whichVolume))).^0.25,[]);
                title({'Undersampled K-space cropped','(without readout zero-padding by Paravision)','Y-Z plan'});

                subplot(1,3,3)
                imshow(abs(squeeze(imageZPReconstructed(whichSliceX,:,:,whichVolume))),[]);
                title({'Zero-filling image','Y-Z plan'});
            
                csMaskRetrieved = squeeze(abs(kSpaceReconstructed) > 0);
                figure;
                imshow(squeeze(csMaskRetrieved(whichSliceX,:,:,whichVolume)),[]);
                title('Compressed Sensing mask retrieved - Y-Z plan - Volume 1');                
                
            else                
                disp('Show results not selected!');
                
            end
            
            
        case 'CS-3D-SE-DTI_DWLoopFirst_BinaryGuideArray'
            disp('Compressed Sensing 3D DTI Spin Echo - One undersampling pattern per 3D volume (Multi3D masks) - DW loop as first one with binary guide array')
        
            
            fileid = fopen(fidPath,'r','ieee-le');                        

            imagData = fread(fileid,[inf],'float64');  % ACQ_word_size vaut 'int16' ou 'int32' en général, sa valeur est écrite dans le fichier ACQP
            fclose(fileid);

            imagData = imagData(1:2:end)+i*imagData(2:2:end); %Séparer les parties réelles et les parties imaginaires

            figure;
            plot(1:size(imagData,1),abs(imagData));
            title('Bruker FID');
            %pause(10)


           
            %% Read corresponding txt file and organization of data
            originalkSpaceReconstructed = zeros(origReadoutDim,finalPhase1Dim,finalPhase2Dim,nVolumes);

            fileID = fopen(acquisitionCSMaskTxtPath,'r');
            allData = fscanf(fileID,'%f');
            fclose(fileID);

            whichCase = allData(1);
            readoutMatrixSize = allData(2);
            phase1MatrixSize = allData(3);
            phase2MatrixSize = allData(4);

            %numberOfSlices = allData(5);
            numberOfSlices = phase2MatrixSize;
            nB0s = allData(5);
            nDiffusionDirections = allData(6);
            nVolumes = nB0s + nDiffusionDirections;
            accelerationFactor = allData(7);
            centerSquareSize = allData(8);

            variableDensityOption = allData(9);
            ellipseOption = allData(10);
            seedValue = allData(11);

            totalNumberOfLines = allData(12);

            posInAllDataOfFirstSlice = 13;
            
            
            %% Data organization
            size(imagData)
            origReadoutDim
            totalNumberOfLines
            
            
            % To do - Make acquisitions always with right number of lines
            % !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            % or less in order to vanish Temp variable and having always the same sumber of lines than "totalNumberOfLines"          
            %originalOrganizedData = zeros(origReadoutDim,totalNumberOfLines,nVolumes);
            %originalOrganizedDataTemp = reshape(imagData, origReadoutDim, totalNumberOfLines, nVolumes);
            %originalOrganizedData(:,1:size(originalOrganizedDataTemp,2),:) = originalOrganizedDataTemp;
            % !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            
            %originalOrganizedData = reshape(imagData, origReadoutDim, nVolumes, totalNumberOfLines);
            
            originalOrganizedData = reshape(imagData, origReadoutDim, totalNumberOfLines);
            %originalOrganizedData = reshape(imagData(1:(end-origReadoutDim)), origReadoutDim, totalNumberOfLines);
           
            %originalOrganizedData = permute(originalOrganizedData,[1 3 2]);
            
            %figure;
            %imagesc(abs(squeeze(originalOrganizedData(:,:,whichVolume))).^0.25);
            %title('Original lines acquired');

            %organizedData = originalOrganizedData(1:finalReadoutDim,:);
            
            
            globalPos = posInAllDataOfFirstSlice;
            posData = 1;
            for whichSlice=1:finalPhase2Dim
                
                for whichLine=1:finalPhase1Dim
                    
                    for volume=1:nVolumes
                        
                        acquisitionMade = allData(globalPos);
                        
                        if (acquisitionMade==1)
                            
                            originalkSpaceReconstructed(:,whichLine,whichSlice,volume) = originalOrganizedData(:,posData);
                            %fprintf("\n Row %d of slice %d of diffusion volume %d acquired.\n", whichLine, whichSlice, volume);
                            posData = posData + 1;
                        end
                        
                        
                        globalPos = globalPos + 1;
                    end
                        
                end
                
            end                                              
            
            kSpaceReconstructed = originalkSpaceReconstructed(1:finalReadoutDim,:,:,:);
            

            %kSpaceReconstructed = originalOrganizedData(1:finalReadoutDim,:,round((origPhase2Dim-finalPhase2Dim)/2):(round((origPhase2Dim-finalPhase2Dim)/2)+finalPhase2Dim-1));
            imageZPReconstructed = fftshift(ifft(fftshift(kSpaceReconstructed,1),[],1),1);
            imageZPReconstructed = fftshift(ifft(fftshift(imageZPReconstructed,2),[],2),2);
            imageZPReconstructed = fftshift(ifft(fftshift(imageZPReconstructed,3),[],3),3);


            
            if (showResults)
                disp('Showing results ...');
                
                figure;
                subplot(1,3,1)
                imagesc(abs(squeeze(originalkSpaceReconstructed(:,:,whichSliceZ,whichVolume))).^0.25);
                title({'Original undersampled k-space','(with readout zero-padding by Paravision)','X-Y plan'});

                subplot(1,3,2)
                imshow(abs(squeeze(kSpaceReconstructed(:,:,whichSliceZ,whichVolume))).^0.25,[]);
                title({'Undersampled K-space cropped','(without readout zero-padding by Paravision)','X-Y plan'});

                subplot(1,3,3)
                imshow(abs(squeeze(imageZPReconstructed(:,:,whichSliceZ,whichVolume))),[]);
                title({'Zero-filling image','X-Y plan'});
                
                
                figure;
                subplot(1,3,1)
                imagesc(abs(squeeze(originalkSpaceReconstructed(:,whichSliceY,:,whichVolume))).^0.25);
                title({'Original undersampled k-space','(with readout zero-padding by Paravision)','X-Z plan'});

                subplot(1,3,2)
                imshow(abs(squeeze(kSpaceReconstructed(:,whichSliceY,:,whichVolume))).^0.25,[]);
                title({'Undersampled K-space cropped','(without readout zero-padding by Paravision)','X-Z plan'});

                subplot(1,3,3)
                imshow(abs(squeeze(imageZPReconstructed(:,whichSliceY,:,whichVolume))),[]);
                title({'Zero-filling image','X-Z plan'});
                
                
                figure;
                subplot(1,3,1)
                imagesc(abs(squeeze(originalkSpaceReconstructed(whichSliceX,:,:,whichVolume))).^0.25);
                title({'Original undersampled k-space','(with readout zero-padding by Paravision)','Y-Z plan'});

                subplot(1,3,2)
                imshow(abs(squeeze(kSpaceReconstructed(whichSliceX,:,:,whichVolume))).^0.25,[]);
                title({'Undersampled K-space cropped','(without readout zero-padding by Paravision)','Y-Z plan'});

                subplot(1,3,3)
                imshow(abs(squeeze(imageZPReconstructed(whichSliceX,:,:,whichVolume))),[]);
                title({'Zero-filling image','Y-Z plan'});
            
                csMaskRetrieved = squeeze(abs(kSpaceReconstructed) ~= 0);
                figure;
                imshow(squeeze(csMaskRetrieved(whichSliceX,:,:,whichVolume)),[]);
                title('Compressed Sensing mask retrieved - Y-Z plan - Volume 1');                
                
            else                
                disp('Show results not selected!');
                
            end
            
            
       
        
        otherwise
            disp('No one selected')
            kSpaceReconstructed = NaN; 
            imageZPReconstructed = NaN;
            csMaskRetrieved = NaN;
    end



end