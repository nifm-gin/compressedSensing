function [undersamplingPattern4D] = makeMonteCarloUndersamplingPattern (readoutDim, phase1Dim, sliceSelectionDim, nChannels, nFullySampledB0s, nTotalVolumes, AF, fullySampledArea, undersamplingPatternStrategy)        

fullySampledRegionRadius = round(sqrt(fullySampledArea/pi));
mono3Dflag=0;

    fullySampledB0Count = nFullySampledB0s;
    for f = 1:nTotalVolumes
        
        if (fullySampledB0Count > 0)
            
            SamplingMask(:,:,f) = ones(phase1Dim,sliceSelectionDim);  
            fullySampledB0Count = fullySampledB0Count - 1;
            
            %figure;
            %imshow(abs(squeeze(SamplingMask(:,:,f))),[]);
        
        elseif (fullySampledB0Count == 0)
            
            
            if (undersamplingPatternStrategy=="Multi3D")
                SamplingMask(:,:,f) = rand_sampleMask4k_ori(phase1Dim,sliceSelectionDim,AF,1,fullySampledRegionRadius);
                %figure;
                %imshow(abs(squeeze(SamplingMask(:,:,f))),[]);
                
            elseif (undersamplingPatternStrategy=="Mono3D")
                % Make only one 3D undersampling pattern ------------------------------------------------------------------
                if (mono3Dflag==0)
                    SamplingMask(:,:,f) = rand_sampleMask4k_ori(phase1Dim,sliceSelectionDim,AF,1,fullySampledRegionRadius);
                    mono3Dflag = 1;
                else
                    SamplingMask(:,:,f) = SamplingMask(:,:,(f-1));
                end
                % ---------------------------------------------------------------------------------------------------------
            
            else
                fprintf("Invalid option!!! Plase choose either Multi3D or Mono3D undersampling strategy!!!")
                return
            end
            
        else
            fprintf("Negative B0 quantity!!! Please correct your choice!!!\n")
            return
        end                            
    end
    %undersamplingPattern4D = permute(repmat(SamplingMask,[1,1,1,finalReadoutDim,nChannels]),[1,2,4,5,3]);
    undersamplingPattern4D = permute(repmat(SamplingMask,[1,1,1,readoutDim,nChannels]),[4,1,2,5,3]);

end