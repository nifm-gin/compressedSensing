function [undersamplingPattern4D] = makePoissonDiskUndersamplingPattern (readoutDim, phase1Dim, sliceSelectionDim, nChannels, nFullySampledB0s, nTotalVolumes, AF, fullySampledArea, variableDensity, ellipse, undersamplingPatternStrategy, showResults)        

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
                SamplingMask(:,:,f) =  poisson2DMask(phase1Dim, sliceSelectionDim, AF, fullySampledArea, variableDensity, ellipse, f, showResults);
                %figure;
                %imshow(abs(squeeze(SamplingMask(:,:,f))),[]);
                
            elseif (undersamplingPatternStrategy=="Mono3D")
                % Make only one 3D undersampling pattern ------------------------------------------------------------------
                if (mono3Dflag==0)
                    SamplingMask(:,:,f) = poisson2DMask(phase1Dim, sliceSelectionDim, AF, fullySampledArea, variableDensity, ellipse, f, showResults);
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