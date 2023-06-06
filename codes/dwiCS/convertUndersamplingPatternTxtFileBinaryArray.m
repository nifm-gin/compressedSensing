function [undersamplingPatternTxtDocPath] = convertUndersamplingPatternTxtFileBinaryArray (undersamplingPatternPath, savePath)
    
    
    fprintf("\nLoading undersampling pattern file...\n");
    load(undersamplingPatternPath);
    %undersamplingPattern = undersamplingPattern5D; %Temp - if "Unrecognized function or variable 'undersamplingPattern'"
    
    fprintf("Undersampling pattern file loaded...\n");
    
    %{
    if (nSlices=="All")
        fprintf("All slices of the 3D volume will be taken.\n");        
        sizePhase2 = size(undersamplingPattern,3);
    
    else
        fprintf("Only %d will be taken.\n",nSlices);
    end
    %}
    
    
    if (undersamplingPatternType=="MonteCarlo-VariableDensity" || undersamplingPatternType=="Poisson-VariableDensity" || undersamplingPatternType=="Poisson-VariableDensity-EllipticalScanning")
        variableDensity = 1;
    else
        variableDensity = 0;
    end
    
    
    if (undersamplingPatternType=="Poisson-UniformDensity-EllipticalScanning" || undersamplingPatternType=="Poisson-VariableDensity-EllipticalScanning")
        ellipse = 1;
    else
        ellipse = 0;
    end
    
    if (undersamplingPatternType=="MonteCarlo-VariableDensity")
        centerDim = 2*round(sqrt(undersamplingPatternFullySampledRegionArea/pi));
    else
        centerDim = round(sqrt(undersamplingPatternFullySampledRegionArea));
    end
    
    whichCase = 3; %3D Compressed Sensing
    seed = 1;
    sizeReadout = size(undersamplingPattern,1);
    sizePhase1 = size(undersamplingPattern,2);
    sizePhase2 = size(undersamplingPattern,3);
    nVolumes = size(undersamplingPattern,5);
    actualAF = (size(undersamplingPattern,1)*size(undersamplingPattern,2)*size(undersamplingPattern,3)*size(undersamplingPattern,4)*size(undersamplingPattern,5)) / (sum(undersamplingPattern(:)));
    numberOfLinesAcquired = sum(sum(sum(sum(undersamplingPattern(1,:,:,1,:)))));
    
    nDiffusionDirections = nVolumes - nB0s;
    
    fprintf("An undersampling pattern for each 3D volume will be saved in the txt file (multi mask).\n");

    timeNow = datevec(now);
    undersamplingPatternTxtDocPath = (strcat(savePath,num2str(timeNow(1)),num2str(timeNow(2),'%.2d'),num2str(timeNow(3),'%.2d'),'_',num2str(timeNow(4),'%.2d'),num2str(timeNow(5),'%.2d'),num2str(timeNow(6),'%02.0f'),'_UndersamplingPattern_',num2str(size(undersamplingPattern,1)),'x',num2str(size(undersamplingPattern,2)),'x',num2str(size(undersamplingPattern,3)),'_AF-',num2str(AF),'_overallAF-',num2str(actual4DAF),'_DWIAF-',num2str(correctedAF),'_nVol-',num2str(size(undersamplingPattern,5)),'_nB0-',num2str(nB0s),'_fsCenterArea-',num2str(undersamplingPatternFullySampledRegionArea),'_',undersamplingPatternType,'_',undersamplingPatternStrategy,'_',B0SamplingStrategy,'_BinAcqGuideArray.txt'));

    disp(undersamplingPatternTxtDocPath);
    %fprintf('%d\n%d\n%d\n%d\n%d\n%f\n%d\n%d\n%d\n%d\n%d\n', whichCase,sizeReadout,sizePhase1,sizePhase2,nVolumes,actualAF,centerDim,ellipse,variableDensity,seed,numberOfLinesAcquired);        
    fid = fopen(undersamplingPatternTxtDocPath, 'wt');
    fprintf(fid, '%d\n%d\n%d\n%d\n%d\n%d\n%f\n%d\n%d\n%d\n%d\n%d\n\n', whichCase,sizeReadout,sizePhase1,sizePhase2,nB0s,nDiffusionDirections,actualAF,centerDim,ellipse,variableDensity,seed,numberOfLinesAcquired);                        
    
    
    binaryGuideArrayForAcquisitions = reshape((permute(squeeze(undersamplingPattern(1,:,:,1,:)), [3 1 2])), [size(undersamplingPattern,5)*size(undersamplingPattern,2)*size(undersamplingPattern,3),1]);       
    
    
    fprintf(fid, '%d\n', binaryGuideArrayForAcquisitions);                
    
    %{
    howManyLines = 0;
    for volume=1:nVolumes
        for slice=1:sizePhase2            
            oneDimMaskBinVector = squeeze(undersamplingPattern(1,:,slice,1,volume));
            %figure,plot(oneDimMaskBinVector,'o'),title('Mask vector');

            acquiredLineIntVector = find(oneDimMaskBinVector);
            %acquiredLineIntMatrix{slice}=(acquiredLineIntVector');

            % Redefine mask according to line chosen
            %poissonMask = zeros(sizeReadout,sizePhase1);    
            %for k = 1:sizePhase1, poissonMask(k,:) = oneDimMaskBinVector ; end


            fprintf(fid, '\n\n%d\n%d\n',slice,size(acquiredLineIntVector,2));                

            if (acquiredLineIntVector)
                fprintf(fid, '\n%d', acquiredLineIntVector);
                fprintf(fid, '\n');                       
                howManyLines = howManyLines + 1;
            %else
                %fprintf(fid, '0\n');  
            end


        end                        
    end
    
    %}
    
    fclose(fid);

    
    fprintf("\nDone!!! Undersampling pattern saved as txt doc in %s\n------------------------------------------------------------------------------------------------------- \n", savePath);
        
            
    
  
    
    
    
end