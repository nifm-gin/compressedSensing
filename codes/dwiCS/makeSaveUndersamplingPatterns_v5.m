function makeSaveUndersamplingPatterns_v5(dataInfo,pathToSaveUndersamplingPatterns,saveUndersamplingPatterns)

    % Script for making undersampling patterns --------------------------------


    %%
    % Read acquisition file in order to reconstruct k-space correctly -------------------------------------------------------------------------------------------------------------------------
    %fprintf("\nAttention: The matrix dimensions of all subjects must be the same, given that the same undersampling patterns will be applied!!!\n")
    %acquisitionFilePath = dataInfo.acquisitionParameters(1);
    %[acquisitionTypeID, origReadoutDim, origPhase1Dim, origPhase2Dim, nEchos, dataInfo.undersamplingPatternDimensions(1), dataInfo.undersamplingPatternDimensions(2), dataInfo.undersamplingPatternDimensions(3), dataInfo.undersamplingPatternDimensions(4), dataInfo.B0Number, nDWIs, dataInfo.undersamplingPatternDimensions(5)] = readAcquisitionFile(acquisitionFilePath);

    % -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    nDWIs = dataInfo.undersamplingPatternDimensions(5) - dataInfo.B0Number;

    firstTimeFlag = 1;
    for AF = dataInfo.AF
        for undersamplingPatternFullySampledRegionArea = dataInfo.undersamplingPatternFullySampledRegionArea
            for undersamplingPatternStrategy = dataInfo.undersamplingPatternStrategy
                for undersamplingPatternType= dataInfo.undersamplingPatternType

                    %flagUndersampledB0MonteCarlo = 0;
                    %flagFullySampledB0MonteCarlo = 0;
                    %flagUndersampledB0PoissonUniformDensity = 0;
                    %flagFullySampledB0PoissonUniformDensity = 0;
                    %flagUndersampledB0PoissonVariableDensity = 0;
                    %flagFullySampledB0PoissonVariableDensity = 0;
                    for B0SamplingStrategy=dataInfo.B0SamplingStrategy



                        % Chossing 4D strategy -----------------------------------------------------------------------------------                    
                        if (undersamplingPatternStrategy=="Multi3D" || undersamplingPatternStrategy=="Mono3D")

                            % Just to show which strategy was selected ------------------------------------------------------------
                            if (undersamplingPatternStrategy=="Multi3D")
                            fprintf('Mono 3D mask - Only one common undersampling pattern will be applied for all 3D volumes.\n\n')                        

                            else
                            fprintf('Multi 3D mask - One undersampling pattern for each 3D volume will be applied.\n\n')
                            end
                            % -undersamplingPattern----------------------------------------------------------------------------------------------------

                                % Choosing type of undersampling pattern -----------------------------------------------------------                
                                if (undersamplingPatternType =="MonteCarlo-VariableDensity")

                                    fprintf('Creating a Monte Carlo variable-density undersampling pattern...\n\n')

                                    %[mask4D] = makeMonteCarloUndersamplingPattern(dataInfo.undersamplingPatternDimensions(1),dataInfo.undersamplingPatternDimensions(2),dataInfo.undersamplingPatternDimensions(3),dataInfo.undersamplingPatternDimensions(4),dataInfo.undersamplingPatternDimensions(5))

                                    %if (((B0SamplingStrategy=="ConventionalCS-UndersampledB0") || (B0SamplingStrategy=="KLRCS-UndersampledB0")) && (not(flagUndersampledB0MonteCarlo)))
                                    if (B0SamplingStrategy=="UndersampledB0")
                                        %undersamplingPattern = zerPoisson-UniformDensity-EllipticalScanning"os(dataInfo.undersamplingPatternDimensions(1),dataInfo.undersamplingPatternDimensions(2),dataInfo.undersamplingPatternDimensions(3),dataInfo.undersamplingPatternDimensions(4),dataInfo.undersamplingPatternDimensions(5));
                                        disp("============== Option 1 ==============")
                                        disp(B0SamplingStrategy)
                                        correctedAF = AF + dataInfo.marginAF;
                                        [undersamplingPattern] = makeMonteCarloUndersamplingPattern (dataInfo.undersamplingPatternDimensions(1), dataInfo.undersamplingPatternDimensions(2), dataInfo.undersamplingPatternDimensions(3), dataInfo.undersamplingPatternDimensions(4), 0, dataInfo.undersamplingPatternDimensions(5), correctedAF, dataInfo.undersamplingPatternFullySampledRegionArea, undersamplingPatternStrategy);
                                        maskDiff = sum(sum(sum((abs(undersamplingPattern(:,:,:,1,(dataInfo.B0Number+1)) - undersamplingPattern(:,:,:,1,(dataInfo.undersamplingPatternDimensions(5))))))));
                                        %flagUndersampledB0MonteCarlo = 1;                                                                                                                                                


                                    %elseif (((B0SamplingStrategy=="ConventionalCS-FullySampledB0") || (B0SamplingStrategy=="KLRCS-FullySampledB0")) && (not(flagFullySampledB0MonteCarlo)))
                                    elseif (B0SamplingStrategy=="FullySampledB0")
                                        disp("================= Option 2 ==============")
                                        disp(B0SamplingStrategy)
                                        correctedAF = (nDWIs)/((dataInfo.undersamplingPatternDimensions(5)/(AF + dataInfo.marginAF))-dataInfo.B0Number);
                                        [undersamplingPattern] = makeMonteCarloUndersamplingPattern (dataInfo.undersamplingPatternDimensions(1), dataInfo.undersamplingPatternDimensions(2), dataInfo.undersamplingPatternDimensions(3), dataInfo.undersamplingPatternDimensions(4), dataInfo.B0Number, dataInfo.undersamplingPatternDimensions(5), correctedAF, dataInfo.undersamplingPatternFullySampledRegionArea, undersamplingPatternStrategy);                                    
                                        maskDiff = sum(sum(sum((abs(undersamplingPattern(:,:,:,1,(dataInfo.B0Number+1)) - undersamplingPattern(:,:,:,1,(dataInfo.undersamplingPatternDimensions(5))))))));
                                        %flagFullySampledB0MonteCarlo = 1;                                                                        

                                    else
                                        fprintf('CS type option chosen does not match any case!!! Please, choose one among the available options!!!\n\n')
                                    end



                                elseif (undersamplingPatternType =="Poisson-UniformDensity-EllipticalScanning")
                                    %if (((B0SamplingStrategy=="ConventionalCS-"Poisson-UniformDensity"UndersampledB0") || (B0SamplingStrategy=="KLRCS-UndersampledB0")) && (not(flagUndersampledB0PoissonUniformDensity)))
                                    if (B0SamplingStrategy=="UndersampledB0")
                                        %undersamplingPattern = zeros(dataInfo.undersamplingPatternDimensions(1),dataInfo.undersamplingPatternDimensions(2),dataInfo.undersamplingPatternDimensions(3),dataInfo.undersamplingPatternDimensions(4),dataInfo.undersamplingPatternDimensions(5));
                                        disp("============== Option 3 ==============")
                                        disp(B0SamplingStrategy)
                                        %undersamplingPattern = 0.25.*ones(dataInfo.undersamplingPatternDimensions(1),dataInfo.undersamplingPatternDimensions(2),dataInfo.undersamplingPatternDimensions(3),dataInfo.undersamplingPatternDimensions(4),dataInfo.undersamplingPatternDimensions(5));                                    
                                        correctedAF = AF + dataInfo.marginAF;
                                        [undersamplingPattern] = makePoissonDiskUndersamplingPattern (dataInfo.undersamplingPatternDimensions(1), dataInfo.undersamplingPatternDimensions(2), dataInfo.undersamplingPatternDimensions(3), dataInfo.undersamplingPatternDimensions(4), 0, dataInfo.undersamplingPatternDimensions(5), correctedAF, dataInfo.undersamplingPatternFullySampledRegionArea, 0, 1, undersamplingPatternStrategy, 0);
                                        maskDiff = sum(sum(sum((abs(undersamplingPattern(:,:,:,1,(dataInfo.B0Number+1)) - undersamplingPattern(:,:,:,1,(dataInfo.undersamplingPatternDimensions(5))))))));
                                        %flagUndersampledB0PoissonUniformDensity = 1;



                                    %elseif (((B0SamplingStrategy=="ConventionalCS-FullySampledB0") || (B0SamplingStrategy=="KLRCS-FullySampledB0")) && (not(flagFullySampledB0PoissonUniformDensity)))
                                    elseif (B0SamplingStrategy=="FullySampledB0")
                                        disp("================= Option 4 ==============")
                                        disp(B0SamplingStrategy)
                                        %undersamplingPattern = 0.5.*ones(dataInfo.undersamplingPatternDimensions(1),dataInfo.undersamplingPatternDimensions(2),dataInfo.undersamplingPatternDimensions(3),dataInfo.undersamplingPatternDimensions(4),dataInfo.undersamplingPatternDimensions(5));
                                        correctedAF = (nDWIs)/((dataInfo.undersamplingPatternDimensions(5)/(AF + dataInfo.marginAF))-dataInfo.B0Number);
                                        [undersamplingPattern] = makePoissonDiskUndersamplingPattern (dataInfo.undersamplingPatternDimensions(1), dataInfo.undersamplingPatternDimensions(2), dataInfo.undersamplingPatternDimensions(3), dataInfo.undersamplingPatternDimensions(4), dataInfo.B0Number, dataInfo.undersamplingPatternDimensions(5), correctedAF, dataInfo.undersamplingPatternFullySampledRegionArea, 0, 1, undersamplingPatternStrategy, 0);
                                        maskDiff = sum(sum(sum((abs(undersamplingPattern(:,:,:,1,(dataInfo.B0Number+1)) - undersamplingPattern(:,:,:,1,(dataInfo.undersamplingPatternDimensions(5))))))));
                                        %flagFullySampledB0PoissonUniformDensity = 1;

                                    else
                                        fprintf('CS type option chosen does not match any case!!! Please, choose one among the available options!!!\n\n')
                                    end




                                elseif (undersamplingPatternType =="Poisson-VariableDensity-EllipticalScanning")
                                    %if (((B0SamplingStrategy=="ConventionalCS-UndersampledB0") || (B0SamplingStrategy=="KLRCS-UndersampledB0")) && (not(flagUndersampledB0PoissonVariableDensity)))
                                    if (B0SamplingStrategy=="UndersampledB0")
                                        %undersamplingPattern = zeros(dataInfo.undersamplingPatternDimensions(1),dataInfo.undersamplingPatternDimensions(2),dataInfo.undersamplingPatternDimensions(3),dataInfo.undersamplingPatternDimensions(4),dataInfo.undersamplingPatternDimensions(5));
                                        disp("============== Option 5 ==============")
                                        disp(B0SamplingStrategy)
                                        %undersamplingPattern = 0.75.*ones(dataInfo.undersamplingPatternDimensions(1),dataInfo.undersamplingPatternDimensions(2),dataInfo.undersamplingPatternDimensions(3),dataInfo.undersamplingPatternDimensions(4),dataInfo.undersamplingPatternDimensions(5));
                                        correctedAF = AF + dataInfo.marginAF;
                                        [undersamplingPattern] = makePoissonDiskUndersamplingPattern (dataInfo.undersamplingPatternDimensions(1), dataInfo.undersamplingPatternDimensions(2), dataInfo.undersamplingPatternDimensions(3), dataInfo.undersamplingPatternDimensions(4), 0, dataInfo.undersamplingPatternDimensions(5), correctedAF, dataInfo.undersamplingPatternFullySampledRegionArea, 1, 1, undersamplingPatternStrategy, 0);
                                        maskDiff = sum(sum(sum((abs(undersamplingPattern(:,:,:,1,(dataInfo.B0Number+1)) - undersamplingPattern(:,:,:,1,(dataInfo.undersamplingPatternDimensions(5))))))));
                                        %flagUndersampledB0PoissonVariableDensity = 1;



                                    %elseif (((B0SamplingStrategy=="ConventionalCS-FullySampledB0") || (B0SamplingStrategy=="KLRCS-FullySampledB0")) && (not(flagFullySampledB0PoissonVariableDensity)))
                                    elseif (B0SamplingStrategy=="FullySampledB0")
                                        disp("================= Option 6 ==============")
                                        disp(B0SamplingStrategy)
                                        %undersamplingPattern = 1.*ones(dataInfo.undersamplingPatternDimensions(1),dataInfo.undersamplingPatternDimensions(2),dataInfo.undersamplingPatternDimensions(3),dataInfo.undersamplingPatternDimensions(4),dataInfo.undersamplingPatternDimensions(5));
                                        correctedAF = (nDWIs)/((dataInfo.undersamplingPatternDimensions(5)/(AF + dataInfo.marginAF))-dataInfo.B0Number);
                                        [undersamplingPattern] = makePoissonDiskUndersamplingPattern (dataInfo.undersamplingPatternDimensions(1), dataInfo.undersamplingPatternDimensions(2), dataInfo.undersamplingPatternDimensions(3), dataInfo.undersamplingPatternDimensions(4), dataInfo.B0Number, dataInfo.undersamplingPatternDimensions(5), correctedAF, dataInfo.undersamplingPatternFullySampledRegionArea, 1, 1, undersamplingPatternStrategy, 0);
                                        maskDiff = sum(sum(sum((abs(undersamplingPattern(:,:,:,1,(dataInfo.B0Number+1)) - undersamplingPattern(:,:,:,1,(dataInfo.undersamplingPatternDimensions(5))))))));
                                        %flagFullySampledB0PoissonVariableDensity = 1;

                                    else
                                        fprintf('CS type option chosen does not match any case!!! Please, choose one among the available options!!!\n\n')
                                    end



                                elseif (undersamplingPatternType =="Poisson-UniformDensity")
                                    %if (((B0SamplingStrategy=="ConventionalCS-UndersampledB0") || (B0SamplingStrategy=="KLRCS-UndersampledB0")) && (not(flagUndersampledB0PoissonUniformDensity)))
                                    if (B0SamplingStrategy=="UndersampledB0")
                                        %undersamplingPattern = zeros(dataInfo.undersamplingPatternDimensions(1),dataInfo.undersamplingPatternDimensions(2),dataInfo.undersamplingPatternDimensions(3),dataInfo.undersamplingPatternDimensions(4),dataInfo.undersamplingPatternDimensions(5));
                                        disp("============== Option 3 ==============")
                                        disp(B0SamplingStrategy)
                                        %undersamplingPattern = 0.25.*ones(dataInfo.undersamplingPatternDimensions(1),dataInfo.undersamplingPatternDimensions(2),dataInfo.undersamplingPatternDimensions(3),dataInfo.undersamplingPatternDimensions(4),dataInfo.undersamplingPatternDimensions(5));                                    
                                        correctedAF = AF + dataInfo.marginAF;
                                        [undersamplingPattern] = makePoissonDiskUndersamplingPattern (dataInfo.undersamplingPatternDimensions(1), dataInfo.undersamplingPatternDimensions(2), dataInfo.undersamplingPatternDimensions(3), dataInfo.undersamplingPatternDimensions(4), 0, dataInfo.undersamplingPatternDimensions(5), correctedAF, dataInfo.undersamplingPatternFullySampledRegionArea, 0, 0, undersamplingPatternStrategy, 0);
                                        maskDiff = sum(sum(sum((abs(undersamplingPattern(:,:,:,1,(dataInfo.B0Number+1)) - undersamplingPattern(:,:,:,1,(dataInfo.undersamplingPatternDimensions(5))))))));
                                        %flagUndersampledB0PoissonUniformDensity = 1;



                                    %elseif (((B0SamplingStrategy=="ConventionalCS-FullySampledB0") || (B0SamplingStrategy=="KLRCS-FullySampledB0")) && (not(flagFullySampledB0PoissonUniformDensity)))
                                    elseif (B0SamplingStrategy=="FullySampledB0")
                                        disp("================= Option 4 ==============")
                                        disp(B0SamplingStrategy)
                                        %undersamplingPattern = 0.5.*ones(dataInfo.undersamplingPatternDimensions(1),dataInfo.undersamplingPatternDimensions(2),dataInfo.undersamplingPatternDimensions(3),dataInfo.undersamplingPatternDimensions(4),dataInfo.undersamplingPatternDimensions(5));
                                        correctedAF = (nDWIs)/((dataInfo.undersamplingPatternDimensions(5)/(AF + dataInfo.marginAF))-dataInfo.B0Number);
                                        [undersamplingPattern] = makePoissonDiskUndersamplingPattern (dataInfo.undersamplingPatternDimensions(1), dataInfo.undersamplingPatternDimensions(2), dataInfo.undersamplingPatternDimensions(3), dataInfo.undersamplingPatternDimensions(4), dataInfo.B0Number, dataInfo.undersamplingPatternDimensions(5), correctedAF, dataInfo.undersamplingPatternFullySampledRegionArea, 0, 0, undersamplingPatternStrategy, 0);
                                        maskDiff = sum(sum(sum((abs(undersamplingPattern(:,:,:,1,(dataInfo.B0Number+1)) - undersamplingPattern(:,:,:,1,(dataInfo.undersamplingPatternDimensions(5))))))));
                                        %flagFullySampledB0PoissonUniformDensity = 1;

                                    else
                                        fprintf('CS type option chosen does not match any case!!! Please, choose one among the available options!!!\n\n')
                                    end




                                elseif (undersamplingPatternType =="Poisson-VariableDensity")
                                    %if (((B0SamplingStrategy=="ConventionalCS-UndersampledB0") || (B0SamplingStrategy=="KLRCS-UndersampledB0")) && (not(flagUndersampledB0PoissonVariableDensity)))
                                    if (B0SamplingStrategy=="UndersampledB0")
                                        %undersamplingPattern = zeros(dataInfo.undersamplingPatternDimensions(1),dataInfo.undersamplingPatternDimensions(2),dataInfo.undersamplingPatternDimensions(3),dataInfo.undersamplingPatternDimensions(4),dataInfo.undersamplingPatternDimensions(5));
                                        disp("============== Option 5 ==============")
                                        disp(B0SamplingStrategy)
                                        %undersamplingPattern = 0.75.*ones(dataInfo.undersamplingPatternDimensions(1),dataInfo.undersamplingPatternDimensions(2),dataInfo.undersamplingPatternDimensions(3),dataInfo.undersamplingPatternDimensions(4),dataInfo.undersamplingPatternDimensions(5));
                                        correctedAF = AF + dataInfo.marginAF;
                                        [undersamplingPattern] = makePoissonDiskUndersamplingPattern (dataInfo.undersamplingPatternDimensions(1), dataInfo.undersamplingPatternDimensions(2), dataInfo.undersamplingPatternDimensions(3), dataInfo.undersamplingPatternDimensions(4), 0, dataInfo.undersamplingPatternDimensions(5), correctedAF, dataInfo.undersamplingPatternFullySampledRegionArea, 1, 0, undersamplingPatternStrategy, 0);
                                        maskDiff = sum(sum(sum((abs(undersamplingPattern(:,:,:,1,(dataInfo.B0Number+1)) - undersamplingPattern(:,:,:,1,(dataInfo.undersamplingPatternDimensions(5))))))));
                                        %flagUndersampledB0PoissonVariableDensity = 1;



                                    %elseif (((B0SamplingStrategy=="ConventionalCS-FullySampledB0") || (B0SamplingStrategy=="KLRCS-FullySampledB0")) && (not(flagFullySampledB0PoissonVariableDensity)))
                                    elseif (B0SamplingStrategy=="FullySampledB0")
                                        disp("================= Option 6 ==============")
                                        disp(B0SamplingStrategy)
                                        %undersamplingPattern = 1.*ones(dataInfo.undersamplingPatternDimensions(1),dataInfo.undersamplingPatternDimensions(2),dataInfo.undersamplingPatternDimensions(3),dataInfo.undersamplingPatternDimensions(4),dataInfo.undersamplingPatternDimensions(5));
                                        correctedAF = (nDWIs)/((dataInfo.undersamplingPatternDimensions(5)/(AF + dataInfo.marginAF))-dataInfo.B0Number);
                                        [undersamplingPattern] = makePoissonDiskUndersamplingPattern (dataInfo.undersamplingPatternDimensions(1), dataInfo.undersamplingPatternDimensions(2), dataInfo.undersamplingPatternDimensions(3), dataInfo.undersamplingPatternDimensions(4), dataInfo.B0Number, dataInfo.undersamplingPatternDimensions(5), correctedAF, dataInfo.undersamplingPatternFullySampledRegionArea, 1, 0, undersamplingPatternStrategy, 0);
                                        maskDiff = sum(sum(sum((abs(undersamplingPattern(:,:,:,1,(dataInfo.B0Number+1)) - undersamplingPattern(:,:,:,1,(dataInfo.undersamplingPatternDimensions(5))))))));
                                        %flagFullySampledB0PoissonVariableDensity = 1;

                                    else
                                        fprintf('CS type option chosen does not match any case!!! Please, choose one among the available options!!!\n\n')
                                    end


                                else

                                    fprintf('Option chosen does not match any case!!! Please, choose one among the available options!!!\n\n')


                                end                                                               
                                % --------------------------------------------------------------------------------------------------------                                          


                        else

                            fprintf('Option chosen does not match any case!!! Please, choose one among the available options!!!\n\n')

                        end                

                        % --------------------------------------------------------------------------------------------------------
                        % --------------------------------------------------------------------------------------------------------

                        actual4DAF = (size(undersamplingPattern,1)*size(undersamplingPattern,2)*size(undersamplingPattern,3)*size(undersamplingPattern,4)*size(undersamplingPattern,5) / sum(undersamplingPattern(:)));                    
                        fprintf('\nUndersampling pattern made!!!\n\n')

                        if (dataInfo.showResults == 1)
                            if (dataInfo.B0Number)
                                figure('units','normalized','outerposition',[0 0 1 1])
                                subplot(4,3,1)
                                imshow(squeeze(undersamplingPattern(dataInfo.slicesToShow(1),:,:,1,1)),[])
                                title(strcat("YZ view - slice X=",num2str(dataInfo.slicesToShow(1)),", volume=",num2str(1)))

                                subplot(4,3,2)
                                imshow(squeeze(undersamplingPattern(:,dataInfo.slicesToShow(2),:,1,1)),[])
                                title(strcat("XZ view - slice Y=",num2str(dataInfo.slicesToShow(2)),", volume=",num2str(1)))

                                subplot(4,3,3)
                                imshow(squeeze(undersamplingPattern(:,:,dataInfo.slicesToShow(3),1,1)),[])
                                title(strcat("XY view - slice Z=",num2str(dataInfo.slicesToShow(3)),", volume=",num2str(1)))


                                subplot(4,3,4)
                                imshow(squeeze(undersamplingPattern(dataInfo.slicesToShow(1),:,:,1,dataInfo.B0Number)),[])
                                title(strcat("YZ view - slice X=",num2str(dataInfo.slicesToShow(1)),", volume=",num2str(dataInfo.B0Number)))

                                subplot(4,3,5)
                                imshow(squeeze(undersamplingPattern(:,dataInfo.slicesToShow(2),:,1,dataInfo.B0Number)),[])
                                title(strcat("XZ view - slice Y=",num2str(dataInfo.slicesToShow(2)),", volume=",num2str(dataInfo.B0Number)))

                                subplot(4,3,6)
                                imshow(squeeze(undersamplingPattern(:,:,dataInfo.slicesToShow(3),1,dataInfo.B0Number)),[])
                                title(strcat("XY view - slice Z=",num2str(dataInfo.slicesToShow(3)),", volume=",num2str(dataInfo.B0Number)))


                                subplot(4,3,7)
                                imshow(squeeze(undersamplingPattern(dataInfo.slicesToShow(1),:,:,1,(dataInfo.B0Number+1))),[])
                                title(strcat("YZ view - slice X=",num2str(dataInfo.slicesToShow(1)),", volume=",num2str(dataInfo.B0Number+1)))

                                subplot(4,3,8)
                                imshow(squeeze(undersamplingPattern(:,dataInfo.slicesToShow(2),:,1,(dataInfo.B0Number+1))),[])
                                title(strcat("XZ view - slice Y=",num2str(dataInfo.slicesToShow(2)),", volume=",num2str(dataInfo.B0Number+1)))

                                subplot(4,3,9)                        
                                imshow(squeeze(undersamplingPattern(:,:,dataInfo.slicesToShow(3),1,(dataInfo.B0Number+1))),[])
                                title(strcat("XY view - slice Z=",num2str(dataInfo.slicesToShow(3)),", volume=",num2str(dataInfo.B0Number+1)))


                                subplot(4,3,10)
                                imshow(squeeze(undersamplingPattern(dataInfo.slicesToShow(1),:,:,1,dataInfo.undersamplingPatternDimensions(5))),[])
                                title(strcat("YZ view - slice X=",num2str(dataInfo.slicesToShow(1)),", volume=",num2str(dataInfo.undersamplingPatternDimensions(5))))

                                subplot(4,3,11)
                                imshow(squeeze(undersamplingPattern(:,dataInfo.slicesToShow(2),:,1,dataInfo.undersamplingPatternDimensions(5))),[])
                                title(strcat("XZ view - slice Y=",num2str(dataInfo.slicesToShow(2)),", volume=",num2str(dataInfo.undersamplingPatternDimensions(5))))

                                subplot(4,3,12)
                                imshow(squeeze(undersamplingPattern(:,:,dataInfo.slicesToShow(3),1,dataInfo.undersamplingPatternDimensions(5))),[])
                                title(strcat("XY view - slice Z=",num2str(dataInfo.slicesToShow(3)),", volume=",num2str(dataInfo.undersamplingPatternDimensions(5))))

                            else
                                figure('units','normalized','outerposition',[0 0 1 1])
                                subplot(1,3,1)
                                imshow(squeeze(undersamplingPattern(dataInfo.slicesToShow(1),:,:,1,1)),[])
                                title(strcat("YZ view - slice X=",num2str(dataInfo.slicesToShow(1)),", volume=",num2str(1)))

                                subplot(1,3,2)
                                imshow(squeeze(undersamplingPattern(:,dataInfo.slicesToShow(2),:,1,1)),[])
                                title(strcat("XZ view - slice Y=",num2str(dataInfo.slicesToShow(2)),", volume=",num2str(1)))

                                subplot(1,3,3)
                                imshow(squeeze(undersamplingPattern(:,:,dataInfo.slicesToShow(3),1,1)),[])
                                title(strcat("XY view - slice Z=",num2str(dataInfo.slicesToShow(3)),", volume=",num2str(1)))

                            end
                            sgtitle(strcat(undersamplingPatternType," - ",B0SamplingStrategy," - Nominal AF=",num2str(AF)," - Actual overall AF=",num2str(actual4DAF)," - DW AF=",num2str(correctedAF)," - FS area: ",num2str(undersamplingPatternFullySampledRegionArea)," - ",undersamplingPatternStrategy," - maskDiff: ",num2str(maskDiff)))                    
                            pause(1)



                        end

                        if (saveUndersamplingPatterns)
                            if (firstTimeFlag)                                                                                                
                                timeNow = datevec(now);
                                txtFilename = strcat(num2str(timeNow(1)),num2str(timeNow(2),'%.2d'),num2str(timeNow(3),'%.2d'),'_',num2str(timeNow(4),'%.2d'),num2str(timeNow(5),'%.2d'),num2str(timeNow(6),'%02.0f'),'_',"5D-UndersamplingPatterns_",num2str(size(undersamplingPattern,1)),"x",num2str(size(undersamplingPattern,2)),"x",num2str(size(undersamplingPattern,3)),"x",num2str(size(undersamplingPattern,4)),"x",num2str(size(undersamplingPattern,5)),"_",num2str(dataInfo.B0Number),"-B0s.txt");

                                fileID = fopen(strcat(pathToSaveUndersamplingPatterns,txtFilename),'w');
                                %totalNumberOfUndersamplingPatterns = size(dataInfo.AF,2)*size(dataInfo.undersamplingPatternFullySampledRegionArea,2)*size(dataInfo.undersamplingPatternStrategy,2)*size(dataInfo.undersamplingPatternType,2)*size(dataInfo.B0SamplingStrategy,2);
                                %fprintf(fileID,'%d\r\n',totalNumberOfUndersamplingPatterns); 
                                fclose(fileID);
                            end

                            fprintf('\nSaving undersampling pattern as .mat and image containing some results...')

                            filename = strcat(num2str(timeNow(1)),num2str(timeNow(2),'%.2d'),num2str(timeNow(3),'%.2d'),'_',num2str(timeNow(4),'%.2d'),num2str(timeNow(5),'%.2d'),num2str(timeNow(6),'%02.0f'),'_',"5D-UndersamplingPattern_",num2str(size(undersamplingPattern,1)),"x",num2str(size(undersamplingPattern,2)),"x",num2str(size(undersamplingPattern,3)),"x",num2str(size(undersamplingPattern,4)),"x",num2str(size(undersamplingPattern,5)),"_",num2str(dataInfo.B0Number),"-B0s_NominalAF-",num2str(AF),"_OverallAF-",num2str(actual4DAF),"_DWIAF-",num2str(correctedAF),"_",undersamplingPatternType,"_",undersamplingPatternStrategy,"_FSArea-",num2str(undersamplingPatternFullySampledRegionArea),"_",B0SamplingStrategy);

                            saveas(gcf,strcat(pathToSaveUndersamplingPatterns,filename, ".jpg"));

                            nB0s = dataInfo.B0Number;
                            save(strcat(pathToSaveUndersamplingPatterns,filename,".mat"),'undersamplingPattern','nB0s','AF','actual4DAF','correctedAF','maskDiff','undersamplingPatternFullySampledRegionArea','undersamplingPatternStrategy','undersamplingPatternType','B0SamplingStrategy');

                            fileID = fopen(strcat(pathToSaveUndersamplingPatterns,txtFilename),'a');
                            fprintf(fileID,'%s\r\n',strcat(pathToSaveUndersamplingPatterns,filename,".mat"));
                            fclose(fileID);

                            [undersamplingPatternTxtDocPath] = convertUndersamplingPatternTxtFile(strcat(pathToSaveUndersamplingPatterns,filename,".mat"), 1, pathToSaveUndersamplingPatterns);

                            [undersamplingPatternTxtDocPathBinaryArray] = convertUndersamplingPatternTxtFileBinaryArray(strcat(pathToSaveUndersamplingPatterns,filename,".mat"), pathToSaveUndersamplingPatterns);

                            fprintf("\nUndersampling pattern and image saved with success. \n########################################################################################### \n")
                        end
                        firstTimeFlag = 0;                    
                    end
                end
            end
        end
    end

end
% -------------------------------------------------------------------------