function [poisson2DUndersamplingPattern] = poisson2DMask(sizePhase1, sizePhase2, AF, centerSquareArea, variableDensity, ellipse, seed, showResults)

    centersquareSize = round(sqrt(centerSquareArea));
    phase1AF = 1;    
    deviationAF = 0.005;
    
    flagStart=0;
    actualAF = 1;
    attemptAFVector = [1:0.0001:AF];
    
    %while ((actualAF < (AF+deviationAF)) || (actualAF > (AF-deviationAF)))
    index=[];
    from=1;
    to=length(attemptAFVector);
    fprintf(strcat("\n=============================================== \nCreating 2D undersampling pattern for AF = ", num2str(AF) ,"\n===============================================\n"))
    while from<=to
        fprintf("-------------------------------------------------\n")
    %while ((actualAF < (AF-deviationAF)) || (actualAF > (AF+deviationAF)))        
        if (flagStart==1)
            mid = round((from + to)/2);  
            from
            to
            mid
        end
        %diff = attemptAFVector(mid) - AF;
    
    
        if (flagStart==0)
            phase1AF = 1
        else
            phase1AF = attemptAFVector(mid)
        end
        
    
        phase2AF = phase1AF; % Same AF for both axis
        
        commandString = 'poisson'; % String for entering command                                            

        commandString = strcat(commandString,' -Y',string(sizePhase1),' -Z',string(sizePhase2),' -y',string(phase1AF),' -z',string(phase2AF),' -C',string(centersquareSize));            

        % Adding more options to mask --------------------
        if (variableDensity)
            commandString = strcat(commandString,' -v');
        end

        if (ellipse)
            commandString = strcat(commandString,' -e');
        end
        % ------------------------------------------------

        % Adding seed --------------------------------------------
        commandString = strcat(commandString,' -s',string(seed));
        % ---------------------------------------------------------


        % Convert string to char -----------
        commandString = char(commandString);
        disp(commandString);
        % ----------------------------------


        % Make mask with string command -------------------------
        poisson2DUndersamplingPattern = squeeze(bart(commandString));    
        actualAF = (size(poisson2DUndersamplingPattern,1)*size(poisson2DUndersamplingPattern,2))/sum(poisson2DUndersamplingPattern(:))

            
        if (showResults)
            figure('units','normalized','outerposition',[0 0 1 1]);
            imshow(squeeze(poisson2DUndersamplingPattern),[]);
            title(strcat("Original mask (2D undersampling) - Actual AF=",num2str(actualAF)));
        end
        

        % -------------------------------------------------------  
        
        
    
        
        if (flagStart==1)
            diff = actualAF - AF
            fprintf("-------------------------------------------------\n")



            if ((abs(diff) < (deviationAF)))
                index=mid;
                fprintf("Undersampling pattern matching chosen AF found!!!\n")
                fprintf("==============================================================================\n")
                return
            elseif diff<0   % x(mid) < sval
                from=mid+1;
                fprintf("Undersampling pattern is not matching chosen AF yet!!!\n")
            else              % x(mid) > sval
                to=mid-1;	
                fprintf("Undersampling pattern is not matching chosen AF yet!!!\n")
            end
        end
        
        
        
        if (flagStart==0 && (actualAF > AF))
            fprintf("It is not possible to obtain an undersampling pattern with the specified AF (AF too small)!!!\n")
            flagStart=1;
            return
        elseif (flagStart==0 && (actualAF < AF))
            fprintf("Don't worry. We're gonna obtain un undersampling pattern!!!\n")
            flagStart=1;
        end
    end
    

fprintf("==============================================================================\n")        
        
    
end