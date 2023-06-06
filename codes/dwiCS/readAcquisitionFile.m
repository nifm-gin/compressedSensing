function [acquisitionTypeID, origReadoutDim, origPhase1Dim, origPhase2Dim, nEchos, finalReadoutDim, finalPhase1Dim, finalPhase2Dim, nChannels, nB0s, nDWI, nVolumes] = readAcquisitionFile(txtDocPath)

    fileID = fopen(txtDocPath,'r');
    allData = fscanf(fileID,'%f');
    fclose(fileID);
    
    acquisitionTypeID = allData(1);
    origReadoutDim= allData(2);
    origPhase1Dim= allData(3);
    origPhase2Dim = allData(4); 
    nEchos = allData(5);
    finalReadoutDim = allData(6);
    finalPhase1Dim = allData(7);
    finalPhase2Dim = allData(8);
    nChannels = allData(9);
    nB0s = allData(10);              
    nDWI = allData(11);
    nVolumes = nB0s + nDWI;
    %formatSpec = '%C';
    %T = readtable(txtDocPath,'Format',formatSpec);
    %T
    %head(T,5) 
    
end