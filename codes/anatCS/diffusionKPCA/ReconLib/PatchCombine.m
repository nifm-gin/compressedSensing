%
% Function:
%   [ PatchData PatchNum ] = PatchExtract( Data, PatchSize, Step )
%
% Explanation:
%   Divide original data into patches patch from the original data
%
% Input:
%   Data  ---------- Original data, maximum dimension is 3
%   PatchSize ------ Size of a single patch
%   Step   --------- Step between ajacent patches
%
% Ouput:
%   PatchData ------ The required patch data
%   PatchNum  ------ Total number of patches
%
% Note:
%
%
% Version History:
%   2012-06-27  ---  v1.0 --- by WYH
%     Creat the file
%
% Copyright:
%   Yanhua Wang, Ph.D, EE, UB-SUNY
%

function Data = PatchCombine( PatchData, DataSize, PatchSize, Step )

%% Parameter check
DataDims = length( DataSize );

if DataDims > 3
    error( 'Data dimension exceeds 3!' );
end

if DataDims == 2
    DataDims = 3;
    DataSize(3) = 1;
end

if length( PatchSize ) ~= DataDims
    error( 'Patch size and data dimension are not matched!' );
end

if length( Step ) ~= DataDims
    error( 'Patch step and data dimension are not matched!' );
end

for index = 1 : DataDims
    if PatchSize( index ) > DataSize( index )
        fprintf( 'Patch exceeds data dimension in %d\n', index );
    end
end

PatNum    = ceil( ( DataSize - PatchSize ) ./ ( Step + eps ) ) + 1;
TotPatNum = prod( PatNum );
InPatNum  = size( PatchData );

if InPatNum(2) > TotPatNum
    erros( 'Patch number exceeds the maximum available number!' );
end

%% Combination
ExtSize = ( PatNum - 1 ) .* Step + PatchSize;

if DataDims == 2
    DataSize(3)  = 1;
    ExtSize(3)   = 1;
    PatchSize(3) = 1;
    PatNum(3)    = 1;
    Step(3)      = 0;
end

ExtData = zeros( ExtSize );
CntData = zeros( ExtSize );

Index_D1 = 1 : PatchSize(1);
PatchIdx = 1;

for D1_Index = 1 : PatNum(1)
    
    Index_D2 = 1 : PatchSize(2);
    
    for D2_Index = 1 : PatNum(2)
        
        Index_D3 = 1 : PatchSize(3);
        
        for D3_Index = 1 : PatNum(3)
            
            T1 = PatchData( :, PatchIdx );
            T2 = reshape( T1, PatchSize(1), PatchSize(2), PatchSize(3) );
            T3 = ExtData( Index_D1, Index_D2, Index_D3 );
            T4 = T2 + T3;
            ExtData( Index_D1, Index_D2, Index_D3 ) = T4;
            
            T3 = CntData( Index_D1, Index_D2, Index_D3 );
            T4 = T3 + 1;
            CntData( Index_D1, Index_D2, Index_D3 ) = T4;
            
            Index_D3 = Index_D3 + Step(3);
            PatchIdx = PatchIdx + 1;
            
        end
        
        Index_D2 = Index_D2 + Step(2);
        
    end
    
    Index_D1 = Index_D1 + Step(1);
    
end


%% Finalization
% % ExtData = ExtData ./ CntData;
Data = ExtData( 1:DataSize(1), 1:DataSize(2), 1:DataSize(3) );

