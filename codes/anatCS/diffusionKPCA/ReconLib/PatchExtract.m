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

function [ PatchData TotPatNum ] = PatchExtract( Data, PatchSize, Step )

%% Parameter check
DataSize = size( Data );
DataDims = ndims( Data );

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

%% Preparation
PatNum    = ceil( ( DataSize - PatchSize ) ./ ( Step + eps ) ) + 1;
TotPatNum = prod( PatNum );

ExtSize = ( PatNum - 1 ) .* Step + PatchSize;

if DataDims == 2
    DataSize(3)  = 1;
    ExtSize(3)   = 1;
    PatchSize(3) = 1;
    PatNum(3)    = 1;
    Step(3)      = 0;
end

ExtData = zeros( ExtSize );

ExtData( 1:DataSize(1), 1:DataSize(2), 1:DataSize(3) ) = Data;

for dim = 1 : 3
    if ExtSize(dim) - DataSize(dim) > 0
        for index = 1 : ExtSize(dim) - DataSize(dim)
            ExtData( DataSize(dim) + index, :, : ) = ExtData( DataSize(dim), :, : );
        end
    end
end

%% Extraction
PatchData = zeros( prod( PatchSize ), TotPatNum );

PatchIdx = 1;
Index_D1 = 1 : PatchSize(1);
for D1_Index = 1 : PatNum(1)
    
    Index_D2 = 1 : PatchSize(2);
    for D2_Index = 1 : PatNum(2)
        
        Index_D3 = 1 : PatchSize(3);       
        for D3_Index = 1 : PatNum(3)
            
            Temp = ExtData( Index_D1, Index_D2, Index_D3 );
            PatchData( :, PatchIdx ) = Temp(:);
            Index_D3 = Index_D3 + Step(3);
            PatchIdx = PatchIdx + 1;
            
        end
        
        Index_D2 = Index_D2 + Step(2);
        
    end
    
    Index_D1 = Index_D1 + Step(1);
    
end
