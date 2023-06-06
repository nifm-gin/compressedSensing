
%
% M-Function  explanation:
%   Centralize the training kernel matrix
%
% Input:
%    Kxx ---- non-centered training kernel matrix
%             N by N matrix

function Kc = CentralKernelXX( Kxx )

N  = size( Kxx, 1 );

D1 = sum( Kxx, 1 ) / N;
D1 = repmat( D1, N, 1 );

D2 = sum( Kxx, 2 ) / N;
D2 = repmat( D2, 1, N );

D3 = sum( D1, 2 ) / N;
D3 = repmat( D3, 1, N );

Kc = Kxx - D1 - D2 + D3;
