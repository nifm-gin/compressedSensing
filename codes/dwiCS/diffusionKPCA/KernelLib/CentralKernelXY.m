

function Kc = CentralKernelXY( Kxy, Kxx )

L  = size( Kxy, 1 );
N  = size( Kxx, 1 );

D1 = sum( Kxx, 1 ) / N;
D1 = repmat( D1, L, 1 );

D2 = sum( Kxy, 2 ) / N;
D2 = repmat( D2, 1, N );

D3 = sum( D1, 2 ) / N;
D3 = repmat( D3, 1, N );

Kc = Kxy - D1 - D2 + D3;
