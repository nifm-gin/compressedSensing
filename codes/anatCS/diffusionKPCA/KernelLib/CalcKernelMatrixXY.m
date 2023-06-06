

function Kxy = CalcKernelMatrixXY( X, Y, Para )

[ ~, nX ] = size( X );
[ ~, nY ] = size( Y );

if strcmp( Para.Mode, 'Poly' )
    
    c = Para.c;
    d = Para.d;
    Kxy = ( (Y') * X + c ).^d;
end
