

function Z = CalcPreImage( InitData, TraData, Alpha, PrjCoff, Para )

MaxItrNum = Para.ItrNum;
[ SigDim, SigNum ] = size( TraData );

Z = InitData;

if strcmp( Para.Mode, 'Poly' )
    
elseif strcmp( Para.Mode, 'Gauss' )

	Gamma = Alpha * PrjCoff;
	GaCen = ( 1 - sum( Gamma, 1 ) ) / SigNum;
	GaCen = repmat( GaCen, SigNum, 1 );
	Gamma = Gamma + GaCen;
	
	for ItrIdx = 1 : MaxItrNum
		
		Kxz = CalcKernelMatrixXY( TraData, Z, Para );
		
		VCof = ( Kxz.' ) .* Gamma;
		
		NCof = Kxz * Gamma;
		NCof = diag( NCof ).';
		NCof = repmat( NCof, SigDim, 1 );
		
		Z = TraData * VCof ./ NCof;
		
	end
    
else
    
end

H = eye(SigNum) - 1/SigNum * ones(SigNum);
M = Alpha * Alpha';