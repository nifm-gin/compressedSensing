
function [ PrjCoff, Alpha, D, KxxOut, Kxc, Kxy, Kyc ] = KPCA( PrjData, TraData, Para )

PrjNum = size( PrjData, 2 );
TraNum = size( TraData, 2 );

%% Trainsing kernel matrix
if ~isfield( Para, 'Kxx' )
    KxxOut = CalcKernelMatrixXY( TraData, TraData, Para );
    Kxc = CentralKernelXX( KxxOut );
else
    KxxOut = Para.Kxx;
    Kxc = Para.Kxc;
end

%% eigen value analysis and normailization
if ~isfield( Para, 'Alpha' )
    [ V, D ] = eig( Kxc );
    D = diag( D ).';
    
    Vn = repmat( D, TraNum, 1 );
    Vn = sqrt( abs( Vn ) );
    Mn = double( Vn ~= 0 );
    V  = ( V ./ Vn ) .* Mn;
    
    [ D, IX ] = sort( D, 'descend' );
    Alpha = V( :, IX );
% %     Alpha = V;
else
    Alpha = Para.Alpha;
    D     = Para.D;
end

%% projection
if ( PrjNum == TraNum ) && isequal( PrjData, TraData )
    Kxy = KxxOut;
    Kyc = Kxc;
else
    Kxy = CalcKernelMatrixXY( TraData, PrjData, Para );   % X-Y kernel
    Kyc = CentralKernelXY( Kxy, KxxOut );                    % centralized
end

PrjCoff = Kyc * Alpha;
PrjCoff = PrjCoff.';
