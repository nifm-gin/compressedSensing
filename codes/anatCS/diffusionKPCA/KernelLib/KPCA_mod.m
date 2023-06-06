
function [ PrjCoff, Alpha, D, KxxOut, Kxc, Kxy, Kyc ] = KPCA_mod( PrjData, TraData, Para )
fprintf("\n\n ===============================================\n KPCA modified\n ===============================================\n")

PrjNum = size( PrjData, 2 );
TraNum = size( TraData, 2 );

%% Trainsing kernel matrix
tic
if ~isfield( Para, 'Kxx' )
    %fprintf("-----------------------------------\n Step 1.1: KxxOut = CalcKernelMatrixXY \n-----------------------------------\n")
    KxxOut = CalcKernelMatrixXY( TraData, TraData, Para );
    Kxc = CentralKernelXX( KxxOut );
else
    %fprintf("-----------------------------------\n Step 1.2: KxxOut = Para.Kxx \n-----------------------------------\n")
    KxxOut = Para.Kxx;
    Kxc = Para.Kxc;
end
toc

%% eigen value analysis and normailization
if ~isfield( Para, 'Alpha' )
    %fprintf("-----------------------------------\n Step 2.1: [ V, D ] = eig( Kxc ) \n-----------------------------------\n")
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
    %fprintf("-----------------------------------\n Step 2.2: Alpha = Para.Alpha \n-----------------------------------\n")
    Alpha = Para.Alpha;
    D     = Para.D;
end

%% projection
if ( PrjNum == TraNum ) && isequal( PrjData, TraData )
    %fprintf("-----------------------------------\n Step 3.1: Kxy = KxxOut \n-----------------------------------\n")
    Kxy = KxxOut;
    Kyc = Kxc;
else
    %fprintf("-----------------------------------\n Step 3.2: Kxy = CalcKernelMatrixXY \n-----------------------------------\n")
    Kxy = CalcKernelMatrixXY( TraData, PrjData, Para );   % X-Y kernel
    Kyc = CentralKernelXY( Kxy, KxxOut );                    % centralized
end

%fprintf("-----------------------------------\n Step 4: PrjCoff = Kyc * Alpha \n-----------------------------------\n")
PrjCoff = Kyc * Alpha;
PrjCoff = PrjCoff.';
