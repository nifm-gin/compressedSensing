

function RecData = KSpaceDataComb( AquSpec, FilData, Mask, Coff, Mode )

% mode: 1. noiseless, 2. noise

DataSize = size( AquSpec );
FilSpec  = FFT2_3D( FilData, ones( DataSize ), 0 );

RecSpec = zeros( DataSize );
for index = 1 : DataSize(3)
    FilMax = max( max( FilSpec( :, :, index ) ) );
    AquMax = max( max( AquSpec( :, :, index ) ) );
    
    D1 = ( 1 - Mask( :, :, index ) ) .* FilSpec( :, :, index ) / FilMax;
    D2 = Mask( :, :, index ) .* AquSpec( :, :, index ) / AquMax;
    
    D = ( D1 + D2 ) * AquMax;
    RecSpec( :, :, index ) = D;
end
RecData = IFFT2_3D( RecSpec, ones( DataSize ), 0 );

% % 
% % 
% % AquMax = max
% % 
% % D1 = ( 1 - Mask ) .* FilSpec;
% % 
% % if Mode == 1
% %     D2 = Mask .* AquSpec;
% % else
% %     D2 = Mask .* ( FilSpec + Coff * AquSpec ) / ( 1 + Coff );
% % end
% % 
% % RecSpec = D1 + D2;
% % RecData = IFFT2_3D( RecSpec, ones( DataSize ), 0 );
