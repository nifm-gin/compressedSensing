
function Y = IFFT2_4D( X, Mask, IfShift )

DataSize = size( X );
%DataFFT  = X .* Mask * sqrt( DataSize(1) * DataSize(2) );
DataFFT  = X .* Mask;
if IfShift == 1
    DataFFT = IFFT2Shift_4D( DataFFT );
end

Y = ifft2( DataFFT );