
function Y = FFT2_4D( X, Mask, IfShift )

DataSize = size( X );
%DataFFT  = fft2( X )/sqrt( DataSize(1) * DataSize(2) );
DataFFT  = fft2( X );
if IfShift == 1
    DataFFT = FFT2Shift_4D( DataFFT );
end

Y = DataFFT .* Mask;