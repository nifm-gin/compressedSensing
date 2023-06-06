
function Y = FFT2Shift_4D( X )

ImgSize = size( X );
Y = zeros( ImgSize );

for ii = 1:ImgSize(4)
for index = 1 : ImgSize(3)
    Y( :, :, index,ii ) = fftshift( X( :, :, index,ii ) );
end
end
