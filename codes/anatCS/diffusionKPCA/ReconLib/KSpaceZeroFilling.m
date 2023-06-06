

function [ Data3D, Mask3D ] = KSpaceZeroFilling( ImgSize, PhaseIdx, RawData )

Data3D = zeros( ImgSize );
Mask3D = zeros( ImgSize );

for index = 1 : ImgSize(3)
    Data3D( PhaseIdx( :, index ), :, index ) = RawData( :, :, index );
    Mask3D( PhaseIdx( :, index ), :, index ) = 1;
end