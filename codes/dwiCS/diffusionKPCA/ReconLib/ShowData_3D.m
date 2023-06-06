

function ShowData_3D( Data, ShowFunc, ShowMode )

DataSize = size( Data );

if ShowMode == 1
    
    AxeNum   = ceil( sqrt( DataSize(3) ) );
    
    figure;
    for index = 1 : DataSize(3)
        
        subplot( AxeNum, AxeNum, index );
        if strcmp( ShowFunc, 'imagesc' )
            imagesc( abs( Data( :, :, index ) ) ); colorbar
        elseif strcmp( ShowFunc, 'imshow' )
            imshow( abs( Data( :, :, index ) ), [] );
        else
        end
        
        title( num2str( index ) );
    end
    
else
    for index = 1 : DataSize(3)
        
        figure;
        if strcmp( ShowFunc, 'imagesc' )
            imagesc( abs( Data( :, :, index ) ) ); colorbar
        elseif strcmp( ShowFunc, 'imshow' )
            imshow( abs( Data( :, :, index ) ), [] );
        else
        end
        
        title( num2str( index ) );
    end
end