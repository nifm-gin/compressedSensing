

function [ Err, MSE, NMSE ] = ErrAnalysis( RefData, CmpData, IfShow, ShowLim )

DataSize = size( RefData );

Err  = abs( CmpData ) - abs( RefData );
MSE  = sum( sum( abs( Err ).^2, 1 ), 2 );
MSE  = MSE(:);
NMSE = MSE / ( DataSize(1) * DataSize(2) );

if IfShow == 1
    for index = 1 : DataSize(3)
        figure;
        imagesc( abs( Err( :, :, index ) ), ShowLim ); colorbar;
        title( [ 'Frame:' num2str( index ) ] );
    end
    
    figure;
    plot( MSE, 'r-o', 'LineWidth', 2 );
    title( 'MSE vs Frame' );
    
    figure;
    plot( NMSE, 'r-o', 'LineWidth', 2 );
    title( 'NMSE vs Frame' );
end
