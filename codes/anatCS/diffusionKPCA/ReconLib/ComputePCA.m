
function [ PC, D ] = ComputePCA( TrainData )

DataSize = size( TrainData );

PCASort = permute( TrainData, [ 3 1 2 ] );
PCASort = reshape( PCASort, [ DataSize(3) DataSize(1) * DataSize(2) ] );
PCAMean = mean( PCASort, 2 );
PCASort = PCASort - repmat( PCAMean, [ 1 DataSize(1) * DataSize(2) ] );
CovData = PCASort * PCASort';
[ PC D ] = eig( CovData );
