
%
% M-File  explanation:
%   A part of CS-DL-dMRI algorithm
%   Analyze the performance
%
% Version History:
%   2012-06-28  --  v1.0 -- by WYH
%     Creat the file
%
% Copyright:
%   Yanhua Wang, PhD
%   Postdoc Research Associate
%   Department of Electrical Engineering
%   University at Buffalo, the State University of New York
%   Contact:  wyhlucky@gmail.com
%

%%
RefDataInfo = load( TrDataPath );

RefData = RefDataInfo.ImageData;
CmpData = FinalResult;

ErrAnalysis( RefData, CmpData, 1, [ 0 0.1 ] );
