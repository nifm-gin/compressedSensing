function res = sos_mod(x ,dim, pnorm)
% res = sos(x [,dim, pnorm])
%
% function computes the square root of sum of squares along dimension dim.
% If dim is not specified, it computes it along the last dimension.
%
% (c) Michael Lustig 2009

if nargin < 2
    dim = size(size(x),2);
end

if nargin < 3
    pnorm = 2;
end

% Original ---------------------------------------
%res = squeeze(sum(abs(x.^pnorm),dim)).^(1/pnorm);
% ------------------------------------------------

% Added on 2022.01.07 ----------------------------
res = squeeze(sum((x.^pnorm),dim)).^(1/pnorm);
% ------------------------------------------------