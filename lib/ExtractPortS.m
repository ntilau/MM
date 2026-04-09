function S = ExtractPortS( Stot, SeleDim, index_i ,index_j)
% Extract one block scattering submatrix from a global GSM.
% ------------------------------------------------------------------------
% [IN]
%   Stot    - global scattering matrix arranged in block form
%   SeleDim - size of each block scattering submatrix
%   index_i - output-port block index
%   index_j - input-port block index
%
% [OUT]
%   S       - scattering submatrix from port j to port i
%
S = Stot( (1+(index_i-1)*SeleDim(1)) : (index_i*SeleDim(1)), ...
    (1+(index_j-1)*SeleDim(2)) : (index_j*SeleDim(2)) );
