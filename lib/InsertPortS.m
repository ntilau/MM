function Stot = InsertPortS( Stot, S, nbrModes, index_i ,index_j)
% Insert one block scattering submatrix into a global GSM.
% ------------------------------------------------------------------------
% [IN]
%   Stot     - global scattering matrix arranged in block form
%   S        - scattering submatrix to insert
%   nbrModes - size of each block scattering submatrix
%   index_i  - output-port block index
%   index_j  - input-port block index
%
% [OUT]
%   Stot     - updated global scattering matrix
%
Stot( (1+(index_i-1)*nbrModes(1)) : (index_i*nbrModes(1)), ...
    (1+(index_j-1)*nbrModes(2)) : (index_j*nbrModes(2)) ) = S(:,:);
