function [ A ] = OneModeNormCoeff( a, b, m, n, kx, ky )
% Compute the normalization coefficient for one rectangular-waveguide
% mode.
% ------------------------------------------------------------------------
% [IN]
%   a, b   - waveguide dimensions [m]
%   m, n   - mode indices
%   kx, ky - transverse eigenvalues for the selected mode
%
% [OUT]
%   A      - modal normalization coefficient
%
if (m==0)
    dm = 1;
else
    dm = 0;
end

if (n==0)
    dn = 1;
else
    dn = 0;
end

A = 2 / ( sqrt(a*b) * ...
    sqrt(kx^2*(1+dn)+ky^2*(1+dm)) );
