function [ kx, ky ] = OneModeEigens( a, b, m, n )
% Compute the transverse eigenvalues for mode (m,n) in a rectangular
% waveguide of size a by b.
% ------------------------------------------------------------------------
% [IN]
%   a, b - waveguide dimensions [m]
%   m, n - mode indices
%
% [OUT]
%   kx, ky - transverse eigenvalues along x and y [rad/m]
%
kx = m*pi/a;

ky = n*pi/b;

