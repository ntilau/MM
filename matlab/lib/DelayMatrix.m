function [ WaveGuideSegment ] = DelayMatrix( WaveGuideSegment, k0 )
% Build the diagonal propagation delay matrix for all TE and TM modes in
% one uniform waveguide segment.
% ------------------------------------------------------------------------
% [IN]
%   WaveGuideSegment - WaveGuideSegment struct with modal wavenumbers kh, ke
%   k0               - free-space wavenumber [rad/m]
%
% [OUT]
%   WaveGuideSegment - input struct augmented with delay matrix D for the
%                      current frequency
%
D = [[diag(exp(-1i*WaveGuideSegment.kh.mn*WaveGuideSegment.l)),...
          zeros(WaveGuideSegment.Nh,WaveGuideSegment.Ne)];...
     [zeros(WaveGuideSegment.Ne,WaveGuideSegment.Nh),...
          diag(exp(-1i*WaveGuideSegment.ke.mn*WaveGuideSegment.l))]];

WaveGuideSegment.D = D;
