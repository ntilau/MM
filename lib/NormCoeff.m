function [ WaveGuideSegment ] = NormCoeff( WaveGuideSegment )
% Compute mode-normalization coefficients for all TE and TM modes in one
% waveguide segment.
% ------------------------------------------------------------------------
% [IN]
%   WaveGuideSegment - struct containing geometry and modal eigenvalues
%
% [OUT]
%   WaveGuideSegment - input struct augmented with Ah and Ae
%
Ah = [];
for i=1:WaveGuideSegment.Nh
    Ah(i) = ...
          OneModeNormCoeff(WaveGuideSegment.a, WaveGuideSegment.b, ...
          WaveGuideSegment.mh(i),WaveGuideSegment.nh(i), ...
          WaveGuideSegment.kh.x(i),WaveGuideSegment.kh.y(i));
end
Ae = [];
for i=1:WaveGuideSegment.Ne
    Ae(i) = ...
          OneModeNormCoeff(WaveGuideSegment.a, WaveGuideSegment.b, ...
          WaveGuideSegment.me(i),WaveGuideSegment.ne(i), ...
          WaveGuideSegment.ke.x(i),WaveGuideSegment.ke.y(i));
end

WaveGuideSegment.Ah = Ah;
WaveGuideSegment.Ae = Ae;

