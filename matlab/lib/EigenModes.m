function [ WaveGuideSegment ] = EigenModes( WaveGuideSegment )
% Compute the transverse eigenvalues for all TE and TM modes of a
% rectangular waveguide segment.
% ------------------------------------------------------------------------
% [IN]
%   WaveGuideSegment - struct containing a, b and the mode index lists
%                      mh, nh, me, ne
%
% [OUT]
%   WaveGuideSegment - input struct augmented with kh.x, kh.y, ke.x, ke.y
%                      for all retained modes
%
for i=1:WaveGuideSegment.Nh
    [kh.x(i), kh.y(i)] ...
        = OneModeEigens(WaveGuideSegment.a, WaveGuideSegment.b, ...
          WaveGuideSegment.mh(i),WaveGuideSegment.nh(i));
end
ke.x = [];
ke.y = [];
for i=1:WaveGuideSegment.Ne
    [ke.x(i), ke.y(i)] ...
        = OneModeEigens(WaveGuideSegment.a, WaveGuideSegment.b, ...
          WaveGuideSegment.me(i),WaveGuideSegment.ne(i));
end

WaveGuideSegment.kh = kh;
WaveGuideSegment.ke = ke;

