function [ a, b, xo, yo ] = WaveGuideSegmentGetCrossSection( TwoPortSegment )
% Extract the rectangular cross-section parameters of a waveguide segment.
% ------------------------------------------------------------------------
% [IN]
%   TwoPortSegment - WaveGuideSegment struct
%
% [OUT]
%   a, b   - waveguide width and height [m]
%   xo, yo - centre offsets [m]
%
a = TwoPortSegment.a;
b = TwoPortSegment.b;
xo = TwoPortSegment.xo;
yo = TwoPortSegment.yo;
