function [ z0 ] = ShowSegment( WaveGuideSegment, z0 )
% Draw one rectangular waveguide segment as a set of 3-D patches.
% ------------------------------------------------------------------------
% [IN]
%   WaveGuideSegment - struct containing a, b, l, xo, yo
%   z0               - starting z-position [m]
%
% [OUT]
%   z0               - ending z-position after drawing the segment
%
a = WaveGuideSegment.a;
b = WaveGuideSegment.b;
l = WaveGuideSegment.l;
xo = WaveGuideSegment.xo;
yo = WaveGuideSegment.yo;

patch([xo-a/2, xo+a/2, xo+a/2, xo-a/2],...
    [yo-b/2, yo-b/2, yo-b/2, yo-b/2],...
    [z0, z0, z0+l, z0+l],'b');

patch([xo-a/2, xo+a/2, xo+a/2, xo-a/2],...
    [yo+b/2, yo+b/2, yo+b/2, yo+b/2],...
    [z0, z0, z0+l, z0+l],'b');

patch([xo-a/2, xo-a/2, xo-a/2, xo-a/2],...
    [yo-b/2, yo+b/2, yo+b/2, yo-b/2],...
    [z0, z0, z0+l, z0+l],'b');

patch([xo+a/2, xo+a/2, xo+a/2, xo+a/2],...
    [yo-b/2, yo+b/2, yo+b/2, yo-b/2],...
    [z0, z0, z0+l, z0+l],'b');

z0 = z0+l;
