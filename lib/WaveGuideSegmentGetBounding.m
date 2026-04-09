function [xmin,xmax,ymin,ymax] =...
    WaveGuideSegmentGetBounding(TwoPortSegment)
% Return the axis-aligned bounding box of a waveguide segment cross
% section.
% ------------------------------------------------------------------------
% [IN]
%   TwoPortSegment - WaveGuideSegment struct with a, b, xo, yo
%
% [OUT]
%   xmin, xmax, ymin, ymax - bounding box coordinates [m]
%
[a,b,xo,yo] = WaveGuideSegmentGetCrossSection(TwoPortSegment);

xmin = xo - a/2;
xmax = xo + a/2;
ymin = yo - b/2;
ymax = yo + b/2;
