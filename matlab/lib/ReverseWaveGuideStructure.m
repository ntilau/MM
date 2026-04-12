function WGSout = ReverseWaveGuideStructure(WGSin)
% Reverse the segment order of a WaveGuideStructure cell array.
% ------------------------------------------------------------------------
% [IN]
%   WGSin  - cell array of WaveGuideSegment structs ordered from port 1
%            to port 2
%
% [OUT]
%   WGSout - same segments in reverse order
%
WGSout = {};
NbrSegments = length(WGSin);

Indices = 1:NbrSegments;
newIndices = sort(Indices, 'descend');

for i=1:NbrSegments
    WGSout{i} = WGSin{newIndices(i)};
end
