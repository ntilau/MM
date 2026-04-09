function [WaveGuideSegment, OpenPortSegment,  WaveGuideStructure] = ...
    TwoPortDeviceGetPortSegment(WaveGuideStructure , PortRelated)
% Return the junction-facing segment and the open-port-facing segment for
% one end of a TwoPortDevice.
% ------------------------------------------------------------------------
% [IN]
%   WaveGuideStructure - cell array of WaveGuideSegment structs
%   PortRelated        - 1 for the first port, 2 for the second port
%
% [OUT]
%   WaveGuideSegment   - segment used in the global junction assembly
%   OpenPortSegment    - segment representing the opposite external port
%   WaveGuideStructure - returned unchanged for convenience
%
NSegment = length(WaveGuideStructure);
if(PortRelated == 1)
    WaveGuideSegment = WaveGuideStructure{1};
    OpenPortSegment = WaveGuideStructure{NSegment};
elseif(PortRelated == 2)
    OpenPortSegment = WaveGuideStructure{1};
    WaveGuideSegment = WaveGuideStructure{NSegment};
else
    Error.fatal = sprintf('Could not retrieve Port Segment',i+1);
    return;
end

