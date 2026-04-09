function WaveGuideStructure = ...
    TwoPortDeviceInsertPortSegment( WaveGuideStructureToChange, WaveGuideSegment, Position)
% Insert a waveguide segment at one end of a TwoPortDevice structure.
% ------------------------------------------------------------------------
% [IN]
%   WaveGuideStructureToChange - device structure to extend
%   WaveGuideSegment           - segment to insert
%   Position                   - 1 prepend, 2 append
%
% [OUT]
%   WaveGuideStructure         - updated device structure
%
NSegment = length(WaveGuideStructureToChange);
if(Position == 1)
        WaveGuideStructure.D{2:NSegment+1} = WaveGuideStructureToChange.D{1:NSegment};
        WaveGuideStructure.D{1} = WaveGuideSegment;
elseif(Position == 2)
    WaveGuideStructure.D{1:NSegment} = WaveGuideStructureToChange.D{1:NSegment};
    WaveGuideStructure.D{NSegment+1} = WaveGuideSegment;
else
    Error.fatal = sprintf('Could not insert Port Segment',i+1);
    return;
end
