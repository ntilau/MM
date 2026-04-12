function [ Sf, WaveGuideStructure ] = MultiStep( WaveGuideStructure, k0, Symmetry)
% Compute the GSM of a multi-section WaveGuideStructure at one frequency.
% ------------------------------------------------------------------------
% [IN]
%   WaveGuideStructure - cell array of WaveGuideSegment structs
%   k0                 - free-space wavenumber [rad/m]
%   Symmetry           - symmetry struct used for mode ordering
%
% [OUT]
%   Sf                 - GSM struct for the full structure
%   WaveGuideStructure - input structure updated with frequency data
%
nbrSegments = length(WaveGuideStructure);

for p=1:nbrSegments
%     [WaveGuideStructure{p},Error] = OrderModes(WaveGuideStructure{p},Symmetry);
%     WaveGuideStructure{p} = EigenModes(WaveGuideStructure{p});
%     WaveGuideStructure{p} = NormCoeff(WaveGuideStructure{p});
    WaveGuideStructure{p} = WaveNumbers(WaveGuideStructure{p},k0);
    WaveGuideStructure{p} = DelayMatrix(WaveGuideStructure{p},k0);
end

if(nbrSegments>1)
    for p=1:nbrSegments-1;
        SStep{p} = SingleStep (WaveGuideStructure, p, k0);
    end
else
    SStep{1}.S11 = zeros(WaveGuideStructure{1}.Nh+WaveGuideStructure{1}.Ne);
    SStep{1}.S12 = eye(WaveGuideStructure{1}.Nh+WaveGuideStructure{1}.Ne);
    SStep{1}.S21 = eye(WaveGuideStructure{1}.Nh+WaveGuideStructure{1}.Ne);
    SStep{1}.S22 = zeros(WaveGuideStructure{1}.Nh+WaveGuideStructure{1}.Ne);
end

Stot = Cascade(WaveGuideStructure, SStep);

% Stot = Renormalize(WaveGuideStructure,Stot);

Sf = [ [Stot.S11, Stot.S12] ; [Stot.S21, Stot.S22] ];
