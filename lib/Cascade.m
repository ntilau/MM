function [ Stot ] = Cascade( WaveGuideStructure, S )
% Cascade all step GSMs and uniform-section delay matrices for a complete
% WaveGuideStructure.
% ------------------------------------------------------------------------
% [IN]
%   WaveGuideStructure - cell array of WaveGuideSegment structs with delay
%                        matrices already computed for each section
%   S                  - cell array of step GSM structs
%
% [OUT]
%   Stot               - overall GSM struct for the full cascade
%
S0.S11 = zeros(WaveGuideStructure{1}.Nh+WaveGuideStructure{1}.Ne);
S0.S12 = eye(WaveGuideStructure{1}.Nh+WaveGuideStructure{1}.Ne);
S0.S21 = eye(WaveGuideStructure{1}.Nh+WaveGuideStructure{1}.Ne);
S0.S22 = zeros(WaveGuideStructure{1}.Nh+WaveGuideStructure{1}.Ne);

% Taking this into account
Stot = SingleCascade(S0,WaveGuideStructure{1}.D,S{1});

% Taking all internal steps into account
Np = length(WaveGuideStructure);
if(Np>1)
    for p=2:(Np-1)
        Stot = SingleCascade(Stot,WaveGuideStructure{p}.D,S{p});
    end

    % Last ficticious step
    S0.S11 = zeros(WaveGuideStructure{Np}.Nh+WaveGuideStructure{Np}.Ne);
    S0.S12 = eye(WaveGuideStructure{Np}.Nh+WaveGuideStructure{Np}.Ne);
    S0.S21 = eye(WaveGuideStructure{Np}.Nh+WaveGuideStructure{Np}.Ne);
    S0.S22 = zeros(WaveGuideStructure{Np}.Nh+WaveGuideStructure{Np}.Ne);

    Stot = SingleCascade(Stot,WaveGuideStructure{Np}.D,S0);

end
