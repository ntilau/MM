function [Stot, S] = Nto1Junction( WaveGuideStructure_1, WaveGuideStructure_2, k0)
% Compute the GSM of an N-to-1 rectangular-waveguide junction using the
% modal-matching formulation.
% ------------------------------------------------------------------------
% [IN]
%   WaveGuideStructure_1 - cell array of input-side WaveGuideSegment structs
%   WaveGuideStructure_2 - output-side WaveGuideSegment struct
%   k0                   - free-space wavenumber [rad/m]
%
% [OUT]
%   Stot                 - full block scattering matrix
%   S                    - GSM struct split into S11, S12, S21, S22
%
NSegments = length(WaveGuideStructure_1); % Nbr of input segments
step{2} = WaveGuideStructure_2; % Signal output gate

H=[];
E=[];
LengthofUnitMatrix = 0;
lambda = 2*pi/k0;
Z_lambda = 1i*2*pi*120*pi/lambda;
Y_lambda = 1i*2*pi/(120*pi*lambda);

for i=1:NSegments
    step{1} = WaveGuideStructure_1{i}; % Signal input gates

    [Mhh{i}, Mhe{i}, Meh{i}, Mee{i}] = MxxMatrices(step);
    % No need to save N times Dgamma and U for the output port
    [Dgamma_h_1{i}, Dgamma_h_2, Dgamma_e_1{i}, Dgamma_e_2] = DgammaMatrices(step);
    [U_h_1{i}, U_h_2, U_e_1{i}, U_e_2] = UMatrices(step);

    DHcell{i} = [[Dgamma_h_1{i}, zeros(step{1}.Nh, step{1}.Ne)];...
        [zeros(step{1}.Ne, step{1}.Nh), U_e_1{i}]];
    MHcell{i} = [[Mhh{i} * Dgamma_h_2, Y_lambda * Mhe{i}];...
        [Meh{i} * Dgamma_h_2 / Y_lambda, Mee{i}]];
    DEcell{i} = [[U_h_2, zeros(step{2}.Nh, step{2}.Ne)]; ...
        [zeros(step{2}.Ne, step{2}.Nh),Dgamma_e_2]];
    MEcell{i} = [[transpose(Mhh{i}) , transpose(Meh{i}) * Dgamma_e_1{i} / Z_lambda];...
        [Z_lambda*transpose(Mhe{i}), transpose(Mee{i}) * Dgamma_e_1{i}]];

    Hcell{i} = inv(DHcell{i})*MHcell{i};
    Ecell{i} = inv(DEcell{i})*MEcell{i};

    H = vertcat(H,Hcell{i}); % concatenates vertically -> size(H) = [N*modes x modes]
    E = horzcat(E,Ecell{i}); % size(E) = [modes x N*modes]
    LengthofUnitMatrix = LengthofUnitMatrix + WaveGuideStructure_1{i}.Nh + WaveGuideStructure_1{i}.Ne;
end

HE = H*E; % size(HE) = N*[ modes x modes]
UnitMatrix = eye(LengthofUnitMatrix); % Unit for the inputs
UnitOutMatrix = eye(WaveGuideStructure_2.Nh + WaveGuideStructure_2.Ne); % Unit for the output
UnitMatrix_plus_HE_inverse = inv( UnitMatrix + HE ); % 1/(1+HE )

S.S11 = UnitMatrix_plus_HE_inverse *( UnitMatrix - HE);
S.S12 = 2*UnitMatrix_plus_HE_inverse*H;
S.S21 = 2*E*UnitMatrix_plus_HE_inverse;
S.S22 = E*S.S12 - UnitOutMatrix;

SeleDim = [WaveGuideStructure_1{1}.Nh+WaveGuideStructure_1{1}.Ne,...
    WaveGuideStructure_2.Nh+WaveGuideStructure_2.Ne];

% Controls condensability
if( SeleDim(1) ~= SeleDim(2))
    Warning('Wont be able to condensate');
    return;
end

Stot = [[S.S11, S.S12] ; [S.S21, S.S22]];
