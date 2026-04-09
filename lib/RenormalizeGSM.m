function [Stot] = RenormalizeGSM(S, SeleDim, WaveGuideStructure_1, WaveGuideStructure_2)
% Renormalize every block of a global scattering matrix for an N+1-port
% junction.
% ------------------------------------------------------------------------
% [IN]
%   S                    - global scattering matrix
%   SeleDim              - size of each block scattering submatrix
%   WaveGuideStructure_1 - cell array of input-side waveguide structures
%   WaveGuideStructure_2 - output-side waveguide structure
%
% [OUT]
%   Stot                 - renormalized global scattering matrix
%
SDimTmp = size(S)./SeleDim;
if(SDimTmp(1) == SDimTmp(2))
    SDim = SDimTmp(1);
else
    warning('S is not a square matrix!');
end
% Renormalizes i to j
for i=1:length(WaveGuideStructure_1)
    for j=(i+1):length(WaveGuideStructure_1)
        WGS{i} = WaveGuideStructure_1{i};
        WGS{j} = WaveGuideStructure_1{j};
        Stmp.S11 = ExtractPortS(S, SeleDim, i, i);
        Stmp.S12 = ExtractPortS(S, SeleDim, i, j);
        Stmp.S21 = ExtractPortS(S, SeleDim, j, i);
        Stmp.S22 = ExtractPortS(S, SeleDim, j, j);
        StmpRN = Renormalize(WGS, Stmp);
        S = InsertPortS(S, StmpRN.S11, SeleDim, i, i);
        S = InsertPortS(S, StmpRN.S12, SeleDim, i, j);
        S = InsertPortS(S, StmpRN.S21, SeleDim, j, i);
        S = InsertPortS(S, StmpRN.S22, SeleDim, j, j);
    end
end

%Renormalizes i to N+1
% if(iscell(WaveGuideStructure_2))
%     WGS{2} = WaveGuideStructure_2{1};
% else
%     WGS{2} = WaveGuideStructure_2;
% end
for o=1:length(WaveGuideStructure_2)
    if(iscell(WaveGuideStructure_2))
        WGS{2} = WaveGuideStructure_2{o};
    else
        WGS{2} = WaveGuideStructure_2;
    end
    for i=1:length(WaveGuideStructure_1)
        WGS{1} = WaveGuideStructure_1{i};
        Stmp.S11 = ExtractPortS(S, SeleDim, i, i);
        Stmp.S12 = ExtractPortS(S, SeleDim, i, SDim);
        Stmp.S21 = ExtractPortS(S, SeleDim, SDim, i);
        Stmp.S22 = ExtractPortS(S, SeleDim, SDim, SDim);
        StmpRN = Renormalize(WGS, Stmp);
        S = InsertPortS(S, StmpRN.S11, SeleDim, i, i);
        S = InsertPortS(S, StmpRN.S12, SeleDim, i, SDim);
        S = InsertPortS(S, StmpRN.S21, SeleDim, SDim, i);
        S = InsertPortS(S, StmpRN.S22, SeleDim, SDim, SDim);
    end
end
% if(iscell(WaveGuideStructure_2))
%     WGS{1} = WaveGuideStructure_2{1};
% else
%     WGS{1} = WaveGuideStructure_2;
% end
% for i=1:length(WaveGuideStructure_1)
%     WGS{2} = WaveGuideStructure_1{i};
%     Stmp.S11 = ExtractPortS(S, SeleDim, i, i);
%     Stmp.S12 = ExtractPortS(S, SeleDim, i, SDim);
%     Stmp.S21 = ExtractPortS(S, SeleDim, SDim, i);
%     Stmp.S22 = ExtractPortS(S, SeleDim, SDim, SDim);
%     StmpRN = Renormalize(WGS, Stmp);
%     S = InsertPortS(S, StmpRN.S11, SeleDim, i, i);
%     S = InsertPortS(S, StmpRN.S12, SeleDim, i, SDim);
%     S = InsertPortS(S, StmpRN.S21, SeleDim, SDim, i);
%     S = InsertPortS(S, StmpRN.S21, SeleDim, SDim, SDim);
% end

Stot = S;
