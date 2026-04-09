function [Sout] = CondenseGSM( S, SeleDim, Topology)
% Condense a global scattering matrix by eliminating internally connected
% ports according to the supplied topology.
% ------------------------------------------------------------------------
% [IN]
%   S        - global GSM arranged as block submatrices
%   SeleDim  - size of each block scattering submatrix
%   Topology - topology struct containing the port pairs to condense
%
% [OUT]
%   Sout     - condensed GSM containing only the retained ports
%
nbrTwoPortDevices = (size(S)./SeleDim);
if(nbrTwoPortDevices(1) == nbrTwoPortDevices(2))
    nbrPorts = nbrTwoPortDevices(1);
else
    warning('Matrix is not condensable!!');
    return;
end

TwoPortDevicesLength = length(Topology.PortToCondense.a);

for Loop=1:TwoPortDevicesLength

    indexm = Topology.PortToCondense.a(Loop);
    indexn = Topology.PortToCondense.b(Loop);

    %% building CDEF
    C = ExtractPortS(S, SeleDim, indexm, indexm);
    D = ExtractPortS(S, SeleDim, indexm, indexn) - eye(SeleDim);
    E = ExtractPortS(S, SeleDim, indexn, indexm) - eye(SeleDim);
    F = ExtractPortS(S, SeleDim, indexn, indexn);

    %% blockwise inversion
    warning off all
    % Disables computation warnings
    % They were found some boring warnings when computing the inverse of
    % C D E F matrices

    if(det(C) ~= 0)
        Xi = E*inv(C);
        Phi.nm = inv(Xi*D-F)*Xi;
        Phi.nn = -inv(Xi*D-F);
    elseif(det(E) ~= 0)
        Xi = C*inv(E);
        Phi.nm = inv(D-Xi*F);
        Phi.nn = -inv(D-Xi*F)*Xi;
    end
    if(det(F) ~= 0)
        Psi = D*inv(F);
        Phi.mm = -inv(Psi*E-C);
        Phi.mn = inv(Psi*E-C)*Psi;
    elseif(det(D) ~= 0)
        Psi = F*inv(D);
        Phi.mm = -inv(E-Psi*C)*Psi;
        Phi.mn = inv(E-Psi*C);
    end

    warning on all
    % Re-enables all warnings for following operations (in the calling
    % function)

    %% Rebuilding Matrix
    LoopS=[];
    for i=1:nbrPorts
        if((i ~= indexm) && (i ~= indexn))
            for j=1:nbrPorts
                if( (j ~= indexm) && (j ~= indexn) )
                    Stmp = ExtractPortS(S, SeleDim, i, j) - ...
                        (ExtractPortS(S, SeleDim, i, indexm) * ...
                        (Phi.mm*ExtractPortS(S, SeleDim, indexm, j) + ...
                        Phi.mn*ExtractPortS(S, SeleDim, indexn, j))) - ...
                        (ExtractPortS(S, SeleDim, i, indexn) * ...
                        (Phi.nm*ExtractPortS(S, SeleDim, indexm, j) + ...
                        Phi.nn*ExtractPortS(S, SeleDim, indexn, j)));
                    LoopS = InsertPortS(LoopS, Stmp, SeleDim, i, j);
                end
            end
        end
    end
    Loop = Loop+1;
    S = LoopS;
end

%% Extract OpenPort related scattering matrices
Sout=[];
for i=1:length(Topology.OpenPorts)
    for j=1:length(Topology.OpenPorts)
        Stmp = ExtractPortS(S, SeleDim, Topology.OpenPorts(i), Topology.OpenPorts(j));
        Sout = InsertPortS(Sout, Stmp, SeleDim, i, j);
    end
end
