function [ Topology, Error ] = ...
    MultiPortDeviceTopology( TwoPortDevices, Nto1Connections, OpenPorts)
% Build the local port-numbering and condensation topology used by the
% MultiPortDevice solver.
% ------------------------------------------------------------------------
% [IN]
%   TwoPortDevices  - cell array of TwoPortDevice structs
%   Nto1Connections - cell array of Nto1Connection structs
%   OpenPorts       - cell array of open port pointers
%
% [OUT]
%   Topology        - port numbering and condensation maps
%   Error           - Error struct
%
Error = struct();
nbrTwoPortDevices = length(TwoPortDevices);
nbrNto1Connections = length(Nto1Connections);
nbrOpenPorts = length(OpenPorts);

%% For each Nto1Connection defines a local topology
for index=1:nbrNto1Connections

    nbrSideOne = length(Nto1Connections{index}.SideOne);
    nbrSideTwo = length(Nto1Connections{index}.SideTwo);
    Topology.Nto1{index}.Dimensions = nbrSideOne + nbrSideTwo;
    Topology.Nto1{index}.PortToCondense.a = [];
    Topology.Nto1{index}.PortToCondense.b = [];
    Topology.Nto1{index}.OpenPorts = [];
    indexforNto1Ports = (nbrSideOne + nbrSideTwo)*2;

    % Checks if forward solving or reversed
    if(nbrSideOne > 1 && nbrSideTwo == 1)
        Topology.Nto1{index}.nFurcation = 0;

        for i=1:nbrSideOne
            Topology.Nto1{index}.PortToCondense.a = ...
                horzcat(Topology.Nto1{index}.PortToCondense.a, 2*i);
            Topology.Nto1{index}.PortToCondense.b = ...
                horzcat(Topology.Nto1{index}.PortToCondense.b, indexforNto1Ports+i);
            Topology.Nto1{index}.OpenPorts = horzcat(Topology.Nto1{index}.OpenPorts, 2*i-1);
        end
        for i=1:nbrSideTwo
            Topology.Nto1{index}.PortToCondense.a = horzcat(Topology.Nto1{index}.PortToCondense.a, ...
                indexforNto1Ports + nbrSideOne + i);
            Topology.Nto1{index}.PortToCondense.b = horzcat(Topology.Nto1{index}.PortToCondense.b, ...
                nbrSideOne*2 + i*2-1);
            Topology.Nto1{index}.OpenPorts = horzcat(Topology.Nto1{index}.OpenPorts, ...
                (nbrSideOne+nbrSideTwo)*2);
        end

    elseif(nbrSideTwo > 1 && nbrSideOne == 1)
        Topology.Nto1{index}.nFurcation = 1; % splitting from one to n ports

        for i=1:nbrSideTwo
            Topology.Nto1{index}.PortToCondense.a = ...
                horzcat(Topology.Nto1{index}.PortToCondense.a, 2*i);
            Topology.Nto1{index}.PortToCondense.b = ...
                horzcat(Topology.Nto1{index}.PortToCondense.b, indexforNto1Ports+i);
            Topology.Nto1{index}.OpenPorts = horzcat(Topology.Nto1{index}.OpenPorts, 2*i-1);
        end
        for i=1:nbrSideOne
            Topology.Nto1{index}.PortToCondense.a = horzcat(Topology.Nto1{index}.PortToCondense.a, ...
                indexforNto1Ports + nbrSideTwo + i);
            Topology.Nto1{index}.PortToCondense.b = horzcat(Topology.Nto1{index}.PortToCondense.b, ...
                nbrSideTwo*2 + i*2-1);
            Topology.Nto1{index}.OpenPorts = horzcat(Topology.Nto1{index}.OpenPorts, ...
                (nbrSideOne+nbrSideTwo)*2);
        end
    else
        Error.fatal = 'Nto1 Connection has not been defined correctly';
        return
    end
end
