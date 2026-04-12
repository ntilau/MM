function [ Sf, Sinfo, TwoPortDevices, Nto1Connections,...
    ConnectedPorts, FrequencySweep, Error ] = ...
    MultiPortDevice( TwoPortDevices, Nto1Connections,...
    OpenPorts, ConnectedPorts, FrequencySweep, flag, Options)
% Main entry point for validating, solving, and optionally drawing a
% multiport waveguide device.
% ------------------------------------------------------------------------
% [IN]
%   TwoPortDevices  - cell array of TwoPortDevice definitions
%   Nto1Connections - cell array of N-to-1 connection definitions
%   OpenPorts       - cell array of externally accessible port pointers
%   ConnectedPorts  - cell array of directly connected port pairs
%   FrequencySweep  - frequency sweep definition
%   flag            - 0 no draw, 1 compact draw, 2 exploded draw
%   Options         - solver and drawing options struct
%
% [OUT]
%   Sf, Sinfo       - frequency-dependent condensed GSM data
%   TwoPortDevices, Nto1Connections, ConnectedPorts, FrequencySweep -
%                    validated and enriched inputs
%   Error           - Error struct describing validation or solve failures
%
Sf    = {};
Sinfo = {};
[ TwoPortDevices, Nto1Connections, OpenPorts, ConnectedPorts, ...
    FrequencySweep, Symmetry, Topology, Error ] = ...
    MultiPortDeviceValidate( TwoPortDevices, Nto1Connections, OpenPorts,...
    ConnectedPorts, FrequencySweep, flag, Options);
if (DumpError('|  +--+>',Error))
    return;
end

%
% 2 - PERFORM COMPUTATIONS
%
[ Sf, Sinfo, Error, TwoPortDevices, Nto1Connections, ConnectedPorts] = ...
    MultiPortDeviceSolve( TwoPortDevices, Nto1Connections, OpenPorts,...
    ConnectedPorts, FrequencySweep, Symmetry, Topology, Options);
if (DumpError('|  +--+>',Error))
    return;
end
