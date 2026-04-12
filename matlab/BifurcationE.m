clear; clc
% MM_LIB_PATH_BOOTSTRAP
projectRoot = fileparts(mfilename('fullpath'));
addpath(fullfile(projectRoot, 'lib'));

a      = 0.01905;     % waveguide full width [m]
b      = 0.009525;    % waveguide full height [m]
Nmodes = 12;          % number of modes
l      = 0.01;        % section length [m]

%% Defines the TwoPortDevices
WGS1{1}.D{1}.a = a;
WGS1{1}.D{1}.b = b/2;
WGS1{1}.D{1}.Nmodes = Nmodes;
WGS1{1}.D{1}.l = l;
WGS1{1}.D{1}.xo = 0;
WGS1{1}.D{1}.yo = -0.009525/4;
WGS1{1}.D{1}.zo = -l;

WGS1{2}.D{1}.a = a;
WGS1{2}.D{1}.b = b/2;
WGS1{2}.D{1}.Nmodes = Nmodes;
WGS1{2}.D{1}.l = l;
WGS1{2}.D{1}.xo = 0;
WGS1{2}.D{1}.yo = 0.009525/4;
WGS1{2}.D{1}.zo = -l;

WGS1{3}.D{1}.a = a;
WGS1{3}.D{1}.b = b;
WGS1{3}.D{1}.Nmodes = Nmodes;
WGS1{3}.D{1}.l = l;
WGS1{3}.D{1}.xo = 0;
WGS1{3}.D{1}.yo = 0;
WGS1{3}.D{1}.zo = 0;

%% Defines the Nto1Connections
Nto1{1}.SideOne{1}.TwoPortDeviceIndex = 1;
Nto1{1}.SideOne{1}.TwoPortDevicePort = 2;
Nto1{1}.SideOne{2}.TwoPortDeviceIndex = 2;
Nto1{1}.SideOne{2}.TwoPortDevicePort = 2;
Nto1{1}.SideTwo{1}.TwoPortDeviceIndex = 3;
Nto1{1}.SideTwo{1}.TwoPortDevicePort = 1;
Nto1{1}.zo = 0;

%% Defines the Opened Ports
OpenPorts{1}.TwoPortDeviceIndex = 1;
OpenPorts{1}.TwoPortDevicePort = 1;
OpenPorts{2}.TwoPortDeviceIndex = 2;
OpenPorts{2}.TwoPortDevicePort = 1;
OpenPorts{3}.TwoPortDeviceIndex = 3;
OpenPorts{3}.TwoPortDevicePort = 2;

%% Defines the Frequency Sweep
FS.start = 10*10^9;
FS.end   = 25*10^9;
FS.N     = 51;

%% Defines Connected Ports
ConnectedPorts = {};

%% Enable structure'symmetry solving
% Repeats once the Nto1 device by trasversal mirroring at the open ports at
% side one (N ports) or side two (single port). Use = 0 -> no, = 1 -> yes
Options.DeviceSymmetry.Use = 0;
Options.DeviceSymmetry.Side = 2;
Options.Connections = 0;

%% Enable Device plot
% 2 = Explodes parts; 1 = shows united parts; 0 = no plot
DevicePlotType = 2;

%% Calculus
[ Sf, Sinfo, WGS1, Nto1, ConnectedPorts, FS, Error ] = ...
    MultiPortDevice( WGS1, Nto1, OpenPorts, ConnectedPorts, FS, DevicePlotType, Options);

%% Results
f = FS.f;

ModeStruct = {{1,1,'h',1,0,'h',1,0,'md'};...
              {2,1,'h',1,0,'h',1,0,'md'};...
              {3,1,'h',1,0,'h',1,0,'md'};...
              {3,1,'h',1,1,'h',1,0,'md'};...
              {3,1,'h',1,3,'h',1,0,'md'};...
              {3,1,'e',1,1,'h',1,0,'md'};...
              {3,1,'e',1,3,'h',1,0,'md'};...
              };

GSMDraw(f,Sf,Sinfo,ModeStruct,1);

% Load reference full-wave benchmark data when available.
HFSSe = csvread(fullfile(projectRoot, '..', 'HFSSe.csv'));
hold on;
plot(HFSSe(:,1)*10^9,HFSSe(:,2),'k:',...
    HFSSe(:,1)*10^9,HFSSe(:,3),'r:',...
    HFSSe(:,1)*10^9,HFSSe(:,4),'g:',...
    HFSSe(:,1)*10^9,HFSSe(:,5),'m:',...
    HFSSe(:,1)*10^9,HFSSe(:,6),'b:'...
    );