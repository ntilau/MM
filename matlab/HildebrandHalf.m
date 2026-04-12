clear; clc
% MM_LIB_PATH_BOOTSTRAP
projectRoot = fileparts(mfilename('fullpath'));
addpath(fullfile(projectRoot, 'lib'));

%%  HILDEBRAND COUPLER DIMENSIONS
a = 0.75*0.0254;
b = 0.375*0.0254;
length = 1*0.0254;
t = 0.2*0.0254;
R = 0.040*0.0254;
L1 = 1.282*0.0254;
L2 = 0.572*0.0254;
L3 = 0.172*0.0254;
Lgap = 0.838*0.0254;
W1 = 1.4*0.0254;
W2 = 1.142*0.0254;
W3 = 1.06*0.0254;

%% Define Frequency Sweep
FS.start = 13*10^9;
FS.end   = 15*10^9;
FS.N     = 51;

%% Number of modes for calculation accuracy
Lambda = 299792458/FS.end;
Nmodes = floor(10*sqrt(W1*b)/Lambda + 0.5);
% Nmodes = 16;

%% BUILDING STRUCTURE
WGS{1}.D{1}.a = a;
WGS{1}.D{1}.b = b;
WGS{1}.D{1}.Nmodes = Nmodes;
WGS{1}.D{1}.l = length;
WGS{1}.D{1}.xo = a/2+t/2;
WGS{1}.D{1}.yo = 0.0;
WGS{1}.D{1}.zo = -((L1-Lgap)/2 + (WGS{1}.D{1}.l));

WGS{1}.D{2}.a = W1/2-t/2;
WGS{1}.D{2}.b = b;
WGS{1}.D{2}.Nmodes = Nmodes;
WGS{1}.D{2}.l = (L1-Lgap)/2;
WGS{1}.D{2}.xo = (W1/2+t/2)/2;
WGS{1}.D{2}.yo = 0.0;
WGS{1}.D{2}.zo = -(L1-Lgap)/2;

WGS{2}.D{1}.a = a;
WGS{2}.D{1}.b = b;
WGS{2}.D{1}.Nmodes = Nmodes;
WGS{2}.D{1}.l = length;
WGS{2}.D{1}.xo = -(a/2+t/2);
WGS{2}.D{1}.yo = 0;
WGS{2}.D{1}.zo = -((L1-Lgap)/2 + (WGS{2}.D{1}.l));

WGS{2}.D{2}.a = W1/2-t/2;
WGS{2}.D{2}.b = b;
WGS{2}.D{2}.Nmodes = Nmodes;
WGS{2}.D{2}.l = (L1-Lgap)/2;
WGS{2}.D{2}.xo = -(W1/2+t/2)/2;
WGS{2}.D{2}.yo = 0.0;
WGS{2}.D{2}.zo = -(L1-Lgap)/2;

WGS{3}.D{1}.a = W1;
WGS{3}.D{1}.b = b;
WGS{3}.D{1}.Nmodes = Nmodes;
WGS{3}.D{1}.l = (Lgap-L2)/2;
WGS{3}.D{1}.xo = 0.0;
WGS{3}.D{1}.yo = 0.0;
WGS{3}.D{1}.zo = 0.0;

WGS{3}.D{2}.a = W2;
WGS{3}.D{2}.b = b;
WGS{3}.D{2}.Nmodes = Nmodes;
WGS{3}.D{2}.l = (L2-L3)/2;
WGS{3}.D{2}.xo = 0.0;
WGS{3}.D{2}.yo = 0.0;
WGS{3}.D{2}.zo = WGS{3}.D{1}.l;

WGS{3}.D{3}.a = W3;
WGS{3}.D{3}.b = b;
WGS{3}.D{3}.Nmodes = Nmodes;
WGS{3}.D{3}.l = L3/2;
WGS{3}.D{3}.xo = 0.0;
WGS{3}.D{3}.yo = 0.0;
WGS{3}.D{3}.zo = WGS{3}.D{2}.l;

%% Define Nto1Connections
%% Note that the device solve only Nto1 or 1toM. NtoM can be done by
%% cascading with zero-length segment
Nto1{1}.SideOne{1}.TwoPortDeviceIndex = 1;
Nto1{1}.SideOne{1}.TwoPortDevicePort = 2;
Nto1{1}.SideOne{2}.TwoPortDeviceIndex = 2;
Nto1{1}.SideOne{2}.TwoPortDevicePort = 2;
Nto1{1}.SideTwo{1}.TwoPortDeviceIndex = 3;
Nto1{1}.SideTwo{1}.TwoPortDevicePort = 1;
Nto1{1}.zo = 0;

%% Define Open Ports
OpenPorts{1}.TwoPortDeviceIndex = 1;
OpenPorts{1}.TwoPortDevicePort = 1;
OpenPorts{2}.TwoPortDeviceIndex = 2;
OpenPorts{2}.TwoPortDevicePort = 1;
OpenPorts{3}.TwoPortDeviceIndex = 3;
OpenPorts{3}.TwoPortDevicePort = 2;

%% Define Connected Ports
ConnectedPorts = {};%{1}.TwoPortDeviceIndex = [3, 4];
%ConnectedPorts{1}.TwoPortDevicePort = [2, 1];

%% Define Options
Options.DeviceSymmetry.Use = 1;
Options.DeviceSymmetry.Side = 2;
Options.Connections = 0;

%% Enable Device plot
% 2 = Explodes parts; 1 = shows united parts; 0 = no plot
DevicePlotType = 2;

%% Solve
[ Sf, Sinfo, WGS, Nto1, ConnectedPorts, FS, Error ] = ...
    MultiPortDevice( WGS, Nto1, OpenPorts, ConnectedPorts, FS, DevicePlotType, Options);

f = FS.f;

% Plotting S21 S31 (In Hildebrand's paper 2 and 3 are the half power
% outputs)
ModeStruct = {{3,1,'h',1,0,'h',1,0,'md'};...
              {4,1,'h',1,0,'h',1,0,'md'};...
              };
GSMDraw(f,Sf,Sinfo,ModeStruct,1);
axis([min(f) max(f) -4  -2]);
grid on;
%legend off;

% Plotting S11 S41 reflection coefficients
ModeStruct = {{1,1,'h',1,0,'h',1,0,'md'};...
              {2,1,'h',1,0,'h',1,0,'md'};...
              };
GSMDraw(f,Sf,Sinfo,ModeStruct,1);
axis([min(f) max(f) -60  0]);
grid on;
% legend off;

% Plotting relative phases between S21 and S31
ModeStruct = {{3,1,'h',1,0,'h',1,0,'md'};...
              {4,1,'h',1,0,'h',1,0,'md'};...
              };
% ModeStruct = {{1,1,'h',1,0,'h',1,0,'md'};...
%               {2,1,'h',1,0,'h',1,0,'md'};...
%               };
RelativePhaseDraw(f,Sf,Sinfo,ModeStruct,1);
axis([min(f) max(f) 89 91]);
grid on;
% legend off;
