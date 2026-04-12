clear; clc
% MM_LIB_PATH_BOOTSTRAP
projectRoot = fileparts(mfilename('fullpath'));
addpath(fullfile(projectRoot, 'lib'));

%% Constants
c = 299792458;

%%  RIBLET'S COUPLER DIMENSIONS
a = 0.75*0.0254;
b = 0.375*0.0254;
Lg = c/(10*10^9*1.25);
length = 2*Lg;
Lgap = 5*Lg/4;
t = Lg/500;

%% Define Frequency Sweep
FS.start = 7*10^9;
FS.end   = 15*10^9;
FS.N     = 251;

%% Number of modes for calculation accuracy
Lambda = c/FS.end;
Nmodes = floor(10*sqrt((2*a+t)*b)/Lambda + 0.5);
%Nmodes = 12;

%% BUILDING STRUCTURE
WGS{1}.D{1}.a = a;
WGS{1}.D{1}.b = b;
WGS{1}.D{1}.Nmodes = Nmodes;
WGS{1}.D{1}.l = length;
WGS{1}.D{1}.xo = a/2+t/2;
WGS{1}.D{1}.yo = 0.0;
WGS{1}.D{1}.zo = 0;

WGS{2}.D{1}.a = a;
WGS{2}.D{1}.b = b;
WGS{2}.D{1}.Nmodes = Nmodes;
WGS{2}.D{1}.l = length;
WGS{2}.D{1}.xo = -(a/2+t/2);
WGS{2}.D{1}.yo = 0;
WGS{2}.D{1}.zo = 0;

WGS{3}.D{1}.a = 2*a+t;
WGS{3}.D{1}.b = b;
WGS{3}.D{1}.Nmodes = Nmodes;
WGS{3}.D{1}.l = Lgap;
WGS{3}.D{1}.xo = 0;
WGS{3}.D{1}.yo = 0;
WGS{3}.D{1}.zo = 0;

WGS{4}.D{1}.a = a;
WGS{4}.D{1}.b = b;
WGS{4}.D{1}.Nmodes = Nmodes;
WGS{4}.D{1}.l = length;
WGS{4}.D{1}.xo = a/2+t/2;
WGS{4}.D{1}.yo = 0.0;
WGS{4}.D{1}.zo = 0;

WGS{5}.D{1}.a = a;
WGS{5}.D{1}.b = b;
WGS{5}.D{1}.Nmodes = Nmodes;
WGS{5}.D{1}.l = length;
WGS{5}.D{1}.xo = -(a/2+t/2);
WGS{5}.D{1}.yo = 0.0;
WGS{5}.D{1}.zo = 0;

%% Define Nto1Connections
Nto1{1}.SideOne{1}.TwoPortDeviceIndex = 1;
Nto1{1}.SideOne{1}.TwoPortDevicePort = 2;
Nto1{1}.SideOne{2}.TwoPortDeviceIndex = 2;
Nto1{1}.SideOne{2}.TwoPortDevicePort = 2;
Nto1{1}.SideTwo{1}.TwoPortDeviceIndex = 3;
Nto1{1}.SideTwo{1}.TwoPortDevicePort = 1;
Nto1{1}.zo = 0;

Nto1{2}.SideOne{1}.TwoPortDeviceIndex = 3;
Nto1{2}.SideOne{1}.TwoPortDevicePort = 2;
Nto1{2}.SideTwo{1}.TwoPortDeviceIndex = 4;
Nto1{2}.SideTwo{1}.TwoPortDevicePort = 1;
Nto1{2}.SideTwo{2}.TwoPortDeviceIndex = 5;
Nto1{2}.SideTwo{2}.TwoPortDevicePort = 1;
Nto1{2}.zo = 0.07;

%% Define Open Ports
OpenPorts{1}.TwoPortDeviceIndex = 1;
OpenPorts{1}.TwoPortDevicePort = 1;
OpenPorts{2}.TwoPortDeviceIndex = 2;
OpenPorts{2}.TwoPortDevicePort = 1;
OpenPorts{3}.TwoPortDeviceIndex = 4;
OpenPorts{3}.TwoPortDevicePort = 2;
OpenPorts{4}.TwoPortDeviceIndex = 5;
OpenPorts{4}.TwoPortDevicePort = 2;

%% Define Connected Ports
ConnectedPorts = {};

%% Define Options
% Overall device's symmetry (enable and side)
Options.DeviceSymmetry.Use = 0;
Options.DeviceSymmetry.Side = 2;

%% Enable Device plot
% 2 = Explodes parts; 1 = shows united parts; 0 = no plot
DevicePlotType = 2;

%% Solve
[ Sf, Sinfo, WGS, Nto1, ConnectedPorts, FS, Error ] = ...
    MultiPortDevice( WGS, Nto1, OpenPorts, ConnectedPorts, FS, DevicePlotType, Options);

f = FS.f;

% Plotting S21 S31 (In Hildebrand's paper 2 and 3 are the half power
% outputs)
ModeStruct = {{1,1,'h',1,0,'h',1,0,'md'};...
              {2,1,'h',1,0,'h',1,0,'md'};...
              {3,1,'h',1,0,'h',1,0,'md'};...
              {4,1,'h',1,0,'h',1,0,'md'};...
              };
GSMDraw(f,Sf,Sinfo,ModeStruct,1);
axis([min(f) max(f) -20  0]);
ylabel('Amplitude [dB]','FontSize',12)
xlabel('Frequency [Hz]','FontSize',12)
grid on;
% legend off;

% Plotting S11 S41 reflection coefficients
% ModeStruct = {{1,1,'h',1,0,'h',1,0,'md'};...
%               {2,1,'h',1,0,'h',1,0,'md'};...
%               };
% GSMDraw(f,Sf,Sinfo,ModeStruct,1);
% axis([min(f) max(f) -40  0]);
% ylabel('Amplitude [dB]','FontSize',12)
% xlabel('Frequency [Hz]','FontSize',12)
% grid on;
% legend off;

% Plotting relative phases between S21 and S31
ModeStruct = {{3,1,'h',1,0,'h',1,0,'md'};...
              {4,1,'h',1,0,'h',1,0,'md'};...
              };
RelativePhaseDraw(f,Sf,Sinfo,ModeStruct,1);
% axis([min(f) max(f) 89 91]);
ylabel('Phase [deg]','FontSize',12)
xlabel('Frequency [Hz]','FontSize',12)
grid on;
% legend off;