clear; clc
% MM_LIB_PATH_BOOTSTRAP
projectRoot = fileparts(mfilename('fullpath'));
addpath(fullfile(projectRoot, 'lib'));

c       = 299792458;  % speed of light [m/s]
a       = 0.01905;    % waveguide full width [m]
b       = 0.009525;   % waveguide full height [m]
l       = 0.01;       % section length [m]
Nmodes  = floor(10*sqrt(2*a*b)/(c/(25*10^9)) + 0.5);

%% Defines the TwoPortDevices
WGS1{1}.D{1}.a = a;
WGS1{1}.D{1}.b = b;
WGS1{1}.D{1}.Nmodes = Nmodes;
WGS1{1}.D{1}.l = l;
WGS1{1}.D{1}.xo = -0.01905/2-0.001;
WGS1{1}.D{1}.yo = 0;
WGS1{1}.D{1}.zo = -l;

WGS1{2}.D{1}.a = a;
WGS1{2}.D{1}.b = b;
WGS1{2}.D{1}.Nmodes = Nmodes;
WGS1{2}.D{1}.l = l;
WGS1{2}.D{1}.xo = 0.01905/2+0.001;
WGS1{2}.D{1}.yo = 0;
WGS1{2}.D{1}.zo = -l;

WGS1{3}.D{1}.a = 2*a;
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
FS.N     = 61;

%% Defines Connected Ports
ConnectedPorts = {};

%% Defines Options
Options.DeviceSymmetry.Use = 0;
Options.DeviceSymmetry.Side = 2;
Options.Connections = 0;

%% Enable Device plot
% 2 = Explodes parts; 1 = shows united parts; 0 = no plot
DevicePlotType = 0;

%% Calculus
[ Sf, Sinfo, WGS1, Nto1, ConnectedPorts, FS, Error ] = ...
    MultiPortDevice( WGS1, Nto1, OpenPorts, ConnectedPorts, FS, DevicePlotType, Options);

%% Results
f = FS.f;

ModeStruct = {{1,1,'h',1,0,'h',1,0,'md'};...
              {2,1,'h',1,0,'h',1,0,'md'};...
              {3,1,'h',1,0,'h',1,0,'md'};...
              {3,1,'h',2,0,'h',1,0,'md'};...
              {3,1,'h',3,0,'h',1,0,'md'};...
              {3,1,'h',4,0,'h',1,0,'md'};...
              {3,1,'h',5,0,'h',1,0,'md'};...
              };

GSMDraw(f,Sf,Sinfo,ModeStruct,1);
ylabel('Amplitude [dB]','FontSize',12);
xlabel('Frequency [Hz]','FontSize',12);

axis([min(f) max(f) -70 0]);

% Load reference full-wave benchmark data when available.
load('BifurcationH.mat');
hold on;
plot(HFSSh(:,1)*10^9,HFSSh(:,2),'k:',...
    HFSSh(:,1)*10^9,HFSSh(:,3),'r:',...
    HFSSh(:,1)*10^9,HFSSh(:,4),'g:',...
    HFSSh(:,1)*10^9,HFSSh(:,5),'b:',...
    HFSSh(:,1)*10^9,HFSSh(:,6),'y:',...
    HFSSh(:,1)*10^9,HFSSh(:,7),'m:',...
    HFSSh(:,1)*10^9,HFSSh(:,8),'c:'...
    );