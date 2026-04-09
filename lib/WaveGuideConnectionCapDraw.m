function WaveGuideConnectionCapDraw( TwoPortDevice, Color,  z0 )
% Draw the cap patch corresponding to a TwoPortDevice connection face.
% ------------------------------------------------------------------------
% [IN]
%   TwoPortDevice - segment or device geometry to draw
%   Color         - patch color specification
%   z0            - z-position of the cap plane [m]
%
% [OUT]
%   (none)        - adds patches to the current figure
%
a1  = TwoPortDevice.a;
b1  = TwoPortDevice.b;
xo1 = TwoPortDevice.xo;
yo1 = TwoPortDevice.yo;

x = [xo1 - a1/2, xo1 + a1/2, xo1 + a1/2, xo1 - a1/2];
y = [yo1 - b1/2, yo1 - b1/2, yo1 + b1/2, yo1 + b1/2];

patch(x,y,[z0, z0, z0, z0], Color, 'EdgeColor', 'k', 'Marker', '*');
