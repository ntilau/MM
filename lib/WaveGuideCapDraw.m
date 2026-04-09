function WaveGuideCapDraw( WaveGuideSegmentA, WaveGuideSegmentB, z0 )
% Draw the transverse cap joining two waveguide cross-sections in 3-D.
% ------------------------------------------------------------------------
% [IN]
%   WaveGuideSegmentA - first waveguide segment struct
%   WaveGuideSegmentB - second waveguide segment struct
%   z0                - z-position of the cap plane [m]
%
% [OUT]
%   (none)            - adds patches to the current figure
%
Conn(:,:,1) = [[1,2,2,1];...
               [2,3,3,2];...
               [3,4,4,3];...
               [1,2,2,1];...
               [3,4,4,3];...
               [1,2,2,1];...
               [2,3,3,2];...
               [3,4,4,3];...
               ];

Conn(:,:,2) = [[1,1,2,2];...
               [1,1,2,2];...
               [1,1,2,2];...
               [2,2,3,3];...
               [2,2,3,3];...
               [3,3,4,4];...
               [3,3,4,4];...
               [3,3,4,4];...
               ];

a1  = WaveGuideSegmentA.a;
b1  = WaveGuideSegmentA.b;
xo1 = WaveGuideSegmentA.xo;
yo1 = WaveGuideSegmentA.yo;

a2  = WaveGuideSegmentB.a;
b2  = WaveGuideSegmentB.b;
xo2 = WaveGuideSegmentB.xo;
yo2 = WaveGuideSegmentB.yo;

x = [xo1 - a1/2, xo2 - a2/2, xo2 + a2/2, xo1 + a1/2];
y = [yo1 - b1/2, yo2 - b2/2, yo2 + b2/2, yo1 + b1/2];

for i=1:8
    xp = [x(Conn(i,1,1)),x(Conn(i,2,1)),x(Conn(i,3,1)),x(Conn(i,4,1))];
    yp = [y(Conn(i,1,2)),y(Conn(i,2,2)),y(Conn(i,3,2)),y(Conn(i,4,2))];

    if (min(xp) >= xo1-a1/2-eps & max(xp) <= xo1+a1/2 +eps & ...
        min(yp) >= yo1-b1/2-eps & max(yp) <= yo1+b1/2 + eps)
        in1 = 1;
    else
        in1 = 0;
    end

    if (min(xp) >= xo2-a2/2-eps & max(xp) <= xo2+a2/2+eps  & ...
        min(yp) >= yo2-b2/2-eps & max(yp) <= yo2+b2/2+eps)
        in2 = 1;
    else
        in2 = 0;
    end

    if (in1 | in2)
        patch(xp,yp,[z0, z0, z0, z0],'r')
    end
end

