function  [Position] = TwoPortDeviceDraw( TwoPortDevice, explode )
% Draw one TwoPortDevice either compactly or exploded along z.
% ------------------------------------------------------------------------
% [IN]
%   TwoPortDevice - cell array describing the device waveguide structure
%   explode       - 0 no draw, 1 compact draw, 2 exploded draw
%
% [OUT]
%   Position      - z-coordinates of the two external ports
%
ze = 0;
z0 = TwoPortDevice{1}.zo;
if(explode)
    for i=1:length(TwoPortDevice)
        ze = ze + TwoPortDevice{i}.b;
    end
    ze = ze / length(TwoPortDevice);
    if(sign(z0) == -1)
        z0 = -(ze*(length(TwoPortDevice))) + z0;
        if(length(TwoPortDevice) == 1)
            z0 = z0 - ze/2;
        end
    end
end

for i=1:length(TwoPortDevice)-1
    z0 = z0 + ze/2;
    z0 = ShowSegment(TwoPortDevice{i},z0)+ze/2;
    WaveGuideCapDraw(TwoPortDevice{i},TwoPortDevice{i+1},z0);
end

z0 = z0 + ze/2;
z0 = ShowSegment(TwoPortDevice{length(TwoPortDevice)},z0);
axis equal;

Position = z0 + ze/2;
