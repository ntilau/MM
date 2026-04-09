function Nto1DeviceDraw( TwoPortDevices, Nto1Connection )
% Draw one N-to-1 connection and the attached TwoPortDevices.
% ------------------------------------------------------------------------
% [IN]
%   TwoPortDevices  - cell array of all TwoPortDevice structs
%   Nto1Connection  - Nto1Connection struct to draw
%
% [OUT]
%   (none)          - adds the geometry to the current figure
%
zo = Nto1Connection.zo;

nbrPortsSideOne = length(Nto1Connection.SideOne);
nbrPortsSideTwo = length(Nto1Connection.SideTwo);

if (nbrPortsSideOne > nbrPortsSideTwo)
    for i=1:nbrPortsSideOne
        d1(i) = Nto1Connection.SideOne{i}.TwoPortDeviceIndex;
        p1(i) = Nto1Connection.SideOne{i}.TwoPortDevicePort;

        TwoPortSegment = TwoPortDeviceGetPortSegment(TwoPortDevices{d1(i)}.D,p1(i));
        [a(i),b(i),xo(i),yo(i)] =...
            WaveGuideSegmentGetCrossSection(TwoPortSegment);

        xss(2*i-1) = xo(i) - a(i)/2;
        xss(2*i)   = xo(i) + a(i)/2;
        yss(2*i-1) = yo(i) - b(i)/2;
        yss(2*i)   = yo(i) + b(i)/2;
    end

    for i=1:nbrPortsSideTwo
        d2(i) = Nto1Connection.SideTwo{i}.TwoPortDeviceIndex;
        p2(i) = Nto1Connection.SideTwo{i}.TwoPortDevicePort;
        TwoPortSegment = TwoPortDeviceGetPortSegment(TwoPortDevices{d2(i)}.D,p2(i));
        [xss(2*nbrPortsSideOne+2*i-1),xss(2*nbrPortsSideOne+2*i),...
            yss(2*nbrPortsSideOne+2*i-1),yss(2*nbrPortsSideOne+2*i)] =...
            WaveGuideSegmentGetBounding(TwoPortSegment);
    end
else
    for i=1:nbrPortsSideTwo
        d2(i) = Nto1Connection.SideTwo{i}.TwoPortDeviceIndex;
        p2(i) = Nto1Connection.SideTwo{i}.TwoPortDevicePort;
        TwoPortSegment = TwoPortDeviceGetPortSegment(TwoPortDevices{d2(i)}.D,p2(i));
        [a(i),b(i),xo(i),yo(i)] =...
            WaveGuideSegmentGetCrossSection(TwoPortSegment);

        xss(2*i-1) = xo(i) - a(i)/2;
        xss(2*i)   = xo(i) + a(i)/2;
        yss(2*i-1) = yo(i) - b(i)/2;
        yss(2*i)   = yo(i) + b(i)/2;
    end
    for i=1:nbrPortsSideOne
        d1(i) = Nto1Connection.SideOne{i}.TwoPortDeviceIndex;
        p1(i) = Nto1Connection.SideOne{i}.TwoPortDevicePort;

        TwoPortSegment = TwoPortDeviceGetPortSegment(TwoPortDevices{d1(i)}.D,p1(i));
        [xss(2*nbrPortsSideTwo+2*i-1),xss(2*nbrPortsSideTwo+2*i),...
            yss(2*nbrPortsSideTwo+2*i-1),yss(2*nbrPortsSideTwo+2*i)] =...
            WaveGuideSegmentGetBounding(TwoPortSegment);

    end
end

xs = sort(unique(xss));
ys = sort(unique(yss));

for i=1:(length(xs)-1)
    for j=1:(length(ys)-1)
        if (NotInRect((xs(i)+xs(i+1))/2,(ys(j)+ys(j+1))/2,a,b,xo,yo))
            patch([xs(i),xs(i+1),xs(i+1),xs(i)],...
                  [ys(j),ys(j),ys(j+1),ys(j+1)],...
                  [zo, zo, zo, zo],'r','EdgeColor','k');
        else
            patch([xs(i),xs(i+1),xs(i+1),xs(i)],...
                  [ys(j),ys(j),ys(j+1),ys(j+1)],...
                  [zo, zo, zo, zo],'g','EdgeColor','none');
        end
    end
end

for i=1:length(Nto1Connection.SideOne)
    line([xss(2*i-1),xss(2*i),xss(2*i),xss(2*i-1),xss(2*i-1)],...
        [yss(2*i-1),yss(2*i-1),yss(2*i),yss(2*i),yss(2*i-1)],...
        [zo,zo,zo,zo,zo],'Color','k');
end

for i=1:length(Nto1Connection.SideTwo)
    line([xss(2*i-1),xss(2*i),xss(2*i),xss(2*i-1),xss(2*i-1)],...
        [yss(2*i-1),yss(2*i-1),yss(2*i),yss(2*i),yss(2*i-1)],...
        [zo,zo,zo,zo,zo],'Color','k');
end

axis equal;
