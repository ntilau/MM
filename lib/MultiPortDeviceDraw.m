function MultiPortDeviceDraw(TwoPortDevices, Nto1Connections, flag, Options)
% Draw every TwoPortDevice and N-to-1 connection belonging to a
% MultiPortDevice definition.
% ------------------------------------------------------------------------
% [IN]
%   TwoPortDevices  - cell array of TwoPortDevice structs
%   Nto1Connections - cell array of Nto1Connection structs
%   flag            - 0 compact draw, 1 exploded draw
%   Options         - draw options struct
%
% [OUT]
%   (none)          - opens and populates a figure
%
figure

Position = 0;
EndPosition = 0;

for index=1:length(Nto1Connections)
    CurrentWaveGuideStructure = {};
    for i=1:length(Nto1Connections{index}.SideOne)
        d1(i) = Nto1Connections{index}.SideOne{i}.TwoPortDeviceIndex;
        if(index == 1 && Options.DeviceSymmetry.Use == 1 && Options.DeviceSymmetry.Side == 1)
            SymmetryConnection = TwoPortDeviceGetPortSegment(TwoPortDevices{d1(i)}.D,1);
            WaveGuideConnectionCapDraw(SymmetryConnection, 'm', 0);
        end
        if(index > 1)
            for gamma=1:index-1
                for alpha=1:length(Nto1Connections{gamma}.SideTwo)
                    for beta=1:length(Nto1Connections{index}.SideOne)
                        if (Nto1Connections{gamma}.SideTwo{alpha}.TwoPortDeviceIndex == ...
                                Nto1Connections{index}.SideOne{beta}.TwoPortDeviceIndex)
                            CurrentWaveGuideStructure{1} = ...
                                TwoPortDeviceGetPortSegment(TwoPortDevices{d1(i)}.D,2);
                            CurrentWaveGuideStructure{1}.l = 0;
                        else
                            CurrentWaveGuideStructure = TwoPortDevices{d1(i)}.D;
                        end
                    end
                end
            end
        else
            CurrentWaveGuideStructure = TwoPortDevices{d1(i)}.D;
        end
        CurrentWaveGuideStructure{1}.zo = Position;
        EndPosition = TwoPortDeviceDraw(CurrentWaveGuideStructure, flag);
    end

    Nto1Connections{index}.zo = EndPosition;
    Nto1DeviceDraw(TwoPortDevices, Nto1Connections{index});
    Position = EndPosition;

    for i=1:length(Nto1Connections{index}.SideTwo)
        d2(i) = Nto1Connections{index}.SideTwo{i}.TwoPortDeviceIndex;
        TwoPortDevices{d2(i)}.D{1}.zo = Position;
        EndPosition = TwoPortDeviceDraw(TwoPortDevices{d2(i)}.D, flag);

        if(index < length(Nto1Connections))
            LastSegment = length(TwoPortDevices{d2(i)}.D);
            WaveGuideConnectionCapDraw( TwoPortDevices{d2(i)}.D{LastSegment},...
                'c', EndPosition );
            Position = EndPosition;
        end
        if(index == length(Nto1Connections) && ...
                Options.DeviceSymmetry.Use == 1 && Options.DeviceSymmetry.Side == 2)
            SymmetryConnection = TwoPortDeviceGetPortSegment(TwoPortDevices{d2(i)}.D,2);
            WaveGuideConnectionCapDraw(SymmetryConnection, 'm', EndPosition);
        end
    end
end

axis('equal');
axis('off');
