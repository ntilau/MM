function [ Sf, Sinfo, Error, TwoPortDevices, Nto1Connections, ConnectedPorts] = ...
    MultiPortDeviceSolve( TwoPortDevices, Nto1Connections, OpenPorts,...
    ConnectedPorts, FrequencySweep, Symmetry, Topology, Options)
% Solve a MultiPortDevice across the requested frequency sweep.
% ------------------------------------------------------------------------
% [IN]
%   TwoPortDevices, Nto1Connections - validated device definitions
%   OpenPorts, ConnectedPorts       - external and internal port maps
%   FrequencySweep                  - frequency vector definition
%   Symmetry                       - global symmetry struct
%   Topology                       - condensation topology struct
%   Options                        - solver options struct
%
% [OUT]
%   Sf, Sinfo                      - condensed GSM data for each frequency
%   Error                          - Error struct
%   TwoPortDevices, Nto1Connections, ConnectedPorts - updated inputs
%
nbrNto1Connections = length(Nto1Connections);

% Calculate freq independent parameters for all TwoPortDevices (WaveGuideStructure)
for i=1:length(TwoPortDevices)
    for p=1:length(TwoPortDevices{i}.D)
        [ TwoPortDevices{i}.D{p}, Error ] = OrderModes( TwoPortDevices{i}.D{p}, Symmetry );
        [ TwoPortDevices{i}.D{p} ] = EigenModes( TwoPortDevices{i}.D{p} );
        [ TwoPortDevices{i}.D{p} ] = NormCoeff( TwoPortDevices{i}.D{p} );
    end
end

% Evaluate freq sweeped scattering matrices
f = FrequencySweep.f;
c = 299792458;

for nf=1:length(f)

    lambda = c/f(nf);
    k0 = 2*pi/lambda;

    % For each substructure made of nto1connections and the next
    % waveguidestructures.
    PreviousWGSforSinfo = 0;
    GlobalSideOneIndex = 1;
    GlobalSideTwoIndex = 1;

    for index=1:nbrNto1Connections

        nbrSideOne = length(Nto1Connections{index}.SideOne);
        nbrSideTwo = length(Nto1Connections{index}.SideTwo);
        nbrTwoPortDevices = (nbrSideOne + nbrSideTwo);
        toReverse(index) = Topology.Nto1{index}.nFurcation;

        % Controls if WaveGuideStrutcure have already been solved
        SolvedWaveGuideStructure = [];
        LocalPortstoMerge.a = [];
        LocalPortstoMerge.b = [];
        if(index > 1)
            for gamma=1:index-1
                for alpha=1:length(Nto1Connections{gamma}.SideTwo)
                    for beta=1:length(Nto1Connections{index}.SideOne)
                        if (Nto1Connections{gamma}.SideTwo{alpha}.TwoPortDeviceIndex == ...
                                Nto1Connections{index}.SideOne{beta}.TwoPortDeviceIndex)
                            % Found Segment previously solved
                            % Next operations will be to insert a zerolength
                            % segment to prevent duplication of matrices
                            SolvedWaveGuideStructure = horzcat(SolvedWaveGuideStructure, ...
                                Nto1Connections{index}.SideOne{beta}.TwoPortDeviceIndex);
                            % Define local ports to merge related to Nto1
                            % blocks
                            if(Topology.Nto1{gamma}.nFurcation == 0)
                                alpha = alpha + length(Nto1Connections{gamma}.SideOne);
                            end
                            if(Topology.Nto1{index}.nFurcation == 1)
                                beta = beta + length(Nto1Connections{index}.SideTwo);
                            end
                            % Defines
                            LocalPortstoMerge.a = horzcat(LocalPortstoMerge.a, ...
                                [gamma; alpha]);
                            LocalPortstoMerge.b = horzcat(LocalPortstoMerge.b, ...
                                [index; beta]);
                        end
                    end
                end
            end
        end

        % Assigning Segments location for Nto1Junction computation
        d1 = [];
        d2 = [];
        p1 = [];
        p2 = [];
        TwoPortSegmentSideOne = {};
        OpenPortSegmentSideOne = {} ;
        TwoPortSegmentSideTwo = {};
        OpenPortSegmentSideTwo = {} ;

        % To be improved to different modes
        for i=1:nbrSideOne
            d1(i) = Nto1Connections{index}.SideOne{i}.TwoPortDeviceIndex;
            p1(i) = Nto1Connections{index}.SideOne{i}.TwoPortDevicePort;
            [TwoPortSegmentSideOne{i}, OpenPortSegmentSideOne{i}, TwoPortDevices{d1(i)}.D] = ...
                TwoPortDeviceGetPortSegment(TwoPortDevices{d1(i)}.D, p1(i));
            Sinfo{PreviousWGSforSinfo+i}.mh = OpenPortSegmentSideOne{i}.mh;
            Sinfo{PreviousWGSforSinfo+i}.nh = OpenPortSegmentSideOne{i}.nh;
            Sinfo{PreviousWGSforSinfo+i}.me = OpenPortSegmentSideOne{i}.me;
            Sinfo{PreviousWGSforSinfo+i}.ne = OpenPortSegmentSideOne{i}.ne;
        end
        for i=1:nbrSideTwo
            d2(i) = Nto1Connections{index}.SideTwo{i}.TwoPortDeviceIndex;
            p2(i) = Nto1Connections{index}.SideTwo{i}.TwoPortDevicePort;
            [TwoPortSegmentSideTwo{i}, OpenPortSegmentSideTwo{i}, TwoPortDevices{d2(i)}.D]  = ...
                TwoPortDeviceGetPortSegment(TwoPortDevices{d2(i)}.D, p2);
            Sinfo{PreviousWGSforSinfo+nbrSideOne+i}.mh = OpenPortSegmentSideTwo{i}.mh;
            Sinfo{PreviousWGSforSinfo+nbrSideOne+i}.nh = OpenPortSegmentSideTwo{i}.nh;
            Sinfo{PreviousWGSforSinfo+nbrSideOne+i}.me = OpenPortSegmentSideTwo{i}.me;
            Sinfo{PreviousWGSforSinfo+nbrSideOne+i}.ne = OpenPortSegmentSideTwo{i}.ne;
        end

        PreviousWGSforSinfo = PreviousWGSforSinfo + Topology.Nto1{index}.Dimensions;

        % Scattering matrices dimensions (SeleDim) are evaluated. Values must be
        % equal to perform condensation
        for i=1:length(Sinfo)
            SeleDim{i} = [length(Sinfo{i}.mh)+length(Sinfo{i}.me) ,...
                length(Sinfo{i}.nh)+length(Sinfo{i}.ne)];
        end
        for i=1:length(SeleDim)-1
            if(SeleDim{i} == SeleDim{i+1})
            else
                warning('Different number of modes')
                return;
            end
        end

        % Scattering matrices of each WaveGuideStructure are computed and then
        % put into the future GSM
        FutureGSMSideOne = [];
        for i=1:nbrSideOne
            if(length(SolvedWaveGuideStructure) > 0)
                % add a zerolength segment for Nto1 junction analysis
                for beta=1:length(SolvedWaveGuideStructure)
                    if(d1(i) == SolvedWaveGuideStructure(beta))
                        TwoPortSegmentSideOne{i}.l = 0;
                        OpenPortSegmentSideOne{i} = TwoPortSegmentSideOne{i};
                        [S_TwoPortDeviceSideOne{i}, ActualTwoPortDevices{d1(i)}.D] = MultiStep ...
                            (TwoPortSegmentSideOne, k0, Symmetry);
                    end
                end
            else
                if(toReverse(index) == 1)
                    WaveGuideStructure = TwoPortDevices{d1(i)};
                    ReversedWGS = ReverseWaveGuideStructure(WaveGuideStructure.D);
                    [S_TwoPortDeviceSideOne{i}, ActualTwoPortDevices{d1(i)}.D] = ...
                    MultiStep(ReversedWGS, k0, Symmetry);
                else
                    [S_TwoPortDeviceSideOne{i}, TwoPortDevices{d1(i)}.D] = ...
                    MultiStep(TwoPortDevices{d1(i)}.D, k0, Symmetry);
                end
            end
            FutureGSMSideOne = blkdiag(FutureGSMSideOne, S_TwoPortDeviceSideOne{i});
        end

        FutureGSMSideTwo = [];
        for i=1:nbrSideTwo
            if(toReverse(index) == 1)
                WaveGuideStructure = TwoPortDevices{d2(i)};
                ReversedWGS = ReverseWaveGuideStructure(WaveGuideStructure.D);
                [S_TwoPortDeviceSideTwo{i}, ActualTwoPortDevices{d2(i)}.D] = ...
                MultiStep(ReversedWGS, k0, Symmetry);
            else
                [S_TwoPortDeviceSideTwo{i}, TwoPortDevices{d2(i)}.D] = ...
                    MultiStep( TwoPortDevices{d2(i)}.D, k0, Symmetry);
            end
            FutureGSMSideTwo = blkdiag(FutureGSMSideTwo, S_TwoPortDeviceSideTwo{i});
        end

        % GSM without Nto1 matrices
        if(toReverse(index) == 1)
            FutureGSM = blkdiag(FutureGSMSideTwo, FutureGSMSideOne);
        else
            FutureGSM = blkdiag(FutureGSMSideOne, FutureGSMSideTwo);
        end

        % Updates WaveGuideStructures' parameters
        for i=1:nbrSideOne
            % Controls if a fictuous segment has been inserted
            % before updating
            if(length(SolvedWaveGuideStructure) > 0)
                for beta=1:length(SolvedWaveGuideStructure)
                    if(d1(i) ~= SolvedWaveGuideStructure(beta))
                        if(toReverse(index) == 1)
                            TwoPortDevices{d1(i)}.D = ...
                                ReverseWaveGuideStructure(ActualTwoPortDevices{d1(i)}.D);
                        end
                    end
                end
            else
                if(toReverse(index) == 1)
                    TwoPortDevices{d1(i)}.D = ...
                        ReverseWaveGuideStructure(ActualTwoPortDevices{d1(i)}.D);
                end
            end
            [TwoPortSegmentSideOne{i}, OpenPortSegmentSideOne{i}, TwoPortDevices{d1(i)}.D] = ...
                TwoPortDeviceGetPortSegment(TwoPortDevices{d1(i)}.D, p1(i));
        end
        for i=1:nbrSideTwo
            if(toReverse(index) == 1)
                TwoPortDevices{d2(i)}.D = ReverseWaveGuideStructure(ActualTwoPortDevices{d2(i)}.D);
            end
            [TwoPortSegmentSideTwo{i}, OpenPortSegmentSideTwo{i}, TwoPortDevices{d2(i)}.D]  = ...
                TwoPortDeviceGetPortSegment(TwoPortDevices{d2(i)}.D, p2(i));
        end

        % Nto1 junction's scattering matrix calculus
        if(toReverse(index) == 1)
            S_Nto1Connection = Nto1Junction( TwoPortSegmentSideTwo, TwoPortSegmentSideOne{1}, k0);
        else
            S_Nto1Connection = Nto1Junction( TwoPortSegmentSideOne, TwoPortSegmentSideTwo{1}, k0);
        end

        % Finish the GSM
        GSM{index} = blkdiag(FutureGSM, S_Nto1Connection);

        % Condensation routine
        Stmp{index} = CondenseGSM(GSM{index}, SeleDim{1}, Topology.Nto1{index});

        % Saving data for Renormalization
        for i=1:length(d1)
            for j=1:length(OpenPorts)
                if(OpenPorts{j}.TwoPortDeviceIndex == d1(i))
                    GlobalOpenPortSideOne{GlobalSideOneIndex} = OpenPortSegmentSideOne{i};
                    GlobalSideOneIndex = GlobalSideOneIndex + 1;
                end
            end
        end
        for i=1:length(d2)
            for j=1:length(OpenPorts)
                if(OpenPorts{j}.TwoPortDeviceIndex == d2(i))
                    GlobalOpenPortSideTwo{GlobalSideTwoIndex} = OpenPortSegmentSideTwo{i};
                    GlobalSideTwoIndex = GlobalSideTwoIndex + 1;
                end
            end
        end

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Solving the whole structure

    if(nbrNto1Connections > 1)
        newGSM = [];
        newGSMsize = 0;
        for index=1:nbrNto1Connections
            newGSM = blkdiag(newGSM, Stmp{index});
            newGSMsize = newGSMsize + Topology.Nto1{index}.Dimensions;
        end
        tmp = size(LocalPortstoMerge.a);
        nbrPortstoMerge = tmp(2);

        %% Merging Ports
        if(nbrPortstoMerge > 0)
            % Building newTopology for condensation
            newTopology.PortToCondense.a = [];
            newTopology.PortToCondense.b = [];
            newTopology.OpenPorts = [];

            for i=1:nbrPortstoMerge
                Sum = 0;
                SumIndex = 1;
                while SumIndex < LocalPortstoMerge.a(1,i)
                    Sum = Sum + Topology.Nto1{i}.Dimensions;
                    SumIndex = SumIndex + 1;
                end
                newTopology.PortToCondense.a = horzcat(newTopology.PortToCondense.a, ...
                    Sum + LocalPortstoMerge.a(2,i));

                Sum = 0;
                SumIndex = 1;
                while SumIndex < LocalPortstoMerge.b(1,i)
                    Sum = Sum + Topology.Nto1{i}.Dimensions;
                    SumIndex = SumIndex + 1;
                end
                newTopology.PortToCondense.b = horzcat(newTopology.PortToCondense.b, ...
                    Sum + LocalPortstoMerge.b(2,i));
            end
            newGSMremaingPorts = 1:newGSMsize;
            for i=1:newGSMsize
               for j=1:nbrPortstoMerge
                   if(newGSMremaingPorts(i) == newTopology.PortToCondense.a(j) || ...
                           newGSMremaingPorts(i) == newTopology.PortToCondense.b(j))
                       newGSMremaingPorts(i) = 0;
                   end
               end
            end
            % newGSMremainingPortsIndices = find(newGSMremaingPorts);
            PreviousPorts = 0;
            for index=1:nbrNto1Connections
                dim = Topology.Nto1{index}.Dimensions;
                % To use in forward case
                Indices = 1:dim;
                % To use in reverse case
                revIndices(1) = dim;
                for i=1:dim-1
                    revIndices(i+1) = i;
                end
                loop = 1;
                while (loop <= dim)
                    if (toReverse(index) == 1)
                        if (newGSMremaingPorts(PreviousPorts + revIndices(loop)) ~= 0)
                            newTopology.OpenPorts = horzcat(newTopology.OpenPorts,...
                                PreviousPorts + revIndices(loop));
                        end
                    else
                        if (newGSMremaingPorts(PreviousPorts + Indices(loop)) ~= 0)
                            newTopology.OpenPorts = horzcat(newTopology.OpenPorts, ...
                                PreviousPorts + Indices(loop));
                        end
                    end
                    loop = loop + 1;
                end
                PreviousPorts = PreviousPorts + dim;
            end
            Sout = CondenseGSM(newGSM, SeleDim{1}, newTopology);
        elseif(length(ConnectedPorts) > 0)
            %% TO AUTOMATE!!!
            newTopology.PortToCondense.a = [3];
            newTopology.PortToCondense.b = [6];
            newTopology.OpenPorts = [1, 2, 4, 5];
            Sout = CondenseGSM(newGSM, SeleDim{1}, newTopology);
            Sinfo{4} = Sinfo{3};
        end
    else
        Sout = Stmp{1};
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Symmetry Mirroring
    % At the moment done only for the Hildebrand's Coupler
    %% TO AUTOMATE!!!
    if(Options.DeviceSymmetry.Use == 1)
        %Sout = [];

        % GSM of whole problem
        newGSM = blkdiag(Sout, Sout);

        % Condensation routine
        if(Options.DeviceSymmetry.Side == 1)
            newTopology.PortToCondense.a = [1, 2];
            newTopology.PortToCondense.b = [4, 5];
            newTopology.OpenPorts = [3, 6];
            Sout = CondenseGSM(newGSM, SeleDim{1}, newTopology);
            Sinfo{4} = Sinfo{3};
            GlobalOpenPortSideOne = GlobalOpenPortSideTwo;

        elseif(Options.DeviceSymmetry.Side == 2)
            newTopology.PortToCondense.a = [3];
            newTopology.PortToCondense.b = [6];
            newTopology.OpenPorts = [1, 2, 4, 5];
            Sout = CondenseGSM(newGSM, SeleDim{1}, newTopology);
            Sinfo{4} = Sinfo{3};
            GlobalOpenPortSideTwo = GlobalOpenPortSideOne;
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Final Renormalization

    Sout = RenormalizeGSM(Sout, SeleDim{1}, ...
                GlobalOpenPortSideOne, GlobalOpenPortSideTwo);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Sf{nf} = Sout;

end
