function RelativePhaseDraw(f, Sf, Sinfo, ModeStruct, flag)
% Plot the relative phase between two selected modal S-parameters.
% ------------------------------------------------------------------------
% [IN]
%   f          - frequency vector [Hz]
%   Sf         - cell array of GSMs
%   Sinfo      - port mode-index metadata
%   ModeStruct - descriptor of the two modal traces to compare
%   flag       - 1 creates a new figure, 0 plots on current axes
%
% [OUT]
%   (none)     - draws the phase difference in degrees
%
Color =['k','r','g','b','y','m','c'];
Legend = [];
if(flag)
    figure
    if(length(ModeStruct)>2)
        warning('Cannot evaluate realtive phase between more than two ports!!');
        return;
    end

    for Plot_nbr=1:2
        Plot_infos = ModeStruct{Plot_nbr};
        outPort = Plot_infos{1};
        inPort = Plot_infos{2};
        outMode = Plot_infos{3};
        outMode_m = Plot_infos{4};
        outMode_n = Plot_infos{5};
        inMode = Plot_infos{6};
        inMode_m = Plot_infos{7};
        inMode_n = Plot_infos{8};
        % md = Plot_infos{9};
        for nf=1:length(f)
            S{Plot_nbr}(nf) = ...
                ExtractSingleS(Sf{nf},Sinfo,outPort,inPort,outMode,...
                outMode_m,outMode_n,inMode,inMode_m,inMode_n);
        end
        ModesLegend = char( ...
            ['S',int2str(outPort), int2str(inPort),'_{', ...
            outMode, int2str(outMode_m), int2str(outMode_n),'_{', ...
            inMode, int2str(inMode_m), int2str(inMode_n),'}}'] ...
            );
        Legend = strvcat(Legend, ModesLegend);
    end
    for nf=1:length(f)
        DeltaAngle = (angle(S{1}(nf))-angle(S{2}(nf)))*180/pi;
        if(abs(DeltaAngle) > 180)
            RelativeAngle(nf) = -360*...
                floor((DeltaAngle/360)+ 0.5)+ DeltaAngle;
        else
            RelativeAngle(nf) = DeltaAngle;
        end
    end
    plot(f, RelativeAngle,'b');
    Legend = sprintf('Relative Phase between %s and %s.', Legend(1,:), Legend(2,:));
    legend(Legend);
    axis([min(f) max(f) -180 180]);
end
