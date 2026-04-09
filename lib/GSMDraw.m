function GSMDraw(f, Sf, Sinfo, ModeStruct, flag)
% Plot selected modal S-parameters from a frequency sweep.
% ------------------------------------------------------------------------
% [IN]
%   f          - frequency vector [Hz]
%   Sf         - cell array of GSMs, one per frequency point
%   Sinfo      - port mode-index metadata
%   ModeStruct - cell array describing the modal S-parameters to plot
%   flag       - 1 creates a new figure, 0 plots on current axes
%
% [OUT]
%   (none)     - draws the requested modal traces
%
Color =[' k-',' r-',' g-',' b-',' y-',' m-',' c-'];
Legend = [];
if(flag)
    figure
    for Plot_nbr=1:length(ModeStruct)
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
        hold on;
        plot(f,20*log10(abs(S{Plot_nbr})),char(Color(Plot_nbr*3-2),Color(Plot_nbr*3-1),...
            Color(Plot_nbr*3)));
        ModesLegend = char( ...
            ['S',int2str(outPort), int2str(inPort),'_{', ...
            outMode, int2str(outMode_m), int2str(outMode_n),'_{', ...
            inMode, int2str(inMode_m), int2str(inMode_n),'}}'] ...
            );
        Legend = strvcat(Legend, ModesLegend);
    end
    legend(Legend);
    axis([min(f) max(f) -60 0]);
end
