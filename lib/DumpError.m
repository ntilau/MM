function [halt] = DumpError( pre, Error )
% Print a formatted Error struct and return whether a fatal condition was
% reported.
% ------------------------------------------------------------------------
% [IN]
%   pre   - prefix string prepended to each printed message
%   Error - Error struct with optional warning, info, or fatal fields
%
% [OUT]
%   halt  - 1 if Error contains a fatal condition, 0 otherwise
%
halt = 0;

if (isfield(Error,'AlreadyDumped'))

    if(isfield(Error,'fatal'))
        halt = 1;
    end

else

    if(isfield(Error,'fatal'))
        disp([pre,Error.fatal]);
        halt = 1;
    end

    if(isfield(Error,'warning'))
        for i=1:length(Error.warning)
            disp([pre,Error.warning{i}]);
        end
    end

    if(isfield(Error,'info'))
        for i=1:length(Error.info)
            disp([pre,Error.info{i}]);
        end
    end
end
