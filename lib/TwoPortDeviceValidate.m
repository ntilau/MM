function [ WaveGuideStructure, Symmetry, Error ] = ...
    TwoPortDeviceValidate( WaveGuideStructure )
% Validate a WaveGuideStructure and add zero-length terminal segments when
% needed for junction congruence.
% ------------------------------------------------------------------------
% [IN]
%   WaveGuideStructure - geometry-only TwoPortDevice definition
%
% [OUT]
%   WaveGuideStructure - validated and possibly expanded structure
%   Symmetry           - inferred symmetry struct
%   Error              - Error struct
%
Np = length(WaveGuideStructure);

iw = 1;

Symmetry.x = 1;
Symmetry.y = 1;
Symmetry.H = 1;
Symmetry.E = 1;

Error = struct();

% Check that all fields in 1st segment are there...
if (isfield(WaveGuideStructure{1},'a'))
    a0  = WaveGuideStructure{1}.a;
else
    Error.fatal = ...
        '[ERROR 0110] Waveguide section 1 has no "a" field';
    return
end

if (isfield(WaveGuideStructure{1},'b'))
    b0  = WaveGuideStructure{1}.b;
else
    Error.fatal = ...
        '[ERROR 0111] Waveguide section 1 has no "b" field';
    return
end

if (isfield(WaveGuideStructure{1},'xo'))
    xo0 = WaveGuideStructure{1}.xo;
else
    Error.fatal = ...
        '[ERROR 0112] Waveguide section 1 has no "xo" field';
    return
end

if (isfield(WaveGuideStructure{1},'yo'))
    yo0 = WaveGuideStructure{1}.yo;
else
    Error.fatal = ...
        '[ERROR 0113] Waveguide section 1 has no "yo" field';
    return
end

if (isfield(WaveGuideStructure{1},'zo'))
    zo0 = WaveGuideStructure{1}.zo;
else
    WaveGuideStructure{1}.zo = 0;
    Error.warning{iw} = ...
        '[WARNING 0114] Waveguide section 1 has no "zo" field, added zo=0 Drawings might be incorrect';
    iw = iw+1;
end

if (isfield(WaveGuideStructure{1},'l'))
    l0  = WaveGuideStructure{1}.l;
else
    Error.fatal = ...
        '[ERROR 0115] Waveguide section 1 has no "l" field';
    return
end

if (~isfield(WaveGuideStructure{1},'Nmodes'))
    Error.fatal = ...
        '[ERROR 0116] Waveguide section 1 has no "Nmodes" field';
    return
end

% La guida deve avere larghezza, altezza e dimensioni POSITIVE
if (a0<eps)
    Error.fatal = ...
        '[ERROR 0102] Waveguide dimension a of section 1 is <=0';
    return
end
if (b0<eps)
    Error.fatal = ...
        '[ERROR 0103] Waveguide dimension b of section 1 is <=0';
    return
end
if (l0<-eps)
    Error.fatal = ...
        '[ERROR 0104] Waveguide length l of section 1 is <0';
    return
end
if (WaveGuideStructure{1}.Nmodes<1)
    Error.fatal = ...
        '[ERROR 0105] Waveguide Nmodes of section 1 is <1';
    return
end

for i=1:Np-1
    % Past section is OK
    a1  = WaveGuideStructure{i}.a;
    b1  = WaveGuideStructure{i}.b;
    xo1 = WaveGuideStructure{i}.xo;
    yo1 = WaveGuideStructure{i}.yo;

    % New section needs check...
    % Check that all fields in 1st segment are there...
    if (isfield(WaveGuideStructure{i+1},'a'))
        a2  = WaveGuideStructure{i+1}.a;
    else
        Error.fatal = ...
            sprintf('[ERROR 0110] Waveguide section %d has no "a" field',i+1);
        return
    end

    if (isfield(WaveGuideStructure{i+1},'b'))
        b2  = WaveGuideStructure{i+1}.b;
    else
        Error.fatal = ...
            sprintf('[ERROR 0111] Waveguide section %d has no "b" field',i+1);
        return
    end

    if (isfield(WaveGuideStructure{i+1},'xo'))
        xo2 = WaveGuideStructure{i+1}.xo;
    else
        Error.fatal = ...
            sprintf('[ERROR 0112] Waveguide section %d has no "xo" field',i+1);
        return
    end

    if (isfield(WaveGuideStructure{i+1},'yo'))
        yo2 = WaveGuideStructure{i+1}.yo;
    else
        Error.fatal = ...
            sprintf('[ERROR 0113] Waveguide section %d has no "yo" field',i+1);
        return
    end

    if (isfield(WaveGuideStructure{i+1},'l'))
        l2  = WaveGuideStructure{i+1}.l;
    else
        Error.fatal = ...
            sprintf('[ERROR 0115] Waveguide section %d has no "l" field',i+1);
        return
    end

    if (~isfield(WaveGuideStructure{i+1},'Nmodes'))
        Error.fatal = ...
            sprintf('[ERROR 0116] Waveguide section %d has no "Nmodes" field',i+1);
        return
    end

    % Verifies symmetries;
    if (abs(xo0-xo2)>eps)
        Symmetry.x = 0;
    end
    if (abs(yo0-yo2)>eps)
        Symmetry.y = 0;
    end
    if (abs(a0-a2)>eps)
        Symmetry.H = 0;
    end
    if (abs(b0-b2)>eps)
        Symmetry.E = 0;
    end

    xmin1 = xo1 - a1/2;
    xmax1 = xo1 + a1/2;
    ymin1 = yo1 - b1/2;
    ymax1 = yo1 + b1/2;

    xmin2 = xo2 - a2/2;
    xmax2 = xo2 + a2/2;
    ymin2 = yo2 - b2/2;
    ymax2 = yo2 + b2/2;

    % La guida deve avere larghezza, altezza e dimensioni POSITIVE
    if (a0<eps)
        Error.fatal = ...
            sprintf('[ERROR 0102] Waveguide dimension a of section %d is <=0',i+1);
        return
    end
    if (b0<eps)
        Error.fatal = ...
            sprintf('[ERROR 0103] Waveguide dimension b of section %d is <=0',i+1);
        return
    end
    if (l0<-eps)
        Error.fatal = ...
            sprintf('[ERROR 0104] Waveguide length l of section %d is <0',i+1);
        return
    end
    if (WaveGuideStructure{i+1}.Nmodes<1)
        Error.fatal = ...
            sprintf('[ERROR 0105] Waveguide Nmodes of section %d is <1',i+1);
        return
    end

    % una guda deve essere TUTTA dentro l'ALTRA
    if (xmin1>=xmin2 && xmax1<=xmax2 && ...
        ymin1>=ymin2 && ymax1<=ymax2)
        % OK, boundary enlargment
        flag = 0;
    elseif (xmin1<=xmin2 && xmax1>=xmax2 && ...
        ymin1<=ymin2 && ymax1>=ymax2)
        % OK, boundary reduction
        flag = 0;
    elseif (xmin1>=xmax2 || xmax1<=xmin2 || ...
        ymin1>=ymax2 || ymax1<=ymin2)
        % Questo � un vero problema, le due guide non si intersecano...
        Error.fatal = ...
            sprintf('[ERROR 0100] Wavegude Segments %d and %d do not intersect',i,i+1);
    else
        % Le guide si intersecano ma non � n�
        % un enlargment n� un reduction, occorre
        % inserire un tratto di lunghezza nulla...
        xmina = max(xmin1,xmin2);
        xmaxa = min(xmax1,xmax2);
        ymina = max(ymin1,ymin2);
        ymaxa = min(ymax1,ymax2);
        NewSegment.a  = xmaxa-xmina;
        NewSegment.b  = ymaxa-ymina;
        NewSegment.xo = (xmaxa+xmina)/2;
        NewSegment.yo = (ymaxa+ymina)/2;
        NewSegment.l  = 0.0;
        NewSegment.Nmodes = min(WaveGuideStructure{i}.Nmodes,WaveGuideStructure{i+1}.Nmodes);
        Error.warning{iw} = ...
            sprintf('[WARNING 0500] Waveguide Segments %d and %d needed an interposed 0 length segment',i,i+1);
        iw = iw + 1;
        WaveGuideStructure = {WaveGuideStructure{1:i},NewSegment,...
            WaveGuideStructure{i+1:Np}};
        Np = Np + 1;
    end
end

if (Symmetry.x==0)
    Symmetry.H=0;
end

if (Symmetry.y==0)
    Symmetry.E=0;
end

if (Symmetry.x==1 && Symmetry.y==1 && Symmetry.H==1 && Symmetry.E==1 ...
        && Np>1)
    Error.fatal = ...
        sprintf('[ERROR 0101] Wavegude Segments are all equal!!');
    return
end

