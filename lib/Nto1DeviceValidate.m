function [TwoPortDevices, Nto1Connections, Symmetry, Error] = ...
    Nto1DeviceValidate(TwoPortDevices, Nto1Connections)
% Validate an N-to-1 connection and the associated TwoPortDevices.
% ------------------------------------------------------------------------
% [IN]
%   TwoPortDevices  - cell array of TwoPortDevice structs
%   Nto1Connections - Nto1Connection struct to validate
%
% [OUT]
%   TwoPortDevices, Nto1Connections - validated and completed inputs
%   Symmetry, Error                 - inferred symmetry and diagnostics
%
nbrSideOne = length(Nto1Connections.SideOne);
nbrSideTwo = length(Nto1Connections.SideTwo);
Np = nbrSideOne + nbrSideTwo;

iw = 1;

Symmetry.x = 1;
Symmetry.y = 1;
Symmetry.H = 1;
Symmetry.E = 1;

Error = struct();
if (nbrSideOne > nbrSideTwo)
    for i=1:nbrSideOne
        d1(i) = Nto1Connections.SideOne{i}.TwoPortDeviceIndex;
        p1(i) = Nto1Connections.SideOne{i}.TwoPortDevicePort;
        TwoPortSegment{i} = TwoPortDeviceGetPortSegment(TwoPortDevices{d1(i)}.D, p1(i));
    end
    for i=1:nbrSideTwo
        d2(i) = Nto1Connections.SideTwo{i}.TwoPortDeviceIndex;
        p2(i) = Nto1Connections.SideTwo{i}.TwoPortDevicePort;
        TwoPortSegment{Np} = TwoPortDeviceGetPortSegment(TwoPortDevices{d2(i)}.D, p2(i));
    end
else
    for i=1:nbrSideTwo
        d1(i) = Nto1Connections.SideTwo{i}.TwoPortDeviceIndex;
        p1(i) = Nto1Connections.SideTwo{i}.TwoPortDevicePort;
        TwoPortSegment{i} = TwoPortDeviceGetPortSegment(TwoPortDevices{d1(i)}.D, p1(i));
    end
    for i=1:nbrSideOne
        d2(i) = Nto1Connections.SideOne{i}.TwoPortDeviceIndex;
        p2(i) = Nto1Connections.SideOne{i}.TwoPortDevicePort;
        TwoPortSegment{Np} = TwoPortDeviceGetPortSegment(TwoPortDevices{d2(i)}.D, p2(i));
    end
end
% Check that all fields in SideTwo segment are there...

if (isfield(TwoPortSegment{Np},'a'))
    a0  = TwoPortSegment{Np}.a;
else
    Error.fatal = ...
        '[ERROR 0110] Waveguide section 1 has no "a" field';
    return
end

if (isfield(TwoPortSegment{Np},'b'))
    b0  = TwoPortSegment{Np}.b;
else
    Error.fatal = ...
        '[ERROR 0111] Waveguide section 1 has no "b" field';
    return
end

if (isfield(TwoPortSegment{Np},'xo'))
    xo0 = TwoPortSegment{Np}.xo;
else
    Error.fatal = ...
        '[ERROR 0112] Waveguide section 1 has no "xo" field';
    return
end

if (isfield(TwoPortSegment{Np},'yo'))
    yo0 = TwoPortSegment{Np}.yo;
else
    Error.fatal = ...
        '[ERROR 0113] Waveguide section 1 has no "yo" field';
    return
end

if (isfield(TwoPortSegment{Np},'zo'))
    zo0 = TwoPortSegment{Np}.zo;
else
    TwoPortSegment{Np}.zo = 0;
    Error.warning{iw} = ...
        '[WARNING 0114] Waveguide section 1 has no "zo" field, added zo=0 Drawings might be incorrect';
    iw = iw+1;
end

if (isfield(TwoPortSegment{Np},'l'))
    l0  = TwoPortSegment{Np}.l;
else
    Error.fatal = ...
        '[ERROR 0115] Waveguide section 1 has no "l" field';
    return
end

if (~isfield(TwoPortSegment{Np},'Nmodes'))
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
if (TwoPortSegment{Np}.Nmodes<1)
    Error.fatal = ...
        '[ERROR 0105] Waveguide Nmodes of section 1 is <1';
    return
end

for i=1:Np-1
    % Past section is OK
    a1  = TwoPortSegment{Np}.a;
    b1  = TwoPortSegment{Np}.b;
    xo1 = TwoPortSegment{Np}.xo;
    yo1 = TwoPortSegment{Np}.yo;

    % New section needs check...
    % Check that all fields in 1st segment are there...
    if (isfield(TwoPortSegment{i},'a'))
        a2  = TwoPortSegment{i}.a;
    else
        Error.fatal = ...
            sprintf('[ERROR 0110] Waveguide section %d has no "a" field',i+1);
        return
    end

    if (isfield(TwoPortSegment{i},'b'))
        b2  = TwoPortSegment{i}.b;
    else
        Error.fatal = ...
            sprintf('[ERROR 0111] Waveguide section %d has no "b" field',i+1);
        return
    end

    if (isfield(TwoPortSegment{i},'xo'))
        xo2 = TwoPortSegment{i}.xo;
    else
        Error.fatal = ...
            sprintf('[ERROR 0112] Waveguide section %d has no "xo" field',i+1);
        return
    end

    if (isfield(TwoPortSegment{i},'yo'))
        yo2 = TwoPortSegment{i}.yo;
    else
        Error.fatal = ...
            sprintf('[ERROR 0113] Waveguide section %d has no "yo" field',i+1);
        return
    end

    if (isfield(TwoPortSegment{i},'l'))
        l2  = TwoPortSegment{i}.l;
    else
        Error.fatal = ...
            sprintf('[ERROR 0115] Waveguide section %d has no "l" field',i+1);
        return
    end

    if (~isfield(TwoPortSegment{i},'Nmodes'))
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
    if (TwoPortSegment{i}.Nmodes<1)
        Error.fatal = ...
            sprintf('[ERROR 0105] Waveguide Nmodes of section %d is <1',i+1);
        return
    end

    % una guda deve essere TUTTA dentro l'ALTRA
    if (xmin1>=xmin2 && xmax1<=xmax2 && ...
        ymin1>=ymin2 && ymax1<=ymax2)
        % OK, boundary reduction
        flag = 0;
    elseif (xmin1<=xmin2 && xmax1>=xmax2 && ...
        ymin1<=ymin2 && ymax1>=ymax2)
        % OK, boundary enlargement --> should be THIS one
        flag = 0;
    elseif (xmin1>=xmax2 || xmax1<=xmin2 || ...
        ymin1>=ymax2 || ymax1<=ymin2)
        % Questo � un vero problema, le due guide non si intersecano...
        Error.fatal = ...
            sprintf('[ERROR 0100] Waveguide Segments %d and %d do not intersect',i,i+1);
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
        NewSegment.Nmodes = min(TwoPortSegment{i}.Nmodes,TwoPortSegment{Np}.Nmodes);
        Error.warning{iw} = ...
            sprintf('[WARNING 0500] Waveguide Segments %d and %d needed an interposed 0 length segment',i,i+1);
        iw = iw + 1;
        % The Segment has to be inserted into the the relative
        % TwoPortDevices(WaveguideStructure) after or before the correct D
        % (WaveGuideSegment)
         TwoPortDevices{d1(i)} = TwoPortDeviceInsertPortSegment(TwoPortDevices{d1(i)}, NewSegment, p1(i));
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
        sprintf('[ERROR 0101] Waveguide Segments are all equal!!');
    return
end
