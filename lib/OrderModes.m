function [ WaveGuideSegment, Error ] = OrderModes( WaveGuideSegment, Symmetry )
% Enumerate and sort TE and TM modes for a waveguide segment according to
% the requested symmetry rules.
% ------------------------------------------------------------------------
% [IN]
%   WaveGuideSegment - struct containing geometry and mode-count limits
%   Symmetry         - symmetry struct controlling admissible modes
%
% [OUT]
%   WaveGuideSegment - input struct augmented with ordered mode indices
%   Error            - Error struct
%
if (WaveGuideSegment.a > WaveGuideSegment.b)
    mstart = 1;
    nstart = 0;
else
    mstart = 0;
    nstart = 1;
end

mstep  = 1;
mstop  = WaveGuideSegment.Nmodes+1;

nstep  = 1;
nstop  = WaveGuideSegment.Nmodes+1;

Error = struct;
infi = 1;

if (Symmetry.x == 1)
    mstep = 2;
    Error.info{infi} = '*** INFO *** 840 - Simmetry along x, only odd m';
    infi = infi + 1;
end

if (Symmetry.y == 1)
    nstep = 2;
    Error.info{infi} = '*** INFO *** 841 - Simmetry along y, only even n';
    infi = infi + 1;
end

if (Symmetry.H == 1)
    mstop = mstart;
    nstop = nstop * 2;
    Error.info{infi} = '*** INFO *** 842 - H plane is uniform - m is constant';
    infi = infi + 1;
end

if (Symmetry.E == 1)
    nstop = nstart;
    mstop = mstop * 2;
    Error.info{infi} = '*** INFO *** 843 - E plane is uniform - n is constant';
    infi = infi + 1;
end

idx = 1;
for m = mstart:mstep:mstop
    for n = nstart:nstep:nstop
        [kx,ky] = OneModeEigens(WaveGuideSegment.a,WaveGuideSegment.b,m,n);
        tmp(idx,1) = (kx)^2 + (ky)^2;
        tmp(idx,2) = m;
        tmp(idx,3) = n;
        idx = idx + 1;
    end
end
tmp = sortrows(tmp,1);

% Count the first Nmodes, but those with both m and n not 0 count double
% since they are ALSO TM (e) modes

ih = 1;
ie = 1;
im = 1;
me=[];
ne=[];
while ( ie+ih-2 < WaveGuideSegment.Nmodes )
    mh(ih) = tmp(im,2);
    nh(ih) = tmp(im,3);
    if (tmp(im,2)>0 & tmp(im,3)>0)
        me(ie) = tmp(im,2);
        ne(ie) = tmp(im,3);
        ie = ie+1;
    end
    ih = ih+1;
    im = im+1;
end

WaveGuideSegment.Nh = length(mh);
WaveGuideSegment.mh = mh;
WaveGuideSegment.nh = nh;

WaveGuideSegment.Ne = length(me);
WaveGuideSegment.me = me;
WaveGuideSegment.ne = ne;

Error.info{infi} = sprintf('*** INFO *** 844 - %d modes requested, %d used ',...
    WaveGuideSegment.Nmodes,WaveGuideSegment.Nh+WaveGuideSegment.Ne);

