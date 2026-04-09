function [ StotRN ] = Renormalize( WaveGuideStructure,Stot )
% Renormalize a two-port GSM so every modal sub-block uses the chosen
% waveguide normalization.
% ------------------------------------------------------------------------
% [IN]
%   WaveGuideStructure - 2-element waveguide structure used as reference
%   Stot               - GSM struct to renormalize
%
% [OUT]
%   StotRN             - renormalized GSM struct
%
Np = length(WaveGuideStructure);

Nh1 = WaveGuideStructure{1}.Nh;
Ne1 = WaveGuideStructure{1}.Ne;
Nh2 = WaveGuideStructure{Np}.Nh;
Ne2 = WaveGuideStructure{Np}.Ne;

zeta = 120*pi;

StotRN.S11=zeros(Nh1+Ne1,Nh1+Ne1);
StotRN.S21=zeros(Nh2+Ne2,Nh1+Ne1);
StotRN.S12=zeros(Nh1+Ne1,Nh2+Ne2);
StotRN.S22=zeros(Nh2+Ne2,Nh2+Ne2);

% Renormalizes S11
for i=1:Nh1
    for j=1:Nh1
        StotRN.S11(i,j) = Stot.S11(i,j) * ...
            sqrt(WaveGuideStructure{1}.kh.mn(i)/WaveGuideStructure{1}.kh.mn(j));
    end
end
for i=1:Nh1
    for j=1:Ne1
        StotRN.S11(i,Nh1+j) = Stot.S11(i,Nh1+j) * ...
            sqrt(WaveGuideStructure{1}.kh.mn(i)/WaveGuideStructure{1}.ke.mn(j))*zeta;
    end
end
for i=1:Ne1
    for j=1:Nh1
        StotRN.S11(Nh1+i,j) = Stot.S11(Nh1+i,j) * ...
            sqrt(WaveGuideStructure{1}.ke.mn(i)/WaveGuideStructure{1}.kh.mn(j))/zeta;
    end
end
for i=1:Ne1
    for j=1:Ne1
        StotRN.S11(Nh1+i,Nh1+j) = Stot.S11(Nh1+i,Nh1+j) * ...
            sqrt(WaveGuideStructure{1}.ke.mn(i)/WaveGuideStructure{1}.ke.mn(j));
    end
end

% Renormalizes S21
for i=1:Nh2
    for j=1:Nh1
        StotRN.S21(i,j) = Stot.S21(i,j) * ...
            sqrt(WaveGuideStructure{Np}.kh.mn(i)/WaveGuideStructure{1}.kh.mn(j));
    end
end
for i=1:Nh2
    for j=1:Ne1
        StotRN.S21(i,Nh1+j) = Stot.S21(i,Nh1+j) * ...
            sqrt(WaveGuideStructure{Np}.kh.mn(i)/WaveGuideStructure{1}.ke.mn(j))*zeta;
    end
end
for i=1:Ne2
    for j=1:Nh1
        StotRN.S21(Nh2+i,j) = Stot.S21(Nh2+i,j) * ...
            sqrt(WaveGuideStructure{Np}.ke.mn(i)/WaveGuideStructure{1}.kh.mn(j))/zeta;
    end
end
for i=1:Ne2
    for j=1:Ne1
        StotRN.S21(Nh2+i,Nh1+j) = Stot.S21(Nh2+i,Nh1+j) * ...
            sqrt(WaveGuideStructure{Np}.ke.mn(i)/WaveGuideStructure{1}.ke.mn(j));
    end
end

% Renormalizes S12
for i=1:Nh1
    for j=1:Nh2
        StotRN.S12(i,j) = Stot.S12(i,j) * ...
            sqrt(WaveGuideStructure{1}.kh.mn(i)/WaveGuideStructure{Np}.kh.mn(j));
    end
end
for i=1:Nh1
    for j=1:Ne2
        StotRN.S12(i,Nh1+j) = Stot.S12(i,Nh2+j) * ...
            sqrt(WaveGuideStructure{1}.kh.mn(i)/WaveGuideStructure{Np}.ke.mn(j))*zeta;
    end
end
for i=1:Ne1
    for j=1:Nh2
        StotRN.S12(Nh1+i,j) = Stot.S12(Nh1+i,j) * ...
            sqrt(WaveGuideStructure{1}.ke.mn(i)/WaveGuideStructure{Np}.kh.mn(j))/zeta;
    end
end
for i=1:Ne1
    for j=1:Ne2
        StotRN.S12(Nh1+i,Nh2+j) = Stot.S12(Nh1+i,Nh2+j) * ...
            sqrt(WaveGuideStructure{1}.ke.mn(i)/WaveGuideStructure{Np}.ke.mn(j));
    end
end

% Renormalizes S22
for i=1:Nh2
    for j=1:Nh2
        StotRN.S22(i,j) = Stot.S22(i,j) * ...
            sqrt(WaveGuideStructure{Np}.kh.mn(i)/WaveGuideStructure{Np}.kh.mn(j));
    end
end
for i=1:Nh2
    for j=1:Ne2
        StotRN.S22(i,Nh2+j) = Stot.S22(i,Nh2+j) * ...
            sqrt(WaveGuideStructure{Np}.kh.mn(i)/WaveGuideStructure{Np}.ke.mn(j))*zeta;
    end
end
for i=1:Ne2
    for j=1:Nh2
        StotRN.S22(Nh2+i,j) = Stot.S22(Nh2+i,j) * ...
            sqrt(WaveGuideStructure{Np}.ke.mn(i)/WaveGuideStructure{Np}.kh.mn(j))/zeta;
    end
end
for i=1:Ne2
    for j=1:Ne2
        StotRN.S22(Nh2+i,Nh2+j) = Stot.S22(Nh2+i,Nh2+j) * ...
            sqrt(WaveGuideStructure{Np}.ke.mn(i)/WaveGuideStructure{Np}.ke.mn(j));
    end
end

