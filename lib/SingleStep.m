function [ S ] = SingleStep( WaveGuideStructure, segment, k0 )
% Compute the GSM of one discontinuity between adjacent waveguide
% segments.
% ------------------------------------------------------------------------
% [IN]
%   WaveGuideStructure - cell array of WaveGuideSegment structs with modal
%                        data available at the current frequency
%   segment            - index of the left segment forming the step
%   k0                 - free-space wavenumber [rad/m]
%
% [OUT]
%   S                  - GSM struct for the selected step
%
if (WaveGuideStructure{segment}.a <= WaveGuideStructure{segment+1}.a & ...
    WaveGuideStructure{segment}.b <= WaveGuideStructure{segment+1}.b     )
    % Se questo � vero sto facendo un BOUNDARY ENLARGMENT
    step{1} = WaveGuideStructure{segment};
    step{2} = WaveGuideStructure{segment+1};
    Rflag=0;
else
    % Senn� � un REDUCTION e faccio il problema ENLARGMENT duale...
    step{1} = WaveGuideStructure{segment+1};
    step{2} = WaveGuideStructure{segment};
    Rflag=1;
end

% Matrici M

% Matrice hh
Mhh = [];
for i=1:step{1}.Nh
    for j=1:step{2}.Nh
        I = Integrals(step{1}.kh.x(i), step{1}.kh.y(i), step{2}.kh.x(j), step{2}.kh.y(j), step);
        Mhh(i,j) = step{1}.Ah(i) * step{2}.Ah(j) * ...
            (step{1}.kh.y(i) * step{2}.kh.y(j) * (I(1)-I(2)) * (I(3)+I(4)) + ...
             step{1}.kh.x(i) * step{2}.kh.x(j) * (I(1)+I(2)) * (I(3)-I(4)) );
    end
end

% Matrice he
Mhe = [];
for i=1:step{1}.Nh
    for j=1:step{2}.Ne
        I = Integrals(step{1}.kh.x(i), step{1}.kh.y(i), step{2}.ke.x(j), step{2}.ke.y(j), step);
        Mhe(i,j) = step{1}.Ah(i) * step{2}.Ae(j) * ...
            (-step{1}.kh.y(i) * step{2}.ke.x(j) * (I(1)-I(2)) * (I(3)+I(4)) + ...
             step{1}.kh.x(i) * step{2}.ke.y(j) * (I(1)+I(2)) * (I(3)-I(4)) );
    end
end
if (isempty(Mhe))
    Mhe = zeros(step{1}.Nh,0);
end

% Matrice eh
Meh = [];
for i=1:step{1}.Ne
    for j=1:step{2}.Nh
        Meh(i,j) = 0;
    end
end
if (isempty(Meh))
    Meh = zeros(0,step{2}.Nh);
end

% Matrice ee
Mee = [];
for i=1:step{1}.Ne
    for j=1:step{2}.Ne
        I = Integrals(step{1}.ke.x(i), step{1}.ke.y(i), step{2}.ke.x(j), step{2}.ke.y(j), step);
        Mee(i,j) = step{1}.Ae(i) * step{2}.Ae(j) * ...
            (step{1}.ke.x(i) * step{2}.ke.x(j) * (I(1)-I(2)) * (I(3)+I(4)) + ...
             step{1}.ke.y(i) * step{2}.ke.y(j) * (I(1)+I(2)) * (I(3)-I(4)) );
    end
end

Dgamma_h_1 = 1i*diag(step{1}.kh.mn);
Dgamma_h_2 = 1i*diag(step{2}.kh.mn);
Dgamma_e_1 = 1i*diag(step{1}.ke.mn);
Dgamma_e_2 = 1i*diag(step{2}.ke.mn);
U_h_1 = eye(step{1}.Nh);
U_h_2 = eye(step{2}.Nh);
U_e_1 = eye(step{1}.Ne);
U_e_2 = eye(step{2}.Ne);

lambda = 2*pi / k0;
Z_lambda = 1i*2*pi*120*pi/lambda;
Y_lambda = 1i*2*pi/(120*pi*lambda);

DH = [[Dgamma_h_1, zeros(step{1}.Nh, step{1}.Ne)];...
      [zeros(step{1}.Ne, step{1}.Nh), U_e_1]];
MH = [[Mhh*Dgamma_h_2, Y_lambda*Mhe];...
      [Meh*Dgamma_h_2/Y_lambda, Mee]];

DE = [[U_h_2, zeros(step{2}.Nh, step{2}.Ne)];...
      [zeros(step{2}.Ne, step{2}.Nh),Dgamma_e_2]];
ME = [[Mhh',Meh'*Dgamma_e_1/Z_lambda];...
      [Z_lambda*Mhe',Mee'*Dgamma_e_1]];

H = inv(DH)*MH;
E = inv(DE)*ME;

S11 = inv(eye(step{1}.Nh+step{1}.Ne)+H*E) * (eye(step{1}.Nh+step{1}.Ne)-H*E);
S12 = 2*inv(eye(step{1}.Nh+step{1}.Ne)+H*E)*H;
S21 = E*(eye(step{1}.Nh+step{1}.Ne)+S11);
S22 = E*S12-eye(step{2}.Nh+step{2}.Ne);

if ( Rflag == 0 )
    % Boundary enlargment
    S.S11 = S11;
    S.S12 = S12;
    S.S21 = S21;
    S.S22 = S22;
else
    % Boundary reduction
    S.S11 = S22;
    S.S12 = S21;
    S.S21 = S12;
    S.S22 = S11;
end

