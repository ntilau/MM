function [ Mhh, Mhe, Meh, Mee] = MxxMatrices(step)
% Assemble the TE/TE, TE/TM, TM/TE, and TM/TM coupling matrices for one
% waveguide step.
% ------------------------------------------------------------------------
% [IN]
%   step - 2-element cell containing the two WaveGuideSegment structs
%
% [OUT]
%   Mhh, Mhe, Meh, Mee - modal overlap matrices used by SingleStep
%
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
