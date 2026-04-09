function [ S ] = SingleCascade( SA, D, SB )
% Cascade two step GSMs separated by one uniform-section delay matrix.
% ------------------------------------------------------------------------
% [IN]
%   SA - GSM struct for the first discontinuity
%   D  - delay matrix for the uniform section between the steps
%   SB - GSM struct for the second discontinuity
%
% [OUT]
%   S  - GSM struct for the combined cascade
%
P = inv(eye(size(D)) - D*SB.S11*D*SA.S22);
Q = inv(eye(size(D)) - D*SA.S22*D*SB.S11);

Sg = [[SA.S12,zeros(size(SA.S12,1),size(SB.S21,2))];...
      [zeros(size(SB.S21,1),size(SA.S12,2)),SB.S21]] * ...
    [[P*D*SB.S11,P];[Q,Q*D*SA.S22]] * ...
    [[D*SA.S21,zeros(size(SA.S21,1),size(SB.S12,2))];...
      [zeros(size(SB.S12,1),size(SA.S21,2)),D*SB.S12]] + ...
    [[SA.S11,zeros(size(SA.S11,1),size(SB.S22,2))];
      [zeros(size(SB.S22,1),size(SA.S11,2)),SB.S22]];

[NA,MA] = size (SA.S11);
[NB,MB] = size (SB.S22);

S.S11 =Sg(1:NA,1:MA);
S.S12 =Sg(1:NA,(MA+1):(MA+MB));
S.S21 =Sg((NA+1):(NA+NB),1:MA);
S.S22 =Sg((NA+1):(NA+NB),(MA+1):(MA+MB));

