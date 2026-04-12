function [Dgamma_h_1, Dgamma_h_2, Dgamma_e_1, Dgamma_e_2] = DgammaMatrices(step)
% Build diagonal propagation-constant matrices for both sides of a step.
% ------------------------------------------------------------------------
% [IN]
%   step - 2-element cell containing the two WaveGuideSegment structs
%
% [OUT]
%   Dgamma_h_1, Dgamma_h_2 - TE diagonal i*kz matrices for sides 1 and 2
%   Dgamma_e_1, Dgamma_e_2 - TM diagonal i*kz matrices for sides 1 and 2
%
Dgamma_h_1 = 1i*diag(step{1}.kh.mn);
Dgamma_h_2 = 1i*diag(step{2}.kh.mn);
Dgamma_e_1 = 1i*diag(step{1}.ke.mn);
Dgamma_e_2 = 1i*diag(step{2}.ke.mn);
