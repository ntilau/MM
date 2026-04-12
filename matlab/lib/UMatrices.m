function [U_h_1, U_h_2, U_e_1, U_e_2] = UMatrices(step)
% Build identity matrices of the appropriate TE and TM modal sizes for a
% waveguide step.
% ------------------------------------------------------------------------
% [IN]
%   step - 2-element cell containing the two WaveGuideSegment structs
%
% [OUT]
%   U_h_1, U_h_2 - TE identity matrices for sides 1 and 2
%   U_e_1, U_e_2 - TM identity matrices for sides 1 and 2
%
U_h_1 = eye(step{1}.Nh);
U_h_2 = eye(step{2}.Nh);
U_e_1 = eye(step{1}.Ne);
U_e_2 = eye(step{2}.Ne);
