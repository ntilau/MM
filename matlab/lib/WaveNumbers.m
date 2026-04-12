function [ WaveGuideSegment ] = WaveNumbers( WaveGuideSegment, k0 )
% Compute the longitudinal propagation constants of all retained TE and
% TM modes at one frequency.
% ------------------------------------------------------------------------
% [IN]
%   WaveGuideSegment - WaveGuideSegment struct with transverse eigenvalues
%   k0               - free-space wavenumber [rad/m]
%
% [OUT]
%   WaveGuideSegment - input struct augmented with kh.mn and ke.mn
%
for i=1:WaveGuideSegment.Nh
    kh(i) = -1i*sqrt((WaveGuideSegment.mh(i)*pi/WaveGuideSegment.a)^2 + ...
          (WaveGuideSegment.nh(i)*pi/WaveGuideSegment.b)^2 - k0^2);
end
ke = [];
for i=1:WaveGuideSegment.Ne
    ke(i) = -1i*sqrt((WaveGuideSegment.me(i)*pi/WaveGuideSegment.a)^2 + ...
          (WaveGuideSegment.ne(i)*pi/WaveGuideSegment.b)^2 - k0^2);
end

WaveGuideSegment.kh.mn = kh;
WaveGuideSegment.ke.mn = ke;
