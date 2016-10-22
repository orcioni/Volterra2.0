% Copyright (C) 2006 Massimiliano Pirani
%
%  This program is free software; you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published by
%  the Free Software Foundation; either version 2 of the License, or
%  (at your option) any later version.
%
%  This program is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
%
%  You should have received a copy of the GNU General Public License along
%  with this program; if not, write to the Free Software Foundation, Inc.,
%  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
%
% If you want to contact the authors, please write to s.orcioni@univpm.it,
% or Simone Orcioni, DEIT, Università Politecnica delle Marche,
% via Brecce Bianche, 12 - 60131 Ancona, Italy.
% If you are using this program for a scientific work, we encourage you to cite
% the following paper (the file cite.bib, containing the reference in bibtex
% format is also provided):
% Simone Orcioni, Massimiliano Pirani, and Claudio Turchetti. Advances in 
% Lee-Schetzen method for Volterra filter identification. Multidimensional 
% Systems and Signal Processing, 16(3):265-284, 2005.



%function [R, lags] = frxcorr (X, Y, N, maxlag, method, swap)
%
% frxcorr (fast-real xcorr) it's an alternate function to Matlab's and 
% Octave's xcorr, restricted to real column vectors of equal length case 
% and for causal part of cross-correlation system identification 
% (in Wiener sense), created to be Matlab and Octave compatible.
% The main motivation for frxcorr is speed, for real vectors correlation case,
% traded for the generality and versatility of xcorr. 
% The attempt made within this modification was to trade the general purpose 
% characteristics of the original xcorr for a more efficient memory use and 
% speed: note that, within this context, the Wiener identification process may
% involve very long input/output (real valued) sequences and a huge number of 
% calls to this frxcorr function, so the computing time saving using this
% function can be considerable with respecct to the use of xcorr.
% Error check is at a minimum (in fact no parameter check is made at all), 
% so please take care to following strictly the calling convention instructions
% we advice hereafter.
% Usage:
%
% X and Y, have to be the column vectors (with length N > 1) corresponding
% to the vectors being correlated;
%
% N, have to be the length of X and Y;
%
% maxlag, it's the length of the resulting correlation sequence (i.e. R vector) 
% minus 1 and the R vector will represent the rightmost side of the 
% crosscorrelation sequence (the only part needed in case of causality 
% assumption) with x-axis interval being [0,maxlag] (i.e. the lags output 
% vector).
%
% method, has to be one of the following integers:
%   0, corresponding to 'biased'   for correlation=raw/N, 
%   1, corresponding to 'unbiased' for correlation=raw/(N-|lag|), 
%   2, corresponding to 'coeff'    for correlation=raw/(correlation value 
%   at lag 0),
%   >=3, <0, corresponding to 'none'  for raw correlation
%
% swap, it's a positive integer (use 10 if you don't know what to do) 
% which swaps between two versions of the cross-correlation algorhytm, 
% an FFT version if X and Y length it's greater than or equal to swap value 
% or a plain non FFT version (brute force) for sequences with length less than 
% swap. The swap parameter has to be assessed eurhystically, and the efficiency
% improvement gained with its use should be assessed on a cpu architecture 
% (or whole system) basis, also depending on the choice of using a compiled 
% version of this function or not. You can use the swaptest utility
% as reference script in the demo directory of WVItools.
%
% NOTE: Cxy=rxcorr(x,y...) corresponds to the estimation of causal part (m>=0)
% of E[x(n)y*(n+m)]=E[x(n-m)y*(n)], which corresponds to estimation of an 
% impulse response h(m) of a causal linear system with input x and output y
%

function [R, lags] = frxcorr (X, Y, N, maxlag, method, swap)

if maxlag >=swap
    M = 2^nextpow2(N + maxlag);
    y = fft ([Y; zeros(M-N,1)]);
    x = fft ([zeros(maxlag,1);X;zeros(M-N-maxlag,1)]);
    cor = ifft (x.*conj(y));
    R = cor(1:2*maxlag+1);
    R = R(maxlag+1:-1:1);
    R=real(R);
    switch method
        case 0
            %biased
            R = R / N;
        case 1
            %unbiased
            R = R ./ ([N:-1:N-maxlag]');
        case 2
            %coeff
            R = R ./ ( ones(N,1) * R(maxlag+1) );
    end
else
    R=zeros(maxlag+1,1);
    switch method
        case 0 % biased result, i.e. divide by N for each element
            for i=0:maxlag
                R(i+1)=mean([zeros(i,1);X(1:end-i)].*Y);
            end
        case 1 % unbiased result, i.e. divide by N-abs(lag)
            for i=0:maxlag
                R(i+1)=sum([zeros(i,1);X(1:end-i)].*Y)/(N-i);
            end
        case 2 %coeff
            for i=0:maxlag
                R(i+1)=sum([zeros(i,1);X(1:end-i)].*Y);
            end
            R = R/R(1);
        otherwise
            for i=0:maxlag
                R(i+1)=sum([zeros(i,1);X(1:end-i)].*Y);
            end
    end
end
lags = 0:maxlag;
