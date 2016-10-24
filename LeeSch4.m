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

%function [k4]=LeeSch4(xn, yn, os, R, A, swap, delay)
%
% k4 is the fourth order kernel of Wiener series according to Lee-Schetzen 
% method, which will contain not-a-NaN values ony within the set of 
% fundamental points,  which is the minimum set of non-diagonal points which 
% allows the reconstruction of non-diagonal kernel points by symmetry 
% (see symmetrize function).
% For diagonal points refer to WVdiag4 function.
%
% xn is the input sequence.
%
% yn the output sequence.
%
% os is the input/output sequences index from where the cross-correlation is 
% started, all the sequence values before os are thrown. In can be used when 
% xn and yn have been obtained from an A/D conversion and we the initial 
% transient conditions cut away.
%
% R is the length corresponding to length(k4)-1. The lags domain interval 
% corresponding to k4 is [0,R]x[0,R]x[0,R]x[0,R].
%
% A is the second order moment of xn (i.e. power).
%
% swap is an optional parameter for speed optimization 
% (see fastxcorr for insights), it depends on the machine you're on, and 
% on having compiled this code or not. 
% If you don't know what to do use default value. 
% (We choose 23 in octave environment and 11 in the standalone version)
%
% delay gives the result restricted to the lags domain interval
% [0+delay,R]x[0+delay,R]x[0+delay,R]x[0+delay,R],
% most useful for higher order kernels.

function [k4]=LeeSch4(xn,yn,os,R,A,swap,delay)

    L=length(xn(os:end-delay));
    A4=A*A*A*A;
    k4=NaNmat(R+1,R+1,R+1,R+1);
    % sgm2=2:R-1
    % sgm1=sgm2+1:R
    % sgm3=1:sgm2-1
    % sgm4=0:sgm3-1   
    for sgm2=2:R-1
        p2=[zeros(sgm2,1);xn(os:end-sgm2-delay)];
        for sgm1=sgm2+1:R;
            p1=[zeros(sgm1,1);xn(os:end-sgm1-delay)];
            for sgm3=1:sgm2-1;
                k4(sgm1+1,sgm2+1,sgm3+1,1:sgm3)=frxcorr(xn(os:end-delay),...
                    [zeros(sgm3,1);xn(os:end-sgm3-delay)].*p2.*p1.*yn(os+delay:end),...
                    L,...
                    sgm3-1,0,swap)/24/A4;
            end
        end
    end
end

