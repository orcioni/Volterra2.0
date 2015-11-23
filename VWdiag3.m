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

% function diag3=VWdiag3(xn, yn, os, R, A, delay, k1)
%
% diag3 is the third order kernel of Wiener series according to our method, 
% which will contain not-a-NaN values ony within the set of fundamental 
% diagonal points, which is the minimum set of diagonal points which allows
% the reconstruction of diagonal kernel points by symmetry 
% (see symmetrize function).
% For non-diagonal points refer to LeeSch3 function.
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
% R is the length corresponding to length(diag3)-1. The lags domain interval
% corresponding to diag3 is [0,R]x[0,R]x[0,R].
%
% A is the second order moment of xn (i.e. power).
%
% delay gives the result restricted to the lags domain interval 
% [0+delay,R]x[0+delay,R]x[0+delay,R],
% most useful for higher order kernels.
%
% k1 is the previously obtained Wiener kernel of the first order.

function diag3=VWdiag3(xn,yn,os,R,A,delay,k1)

diag3=NaNmat(R+1,R+1,R+1);   

A3=A*A*A;

for sgm1=0:R
    p1=[zeros(sgm1,1);xn(os:end-sgm1-delay)];
    for sgm2=sgm1:R
        p2=[zeros(sgm2,1);xn(os:end-sgm2-delay)];
        for sgm3=sgm2:R
            if (sgm3~=sgm1) & (sgm3~=sgm2) & (sgm2~=sgm1)
                break
            else 
                diag3(sgm3+1,sgm2+1,sgm1+1)=1/6/A3* mean(...
                p1.*...
                p2.*...
                [zeros(sgm3,1);xn(os:end-sgm3-delay)].*...
                yn(os+delay:end)    )...
                -1/6/A*(   k1(sgm1+1)*(sgm2==sgm3)+k1(sgm2+1)*(sgm1==sgm3)+k1(sgm3+1)*(sgm1==sgm2)  );
            end
        end
    end
end
     
