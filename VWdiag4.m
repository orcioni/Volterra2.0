% function diag4=VWdiag4(xn, yn, os, R, A, delay, k0, k2)
%
% diag4 is the fourth order kernel of Wiener series according to our method, 
% which will contain not-a-NaN values ony within the set of fundamental 
% diagonal points, which is the minimum set of diagonal points which allows
% the reconstruction of diagonal kernel points by symmetry 
% (see symmetrize function).
% For non-diagonal points refer to LeeSch4 function.
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
% R is the length corresponding to length(diag4)-1. The lags domain interval
% corresponding to diag4 is [0,R]x[0,R]x[0,R]x[0,R].
%
% A is the second order moment of xn (i.e. power).
% delay gives the result restricted to the lags domain interval 
% [0+delay,R]x[0+delay,R]x[0+delay,R]x[0+delay,R], most useful for higher
% order kernels.
%
% k0 and k2 are previously obtained Wiener kernels of the zero-th and
% second order respectively.
%
% If you want to contact the authors, please write to s.orcioni@univpm.it,
% or Simone Orcioni, DII, Università Politecnica delle Marche,
% via Brecce Bianche, 12 - 60131 Ancona, Italy.
% If you are using this program for a scientific work, we encourage you to cite
% the following paper (the file cite.bib, containing the reference in bibtex
% format is also provided):
%
% Simone Orcioni, Massimiliano Pirani, and Claudio Turchetti. Advances in
% Lee-Schetzen method for Volterra filter identification. Multidimensional
% Systems and Signal Processing, 16(3):265-284, 2005.
%
% Simone Orcioni. Improving the approximation ability of Volterra series 
% identified with a cross-correlation method. Nonlinear Dynamics, 2014.
%
%﻿Orcioni, S., Terenzi, A., Cecchi, S., Piazza, F., & Carini, A. (2018). 
% Identification of Volterra Models of Tube Audio Devices using 
% Multiple-Variance Method. Journal of the Audio Engineering Society, 
% 66(10), 823–838. https://doi.org/10.17743/jaes.2018.0046

% Copyright (C) 2006 Massimiliano Pirani
% Copyright (C) 2017 Simone Orcioni
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

function diag4=VWdiag4(xn,yn,os,R,A,delay,k0,k2)

if not(isscalar(delay))
    delay4 = delay(2);
    delay42 = delay(2)-delay(1);
else
    delay4 = delay;
    delay42 = delay;
end    
diag4=NaNmat(R+1,R+1,R+1,R+1);   

A4=A*A*A*A;
A2=A*A;

for sgm1=0:R
    p1=[zeros(sgm1,1);xn(os:end-sgm1-delay4)];
    for sgm2=sgm1:R
        p2=[zeros(sgm2,1);xn(os:end-sgm2-delay4)];
        for sgm3=sgm2:R
            p3=[zeros(sgm3,1);xn(os:end-sgm3-delay4)];
            for sgm4=sgm3:R
                ind=sort([sgm1 sgm2 sgm3 sgm4]);
                if (ind(4)>ind(3)) && (ind(3)>ind(2)) && (ind(2)>ind(1)) 
                    break;
                else 
                    Sgm1=sgm1+1;Sgm2=sgm2+1;Sgm3=sgm3+1;Sgm4=sgm4+1;
                     diag4(Sgm4,Sgm3,Sgm2,Sgm1)=1/24/A4* mean(...
                        p1.*...
                        p2.*...
                        p3.*...
                        [zeros(sgm4,1);xn(os:end-sgm4-delay4)].*...
                        yn(os+delay4:end)                          )...
                        -1/12/A*(...
                        k2(Sgm2+delay42,Sgm1+delay42)*(sgm3==sgm4)...
                       +k2(Sgm3+delay42,Sgm1+delay42)*(sgm2==sgm4)...
                       +k2(Sgm4+delay42,Sgm1+delay42)*(sgm2==sgm3)...
                       +k2(Sgm3+delay42,Sgm2+delay42)*(sgm1==sgm4)...
                       +k2(Sgm4+delay42,Sgm2+delay42)*(sgm1==sgm3)...
                       +k2(Sgm4+delay42,Sgm3+delay42)*(sgm1==sgm2)...
                                 )...
                        -1/24/A2*k0*(...
                        (sgm1==sgm2)*(sgm3==sgm4)+...
                        (sgm1==sgm3)*(sgm2==sgm4)+...
                        (sgm1==sgm4)*(sgm2==sgm3)...
                                    );
                 end   
             end
        end
    end
end
