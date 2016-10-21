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

% function diag5=VWdiag5(xn, yn, os, R, A, delay, k1, k3)
%
% diag5 is the fifth order kernel of Wiener series according to our method,
% which will contain not-a-NaN values ony within the set of fundamental
% diagonal points, which is the minimum set of diagonal points which allows
% the reconstruction of diagonal kernel points by symmetry 
% (see symmetrize function).
% For non-diagonal points refer to LeeSch5 function.
%
% xn is the input sequence.
% yn the output sequence.
%
% os is the input/output sequences index from where the cross-correlation is
% started, all the sequence values before os are thrown. In can be used when
% xn and yn have been obtained from an A/D conversion and we the initial
% transient conditions cut away.
%
% R is the length corresponding to length(diag5)-1. The lags domain interval
% corresponding to diag5 is [0,R]x[0,R]x[0,R]x[0,R]x[0,R].
%
% A is the second order moment of xn (i.e. power).
%
% delay gives the result restricted to the lags domain interval 
% [0+delay,R]x[0+delay,R]x[0+delay,R]x[0+delay,R]x[0+delay,R], most useful 
% for higher order kernels.
%
% k1 and k3 are previously obtained Wiener kernels of the first and
% third order respectively.

function diag5=VWdiag5(xn,yn,os,R,A,delay,k1,k3)

diag5=NaNmat(R+1,R+1,R+1,R+1,R+1);   

A5=A*A*A*A*A;
A2=A*A;

for sgm1=0:R
    p1=[zeros(sgm1,1);xn(os:end-sgm1-delay)];
    for sgm2=sgm1:R
        p2=[zeros(sgm2,1);xn(os:end-sgm2-delay)];
        for sgm3=sgm2:R
            p3=[zeros(sgm3,1);xn(os:end-sgm3-delay)];
            for sgm4=sgm3:R
                p4=[zeros(sgm4,1);xn(os:end-sgm4-delay)];
                for sgm5=sgm4:R
                    ind=sort([sgm1 sgm2 sgm3 sgm4 sgm5]);
                    if (ind(5)>ind(4)) && (ind(4)>ind(3)) && (ind(3)>ind(2)) && (ind(2)>ind(1)) 
                        break
                    else 
                        Sgm1=sgm1+1;Sgm2=sgm2+1;Sgm3=sgm3+1;Sgm4=sgm4+1;Sgm5=sgm5+1;
                        diag5(Sgm5,Sgm4,Sgm3,Sgm2,Sgm1)=1/120/A5* mean(...
                            p1.*...
                            p2.*...
                            p3.*...
                            p4.*...
                            [zeros(sgm5,1);xn(os:end-sgm5-delay)].*...
                            yn(os+delay:end)                               )...
                            -1/20/A*(...
                                                  k3(Sgm3,Sgm2,Sgm1)*(sgm4==sgm5)...
                                                 +k3(Sgm4,Sgm2,Sgm1)*(sgm3==sgm5)...
                                                 +k3(Sgm5,Sgm2,Sgm1)*(sgm3==sgm4)...
                                                 +k3(Sgm4,Sgm3,Sgm1)*(sgm2==sgm5)...
                                                 +k3(Sgm5,Sgm3,Sgm1)*(sgm2==sgm4)...
                                                 +k3(Sgm5,Sgm4,Sgm1)*(sgm2==sgm3)...
                                                 +k3(Sgm4,Sgm3,Sgm2)*(sgm1==sgm5)...
                                                 +k3(Sgm5,Sgm3,Sgm2)*(sgm1==sgm4)...
                                                 +k3(Sgm5,Sgm4,Sgm2)*(sgm1==sgm3)...
                                                 +k3(Sgm5,Sgm4,Sgm3)*(sgm1==sgm2)...
                                                           )...
                                 -1/120/A2*(...
                            k1(Sgm1)*((sgm2==sgm5)*(sgm3==sgm4)+(sgm3==sgm5)*(sgm2==sgm4)+(sgm4==sgm5)*(sgm2==sgm3))+...
                            k1(Sgm2)*((sgm1==sgm5)*(sgm3==sgm4)+(sgm3==sgm5)*(sgm1==sgm4)+(sgm4==sgm5)*(sgm1==sgm3))+...
                            k1(Sgm3)*((sgm1==sgm5)*(sgm2==sgm4)+(sgm2==sgm5)*(sgm1==sgm4)+(sgm4==sgm5)*(sgm1==sgm2))+...
                            k1(Sgm4)*((sgm1==sgm5)*(sgm2==sgm3)+(sgm2==sgm5)*(sgm1==sgm3)+(sgm3==sgm5)*(sgm1==sgm2))+...
                            k1(Sgm5)*((sgm1==sgm2)*(sgm3==sgm4)+(sgm1==sgm3)*(sgm2==sgm4)+(sgm1==sgm4)*(sgm2==sgm3))...
                                                           );
                    end                       
                end
            end
        end
    end
end   
