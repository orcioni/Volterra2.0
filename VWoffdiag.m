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

% function k=VWoff(xn, yn, ord, os, R, A, swap, delay)
% It's a wrapper function to easily call a kernel identification of order
% from 0 to 5.
% The input parameters are passed to the function LeeSchN of order N.
% See LeeSchN functions help for details.

function k=VWoffdiag(xn,yn,ord,os,R,A,swap,delay)

try swap; catch swap=23;end
try delay; catch delay=0;end



switch ord
case 0
    [k]=LeeSch0(yn,os);
case 1
    [k]=LeeSch1(xn,yn,os,R,A,swap,delay);
case 2
    [k]=LeeSch2(xn,yn,os,R,A,swap,delay);
case 3 
    [k]=LeeSch3(xn,yn,os,R,A,swap,delay);
case 4
    [k]=LeeSch4(xn,yn,os,R,A,swap,delay);
case 5
    [k]=LeeSch5(xn,yn,os,R,A,swap,delay);
end
