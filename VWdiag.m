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

% function kdiag=VWdiag(xn, yn, ord, os, R, delay, varargin)
% Is a wrapper function to easily call a kernel diagonal identification of 
% order from 2 to 5.
% The input parameters are passed to the function WVdiagN of order N.
% See VWdiagN functions help for details.
% varagin should contain the list of lower order prevously identified Wiener
% kernels, e.g.:
% varargin={k0,k2,... k2n} or varargin={k1, .. k(2n+1)} where order 
% N=2,3,4,5 and N=2n or N=2n+1

function kdiag=VWdiag(xn,yn,ord,os,R,A,delay,varargin)

switch ord
case 2
    kdiag=VWdiag2(xn,yn,os,R,A,delay,varargin{1});
case 3 
    kdiag=VWdiag3(xn,yn,os,R,A,delay,varargin{1});
case 4
    kdiag=VWdiag4(xn,yn,os,R,A,delay,varargin{1},varargin{2});
case 5
    kdiag=VWdiag5(xn,yn,os,R,A,delay,varargin{1},varargin{2});
end

    
