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
    
% function [h0,h1,h2,...,hp]=Wiener2Volterra(A, k0, k1, k2, k3, k4, k5)
%
% Wiener kernels h0,h1,h2,...,h5 are obtained from Wiener kernels. 
% A is the power of the zero-mean white gaussian input used for Wiener
% identification.
% See Documentation and References provided.

function varargout=Wiener2Volterra(A,k0,k1,k2,k3,k4,k5)

if exist('k1')
    if size(k1,1)==1,k1=k1';end
end


switch nargin
case 2 % zero 
    varargout{1}=k0;
    
case 3 % primo ordine
    varargout{1}=k0;
    varargout{2}=k1;
    
case 4 % secondo ordine
        
    varargout{1}=k0+k02f(k2,A);
    varargout{2}=k1;
    varargout{3}=k2;
    
case 5 % terzo ordine
    
    R1=size(k1,1);
    R3=size(k3,1);
    varargout{1}=k0+k02f(k2,A);
    varargout{2}=[k1;zeros(R3-R1,1)]+[k13f(k3,A);zeros(R1-R3,1)];%k1
    varargout{3}=k2;%h2
    varargout{4}=k3;%h3

case 6 % quarto ordine
    
    R1=size(k1,1);
    R2=size(k2,1);
    R3=size(k3,1);
    R4=size(k4,1);
    
    varargout{1}=k0+k02f(k2,A)+k04f(k4,A);%h0
    varargout{2}=[k1;zeros(R3-R1,1)]+[k13f(k3,A);zeros(R1-R3,1)];%h1
    kz=zeros(max([R4 R2]));
    kz(1:R2,1:R2)=k2;k2=kz;
    kz=zeros(max([R4 R2]));
    k24=k24f(k4,A);
    kz(1:R4,1:R4)=k24;k24=kz;
    varargout{3}=k2+k24;%h2
    varargout{4}=k3;%h3
    varargout{5}=k4;%h4

    
case 7 % quinto ordine
    
    R1=size(k1,1);
    R2=size(k2,1);
    R3=size(k3,1);
    R4=size(k4,1);
    R5=size(k5,1);
    
    R135=max([R1 R3 R5]);
    R24=max([R2 R4]);
    R35=max([R3 R5]);
    
    varargout{1}=k0+k02f(k2,A)+k04f(k4,A);%h0
    varargout{2}=[k1;zeros(R135-R1,1)]+[k13f(k3,A);zeros(R135-R3,1)]+[k15f(k5,A);zeros(R135-R5,1)];%h1
    kz=zeros(R24);
    kz(1:R2,1:R2)=k2;k2=kz;
    kz=zeros(R24);
    k24=k24f(k4,A);
    kz(1:R4,1:R4)=k24;k24=kz;
    varargout{3}=k2+k24;%h2
    kz=zeros(R35,R35,R35);
    kz(1:R3,1:R3,1:R3)=k3;k3=kz;
    k35=k35f(k5,A);
    kz=zeros(R35,R35,R35);
    kz(1:R5,1:R5,1:R5)=k35;k35=kz;
    varargout{4}=k3+k35;%h3
    varargout{5}=k4;%h4
    varargout{6}=k5;%h5
    
otherwise
    error('il numero dei parametri in ingresso non è valido');
end
