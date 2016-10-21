% Copyright (C) 2006 Massimiliano Pirani
% Copyright (C) 2014-2016 Simone Orcioni
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
% the following papers:
% Simone Orcioni. Improving the approximation ability of Volterra series identified
% with a cross-correlation method. Nonlinear Dynamics, 2014.
%
% Simone Orcioni, Massimiliano Pirani, and Claudio Turchetti. Advances in
% Lee-Schetzen method for Volterra filter identification. Multidimensional
% Systems and Signal Processing, 16(3):265-284, 2005.

%% oldstyle function format
%function [Vkernel, Wkernel] = ident_volt(order,memspan,sigma_noise,dim_input,des_system)
% sigmanoise is the vector of the sigma of the noise to be used in the identification
%
% dim_input is the length of the input vector to be used in the identification
%
% des_system is a string containing the name of the funtion implementing the system to be
%	indentified

%% other compatible function formats
%function [Vkernel, Wkernel] = ident_volt(order,memspan,input_vector,des_system)
%function [Vkernel, Wkernel] = ident_volt(order,memspan,input_vector,output_vector)
%
% order is the order of the kernel to be identified
%
% memspan is the memory span of the kernel
%
% input_vector is the input vector to be used in the identification
%
% output_vector is the output vector to be used in the identification
%
% des_system is a string containing the name of the funtion implementing the system to be
%	indentified

function [Vkernel, Wkernel] = ident_volt(order,memspan,varargin)

switch nargin
    case 5
        sigma_noise = varargin{1};
        dim_input   = varargin{2};
        des_system  = varargin{3};
        if is_octave
          randn('seed',1);
        else
        %% Added for Replacing Discouraged Syntaxes of rand and randn
         % version that ensure same results as published
        rng(1,'v4');
        end
        xn = sigma_noise*randn(dim_input,1);
        A = sigma_noise^2;
        yn = feval(des_system,xn);
    case 4
         if ischar(varargin{2})
            xn = varargin{1};
            yn = feval(varargin{2},xn);
         else
            xn = varargin{1};
            yn = varargin{2};
         end
         A = var(xn);
    otherwise
        error('Unexpected number of arguments');
end

swaptable = [100,2; 500,4; 1e3,4; 5e3,25; 1e4,31; 5e4,17; 1e5,10; 5e5,10; 1e6,9];

%%%%%%%%%%%%%%%% disp('if we have the swaptable.mat file we use it')
try
    swap=swapassess(length(xn),'swaptable.mat');
catch
    swap=swapassess(length(xn),swaptable);
end

%%%%%%%%%%%%%%%% disp('Identify kernels')

if order >= 0
	Wkernel.k0 = LeeSch0(yn,1);
end

if order >= 1
	Wkernel.k1 = LeeSch1(xn,yn,1,memspan-1,A,swap,0);
end

if order >= 2
	Wkernel.k2=zeros(memspan,memspan);
	koff2     = LeeSch2(xn,yn,1,memspan-1,A,swap,0);
	kdiag2    = VWdiag2(xn,yn,1,memspan-1,A,0,Wkernel.k0);
	Wkernel.k2= NaN2zero(symmetrize(koff2)) + NaN2zero(symmetrize(kdiag2));
end

if order >= 3
	Wkernel.k3=zeros(memspan,memspan,memspan);
	koff3     = LeeSch3(xn,yn,1,memspan-1,A,swap,0);
	kdiag3    = VWdiag3(xn,yn,1,memspan-1,A,0,Wkernel.k1);
	Wkernel.k3= NaN2zero(symmetrize(koff3)) + NaN2zero(symmetrize(kdiag3));
end

if order >= 4
	Wkernel.k4=zeros(memspan,memspan,memspan,memspan);
	koff4     = LeeSch4(xn,yn,1,memspan-1,A,swap,0);
	kdiag4    = VWdiag4(xn,yn,1,memspan-1,A,0,Wkernel.k0,Wkernel.k2);
	Wkernel.k4= NaN2zero(symmetrize(koff4)) + NaN2zero(symmetrize(kdiag4));
end

if order >= 5
	Wkernel.k5=zeros(memspan,memspan,memspan,memspan,memspan);
	koff5 = VWoffdiag(xn,yn,5,1,length(Wkernel.k5)-1,A,swap,0);
	kdiag5     = VWdiag(xn,yn,5,1,length(Wkernel.k5)-1,A,0,Wkernel.k1,Wkernel.k3);
	Wkernel.k5= NaN2zero(symmetrize(koff5)) + NaN2zero(symmetrize(kdiag5));
end


%%%%%%%%disp('Converting the Wiener system in the equivalent Volterra filter')
%%%%%%%%

switch(order)
case 0
	[Vkernel.h0] = Wiener2Volterra(A,Wkernel.k0);
case 1
	[Vkernel.h0,Vkernel.h1] = Wiener2Volterra(A,Wkernel.k0,Wkernel.k1);
case 2
	[Vkernel.h0,Vkernel.h1,Vkernel.h2] = Wiener2Volterra(A,Wkernel.k0,Wkernel.k1,Wkernel.k2);
case 3
	[Vkernel.h0,Vkernel.h1,Vkernel.h2,Vkernel.h3] = Wiener2Volterra(A,Wkernel.k0,Wkernel.k1,Wkernel.k2,Wkernel.k3);
case 4
	[Vkernel.h0,Vkernel.h1,Vkernel.h2,Vkernel.h3,Vkernel.h4] = Wiener2Volterra(A,Wkernel.k0,Wkernel.k1,Wkernel.k2,Wkernel.k3,Wkernel.k4);
case 5
	[Vkernel.h0,Vkernel.h1,Vkernel.h2,Vkernel.h3,Vkernel.h4,Vkernel.h5] =Wiener2Volterra(A,Wkernel.k0, Wkernel.k1, Wkernel.k2, Wkernel.k3, Wkernel.k4, Wkernel.k5);
otherwise disp('Not a valid order');
endswitch
end
