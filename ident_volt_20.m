%function [Vkernel, Wkernel] = ident_volt_20(order,memspan,input_matrix,output_matrix)
%function [Vkernel, Wkernel] = ident_volt_20(order,memspan,input_matrix,des_system)
%function [Vkernel, Wkernel] = ident_volt_20(order,memspan,sigma_noise,dim_input,des_system)
%
% order is the order of the kernel to be identified
%
% memspan is the memory span of the kernels
% if a vector, its values will be assigned to each kernel with this order
% memspan(i) -> kernel_i
%
% input_matrix is the matrix with the vectors to be used in the identification
%
% output_matrix is the matrix with the output vectors to be used in the identification
%
% des_system is a string containing the name of the funtion implementing the system to be
%	indentified
%
% sigmanoise is the vector of the sigma of the noise to be used in the identification
% the noise is generated insede the function instead of being passed as input argument
%
% dim_input is the length of the input vector to be used in the identification
%
% If you want to contact the authors, please write to s.orcioni@univpm.it,
% or Simone Orcioni, DEIT, Università Politecnica delle Marche,
% via Brecce Bianche, 12 - 60131 Ancona, Italy.
% If you are using this program for a scientific work, we encourage you to cite
% the following paper:
%
% Simone Orcioni. Improving the approximation ability of Volterra series 
% identified with a cross-correlation method. Nonlinear Dynamics, 2014.
%
%﻿Orcioni, S., Terenzi, A., Cecchi, S., Piazza, F., & Carini, A. (2018). 
% Identification of Volterra Models of Tube Audio Devices using 
% Multiple-Variance Method. Journal of the Audio Engineering Society, 
% 66(10), 823–838. https://doi.org/10.17743/jaes.2018.0046


% Copyright (C) 2014-2017 Simone Orcioni
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

function [Vkernel, Wkernel] = ident_volt_20(order,memspan,varargin)

if isscalar(memspan)
        memspan = repmat(memspan,order,1);
        delays = zeros(order,1);
elseif isvector(memspan)
        delays = zeros(order,1);
else
    delays = memspan(2,:);
    memspan = memspan(1,:);
end

  switch nargin
    case 5
        %% old style format: internal generated noise
        %%%%%%%%%%%%%%%% disp('Identify kernels')
        sigma_noise = varargin{1};
        dim_input = varargin{2};
        des_system = varargin{3};

        if is_octave
          randn('seed',1);
        else
            %% Added for Replacing Discouraged Syntaxes of rand and randn
            rng(1);
        end
        x0=randn(dim_input,1);
        xn = repmat(x0,order+1);
        xn = xn.*sigma_noise;
        %% VoltFilt works only with vector argument
        %% so we make feval working alsways on vector
        for i = 1:order+1
          yn(:,i) = feval(des_system, xn(:,i));
        end
    case 4     
        %% new format: noise in input function parameters
        if ischar(varargin{2})
            xn = varargin{1};
        %% VoltFilt works only with vector argument
        %% so we make feval working alsways on vector
        for i = 1:order+1
          yn(:,i) = feval(des_system, xn(:,i));
        end
        else
            xn = varargin{1};
            yn = varargin{2};
        end
     otherwise
        error('wrong number of function arguments');
  end

  swaptable = [100,2; 500,4; 1e3,4; 5e3,25; 1e4,31; 5e4,17; 1e5,10; 5e5,10; 1e6,9];
  dim_input = size(xn,1);
  %%%%%%%%%%%%%%%% disp('if we have the swaptable.mat file we use it')
  try
      swap=swapassess(dim_input,'swaptable.mat');
  catch
      swap=swapassess(dim_input,swaptable);
  end
  
  A = var(xn);
  if order >= 0
      Wkernel.k0 = LeeSch0(yn(:,1),1);
      Vkernel.h0 = Wkernel.k0;
  end

  if order >= 1
      k0 = LeeSch0(yn(:,2),1);
      Wkernel.k1 = LeeSch1(xn(:,2),yn(:,2),1,memspan(1)-1,A(2),swap,delays(1));
      [Vkernel.h0, Vkernel.h1] = Wiener2Volterra_20(delays, A, Wkernel.k0, Wkernel.k1);
  end

  if order >= 2
      Wkernel.k2=zeros(memspan(2),memspan(2));
      k0 = LeeSch0(yn(:,3),1);
      k1 = LeeSch1(xn(:,3),yn(:,3),1,memspan(1)-1,A(3),swap,delays(1));
      koff2     = LeeSch2(xn(:,3),yn(:,3),1,memspan(2)-1,A(3),swap,delays(2));
      kdiag2    = VWdiag2(xn(:,3),yn(:,3),1,memspan(2)-1,A(3),delays(2), k0);
      Wkernel.k2 = NaN2zero(symmetrize(koff2)) + NaN2zero(symmetrize(kdiag2));
      [Vkernel.h0, Vkernel.h1, Vkernel.h2] = Wiener2Volterra_20(delays, A, Wkernel.k0, Wkernel.k1, Wkernel.k2);
  end

  if order >= 3
      Wkernel.k3=zeros(memspan(3),memspan(3),memspan(3));

      k0 = LeeSch0(yn(:,4),1);

      k1 = LeeSch1(xn(:,4),yn(:,4),1,memspan(1)-1,A(4),swap,delays(1));

      koff2     = LeeSch2(xn(:,4),yn(:,4),1,memspan(2)-1,A(4),swap,delays(2));
      kdiag2    = VWdiag2(xn(:,4),yn(:,4),1,memspan(2)-1,A(4),delays(2), k0);
      k2 = NaN2zero(symmetrize(koff2)) + NaN2zero(symmetrize(kdiag2));

      koff3      = LeeSch3(xn(:,4),yn(:,4),1,memspan(3)-1,A(4),swap,delays(3));
      kdiag3     = VWdiag3(xn(:,4),yn(:,4),1,memspan(3)-1,A(4),[delays(1) delays(3)],k1);
      Wkernel.k3 = NaN2zero(symmetrize(koff3)) + NaN2zero(symmetrize(kdiag3));

      [Vkernel.h0, Vkernel.h1, Vkernel.h2, Vkernel.h3] = Wiener2Volterra_20(delays, A, Wkernel.k0, Wkernel.k1, Wkernel.k2, Wkernel.k3);
  end

  if order >= 4
      Wkernel.k4=zeros(memspan(4),memspan(4),memspan(4),memspan(4));

      k0 = LeeSch0(yn(:,5),1);

      k1 = LeeSch1(xn(:,5),yn(:,5),1,memspan(1)-1,A(5),swap,delays(1));

      koff2     = LeeSch2(xn(:,5),yn(:,5),1,memspan(2)-1,A(5),swap,delays(2));
      kdiag2    = VWdiag2(xn(:,5),yn(:,5),1,memspan(2)-1,A(5),delays(2), k0);
      k2 = NaN2zero(symmetrize(koff2)) + NaN2zero(symmetrize(kdiag2));

      koff3     = LeeSch3(xn(:,5),yn(:,5),1,memspan(3)-1,A(5),swap,delays(3));
      kdiag3    = VWdiag3(xn(:,5),yn(:,5),1,memspan(3)-1,A(5),[delays(1) delays(3)],k1);
      k3 = NaN2zero(symmetrize(koff3)) + NaN2zero(symmetrize(kdiag3));

      koff4     = LeeSch4(xn(:,5),yn(:,5),1,memspan(4)-1,A(5),swap,delays(4));
      kdiag4    = VWdiag4(xn(:,5),yn(:,5),1,memspan(4)-1,A(5),[delays(2) delays(4)],k0,k2);
      Wkernel.k4= NaN2zero(symmetrize(koff4)) + NaN2zero(symmetrize(kdiag4));

      [Vkernel.h0, Vkernel.h1, Vkernel.h2, Vkernel.h3, Vkernel.h4] = Wiener2Volterra_20(delays, A, Wkernel.k0, Wkernel.k1, Wkernel.k2, Wkernel.k3, Wkernel.k4);
  end

  if order >= 5
      Wkernel.k5=zeros(memspan(5),memspan(5),memspan(5),memspan(5),memspan(5));
      k0 = LeeSch0(yn(:,6),1);

      k1 = LeeSch1(xn(:,6),yn(:,6),1,memspan(1)-1,A(6),swap,delays(1));

      koff2     = LeeSch2(xn(:,6),yn(:,6),1,memspan(2)-1,A(6),swap,delays(2));
      kdiag2    = VWdiag2(xn(:,6),yn(:,6),1,memspan(2)-1,A(6),delays(2), k0);
      k2 = NaN2zero(symmetrize(koff2)) + NaN2zero(symmetrize(kdiag2));

      koff3     = LeeSch3(xn(:,6),yn(:,6),1,memspan(3)-1,A(6),swap,delays(3));
      kdiag3    = VWdiag3(xn(:,6),yn(:,6),1,memspan(3)-1,A(6),[delays(1) delays(3)],k1);
      k3= NaN2zero(symmetrize(koff3)) + NaN2zero(symmetrize(kdiag3));

      koff4     = LeeSch4(xn(:,6),yn(:,6),1,memspan(4)-1,A(6),swap,delays(4));
      kdiag4    = VWdiag4(xn(:,6),yn(:,6),1,memspan(4)-1,A(6),[delays(2) delays(4)],k0,k2);

      koff5     = LeeSch5(xn(:,6),yn(:,6),1,memspan(5)-1,A(6),swap,delays(5));
      kdiag5    = VWdiag5(xn(:,6),yn(:,6),1,memspan(5)-1,A(6),[delays(1) delays(3) delays(5)],k1,k3);
      Wkernel.k5= NaN2zero(symmetrize(koff5)) + NaN2zero(symmetrize(kdiag5));
      [Vkernel.h0, Vkernel.h1, Vkernel.h2, Vkernel.h3, Vkernel.h4,Vkernel.h5] = Wiener2Volterra_20(delays, A, Wkernel.k0, Wkernel.k1, Wkernel.k2, Wkernel.k3, Wkernel.k4, Wkernel.k5);
  end
end
