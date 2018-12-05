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
%
% If you want to contact the authors, please write to s.orcioni@univpm.it,
% or Simone Orcioni, DEIT, Università Politecnica delle Marche,
% via Brecce Bianche, 12 - 60131 Ancona, Italy.
% If you are using this program for a scientific work, we encourage you to cite
% the following papers:
%
% Simone Orcioni. Improving the approximation ability of Volterra series identified
% with a cross-correlation method. Nonlinear Dynamics, 2014.
%
% Simone Orcioni, Massimiliano Pirani, and Claudio Turchetti. Advances in
% Lee-Schetzen method for Volterra filter identification. Multidimensional
% Systems and Signal Processing, 16(3):265-284, 2005.

addpath('./examples');
 memspan=10;
 sigma_noise=[0.2,0.4,0.8];

echo on
fprintf('Identification with multiple variances;')
tic
Vkernel_n2 = ident_volt_20(2,memspan,sigma_noise,1e6,'volt2_system');
toc


fprintf('\n Identification with sigma = 0.2;')
tic
Vkernel02 = ident_volt(2,memspan,sigma_noise(1),1e6,'volt2_system');
toc

fprintf('\n Identification with sigma = 0.4;')
tic
Vkernel04 = ident_volt(2,memspan,sigma_noise(2),1e6,'volt2_system');
toc

fprintf('\n Identification with sigma = 0.8;')
tic
Vkernel08 = ident_volt(2,memspan,sigma_noise(3),1e6,'volt2_system');
toc

xn = randn(1e6,1);

fprintf('\n Test of Volterra system  identified whith sigma = 0.2;')
mseyn02 = test_sigma(Vkernel02, 2, 0.2,1.6, 'volt2_system');


fprintf('\n Test of Volterra system  identified whith sigma = 0.4;')
mseyn04 = test_sigma(Vkernel04, 2, 0.2,1.6, 'volt2_system');


fprintf('\n Test of Volterra system  identified whith sigma = 0.8;')
mseyn08 = test_sigma(Vkernel08, 2, 0.2,1.6, 'volt2_system');


fprintf('\n Test of Volterra system  identified whith multiple variances;')
mseyn_n2 = test_sigma(Vkernel_n2, 2, 0.2,1.6, 'volt2_system');



 h1  = [0.2264190   0.8539435   1.0243269   0.1957670  -0.3426567  -0.0456011 0.1097026  -0.0088268  -0.0177919   0.0047174];
 h2n = h1'*h1;
 h2 =  9/54*h2n;

fprintf('Error on kernels')
 msenh3(1)=mse(Vkernel02.h2,h2)/mse(h2);
 msenh3(2)=mse(Vkernel04.h2,h2)/mse(h2);
 msenh3(3)=mse(Vkernel08.h2,h2)/mse(h2);

 msenh3_n2=mse(Vkernel_n2.h2,h2)/mse(h2);

 msenh1(1)=mse(Vkernel02.h1,h1)/mse(h1);
 msenh1(2)=mse(Vkernel04.h1,h1)/mse(h1);
 msenh1(3)=mse(Vkernel08.h1,h1)/mse(h1);

 msenh1_n2=mse(Vkernel_n2.h1,h1)/mse(h1);



 figure;
 loglog(mseyn02(:,1),sqrt(mseyn02(:,2)),'-or');hold;
 xl = xlim;
 xl(1) = mseyn02(1,1);
 xl(2) = mseyn02(end,1);
 xlim(xl);
 loglog(mseyn04(:,1),sqrt(mseyn04(:,2)),'-or');
 loglog(mseyn08(:,1),sqrt(mseyn08(:,2)),'-or');
 loglog(mseyn_n2(:,1),sqrt(mseyn_n2(:,2)),'-xb');
 xlabel('\sigma_x');
 ax = gca;
 ax.FontSize = 14;
 ylabel('RMSE');