% Copyright (C) 2014 Simone Orcioni
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
% Simone Orcioni. Improving the approximation ability of Volterra series identified
% with a cross-correlation method. Nonlinear Dynamic, 2014.
%
% Simone Orcioni, Massimiliano Pirani, and Claudio Turchetti. Advances in
% Lee-Schetzen method for Volterra filter identification. Multidimensional
% Systems and Signal Processing, 16(3):265-284, 2005.

memspan=10;
sigma_noise=[0.2,0.4,0.8];

disp('Identification with multiple variances')
Vkernel_n2 = ident_volt_20(2,memspan,[sigma_noise(1) sigma_noise(2) sigma_noise(3)],4e6,'nl_system');
disp(['Identification with sigma = ' num2str(sigma_noise(1))]);
Vkernel03 = ident_volt(2,memspan,sigma_noise(1),4e6,'nl_system');
disp(['Identification with sigma = ' num2str(sigma_noise(2))]);
Vkernel06 = ident_volt(2,memspan,sigma_noise(2),4e6,'nl_system');
disp(['Identification with sigma = ' num2str(sigma_noise(3))]);
Vkernel12 = ident_volt(2,memspan,sigma_noise(3),4e6,'nl_system');

xn = randn(1e6,1);

disp(['Test of Volterra system  identified whith sigma = ' num2str(sigma_noise(1))]);
mseyn03 = test_sigma(Vkernel03, 2, 0.2,1.6, 'nl_system');
disp(['Test of Volterra system  identified whith sigma = ' num2str(sigma_noise(2))]);
mseyn06 = test_sigma(Vkernel06, 2, 0.2,1.6, 'nl_system');
disp(['Test of Volterra system  identified whith sigma = ' num2str(sigma_noise(3))]);
mseyn12 = test_sigma(Vkernel12, 2, 0.2,1.6, 'nl_system');

disp('Test of Volterra system  identified whith multiple variances');
mseyn_n2 = test_sigma(Vkernel_n2, 2, 0.2,1.6, 'nl_system');

figure;
loglog(mseyn03(:,1),sqrt(mseyn03(:,2)),'-or');hold
xl = xlim;
xl(1) = mseyn03(1,1);
xl(2) = mseyn03(end,1);
xlim(xl);
loglog(mseyn06(:,1),sqrt(mseyn06(:,2)),'-or');
loglog(mseyn12(:,1),sqrt(mseyn12(:,2)),'-or');

loglog(mseyn_n2(:,1),mseyn_n2(:,2),'-xb');

%save 'ident_nonlin2_kernel.mat' Vkernel_n2 Vkernel03 Vkernel06 Vkernel12 Vkernel24
%save 'ident_nonlin2_err.mat' mseyn03 mseyn06 mseyn12 mseyn24 mseyn_n2