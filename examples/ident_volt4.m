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
sigma_noise=[0.1,0.2,0.4,0.8,1.6];

disp('Identification with multiple variances')
Vkernel_n2 = ident_volt_20(4,memspan,sigma_noise,1e6,'volt4_system');
disp(['Identification with sigma = ' num2str(sigma_noise(1))]);
Vkernel01 = ident_volt(4,memspan,sigma_noise(1),1e6,'volt4_system');
disp(['Identification with sigma = ' num2str(sigma_noise(2))]);
Vkernel02 = ident_volt(4,memspan,sigma_noise(2),1e6,'volt4_system');
disp(['Identification with sigma = ' num2str(sigma_noise(3))]);
Vkernel04 = ident_volt(4,memspan,sigma_noise(3),1e6,'volt4_system');
disp(['Identification with sigma = ' num2str(sigma_noise(4))]);
Vkernel08 = ident_volt(4,memspan,sigma_noise(4),1e6,'volt4_system');
disp(['Identification with sigma = ' num2str(sigma_noise(5))]);
Vkernel16 = ident_volt(4,memspan,sigma_noise(5),1e6,'volt4_system');


xn = randn(1e6,1);

disp(['Test of Volterra system  identified whith sigma = ' num2str(sigma_noise(1))]);
mseyn01 = test_sigma(Vkernel01, 4, 0.1,1.6, 'volt4_system');
disp(['Test of Volterra system  identified whith sigma = ' num2str(sigma_noise(2))]);
mseyn02 = test_sigma(Vkernel02, 4, 0.1,1.6, 'volt4_system');
disp(['Test of Volterra system  identified whith sigma = ' num2str(sigma_noise(3))]);
mseyn04 = test_sigma(Vkernel04, 4, 0.1,1.6, 'volt4_system');
disp(['Test of Volterra system  identified whith sigma = ' num2str(sigma_noise(4))]);
mseyn08 = test_sigma(Vkernel08, 4, 0.1,1.6, 'volt4_system');
disp(['Test of Volterra system  identified whith sigma = ' num2str(sigma_noise(5))]);
mseyn16 = test_sigma(Vkernel16, 4, 0.1,1.6, 'volt4_system');toc;

disp('Test of Volterra system  identified whith multiple variances');
mseyn_n2 = test_sigma(Vkernel_n2, 4, 0.1,1.6, 'volt4_system');

h1  = [0.2264190, 0.8539435, 1.0243269, 0.1957670,  -0.3426567,  -0.0456011, 0.1097026,  -0.0088268,  -0.0177919,   0.0047174];
h2 = h1'*h1;
h3  = zeros (10, 10, 10);
h4  = zeros (10, 10, 10, 10);
h5  = zeros (10, 10, 10, 10, 10);
for i = 1:10
h3(i,:,:)     = h1(i)*h2;
end
for i = 1:10
h4(i,:,:,:)   = h1(i)*h3;
end
for i = 1:10
h5(i,:,:,:,:) = h1(i)*h4;
end
h2 =  9/54*h2;
h3 = -27/486*h3;
h4 = -135/5832*h4;

msenh4(1)=mse(Vkernel02.h4,h4)/mse(h4);
msenh4(2)=mse(Vkernel04.h4,h4)/mse(h4);
msenh4(3)=mse(Vkernel08.h4,h4)/mse(h4);
msenh4(4)=mse(Vkernel16.h4,h4)/mse(h4);

msenh4_n2=mse(Vkernel_n2.h4,h4)/mse(h4);

msenh1(1)=mse(Vkernel02.h1,h1)/mse(h1);
msenh1(2)=mse(Vkernel04.h1,h1)/mse(h1);
msenh1(3)=mse(Vkernel08.h3,h1)/mse(h1);
msenh1(4)=mse(Vkernel16.h3,h1)/mse(h1);

msenh1_n2=mse(Vkernel_n2.h1,h1)/mse(h1);



figure;
loglog(mseyn01(:,1),sqrt(mseyn01(:,2),'-or');hold
xl = xlim;
xl(1) = mseyn02(1,1);
xl(2) = mseyn02(end,1);
xlim(xl);
loglog(mseyn02(:,1),sqrt(mseyn02(:,2)),'-or');
loglog(mseyn04(:,1),sqrt(mseyn04(:,2)),'-or');
loglog(mseyn08(:,1),sqrt(mseyn08(:,2)),'-or');
loglog(mseyn16(:,1),sqrt(mseyn16(:,2)),'-or');
loglog(mseyn_n2(:,1),sqrt(mseyn_n2(:,2)),'-xb');
