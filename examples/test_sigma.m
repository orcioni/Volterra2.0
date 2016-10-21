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

function mseyn = test_sigma(Vkernel,orderid,sigmastart, sigmaend,test_system)

sigma_index_start=log10(sigmastart/4)/log10(2);
sigma_index_stop=log10(sigmaend*4)/log10(2);
sigmav=2.^[sigma_index_start:(sigma_index_stop-sigma_index_start)/(5*log10(sigmaend*4/sigmastart)/log10(2)):sigma_index_stop];
if is_octave
    randn('seed',1);
else
rng(3);
end
xn = randn(0.5e6,1);
y = zeros(size(xn));
mseyn = zeros(length(sigmav),2);
mseyn(:,1) = sigmav;

for i = 1:length(sigmav)
    y = BasicVoltFilt(sigmav(i)*xn,Vkernel,orderid,0);
    ydes = feval(test_system,sigmav(i)*xn);
    mseyn(i,2) = mse(y,ydes)/mse(ydes);
end
end
