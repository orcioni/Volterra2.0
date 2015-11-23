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

function  out = volt3_system(x)
h1  = [0.2264190   0.8539435   1.0243269   0.1957670  -0.3426567  -0.0456011 0.1097026  -0.0088268  -0.0177919   0.0047174];
h2n = h1'*h1;
h3  = zeros (10, 10, 10);
h3n = zeros (10, 10, 10);
for i = 1:10
	h3n(i,:,:)     = h1(i)*h2n;
end
h2 =  9/54*h2n;
h3 = -27/486*h3n;
y1 = VoltFilt (x, h1, 1, 0);
y2 = VoltFilt (x, h2, 2, 0);
y3 = VoltFilt (x, h3, 3, 0);
out = y1+y2+y3;
end

