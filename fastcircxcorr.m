%function [R, lags] = fastcircxcorr (X, Y, N, maxlag, method, swap)
%
% If you want to contact the authors, please write to s.orcioni@univpm.it,
% or Simone Orcioni, DII, Università Politecnica delle Marche,
% via Brecce Bianche, 12 - 60131 Ancona, Italy.
% If you are using this program for a scientific work, we encourage you to cite
% the following paper (the file cite.bib, containing the reference in bibtex
% format is also provided):
%
% Simone Orcioni. Improving the approximation ability of Volterra series 
% identified with a cross-correlation method. Nonlinear Dynamics, 2014.
%
%﻿Orcioni, S., Terenzi, A., Cecchi, S., Piazza, F., & Carini, A. (2018). 
% Identification of Volterra Models of Tube Audio Devices using 
% Multiple-Variance Method. Journal of the Audio Engineering Society, 
% 66(10), 823–838. https://doi.org/10.17743/jaes.2018.0046


% Copyright (C) 2006 Massimiliano Pirani
% Copyright (C) 2018 Simone Orcioni
% Copyright (C) 2018 Alberto Carini
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


function [R, lags] = fastcircxcorr (X, Y, N, maxlag, method, swap)

R=zeros(maxlag+1,1);
switch method
    case 0 % biased result, i.e. divide by N for each element
        for i=0:maxlag
            R(i+1) = mean(Y.*circshift(X,i));
        end
    case 1 % unbiased result, i.e. divide by N-abs(lag)
        for i=0:maxlag
            R(i+1) = sum(Y.*circshift(X,i))/(N-i);
        end
    case 2 %coeff
        for i=0:maxlag
            R(i+1) = sum(Y.*circshift(X,i));
        end
        R = R/R(1);
    otherwise
        for i=0:maxlag
            R(i+1) = sum(Y.*circshift(X,i));
        end
end
lags = 0:maxlag;