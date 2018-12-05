% function [f2,f3, ...,f10]=FondExtract(h2,h3,...,h10)
%
% The fundamental points constitute the smallest set of kernel points which 
% permits the reconstruction of the kernel by symmetry 
% (i.e. f10=FondExtract(h10); h10=symmetrize(f10), see also symmetrize).
% This function puts NaN values on non-fundamental points of h2 whatever 
% the dimension of the array is (but less than 11), if it is provided as 
% input alone.
% If a list of crescent order arrays/kernels is provided as argument the 
% output will be the corresponding resulting array list with diagonal points
% marked as NaN.
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

function varargout=FondExtract(h2,h3,h4,h5,h6,h7,h8,h9,h10)

single=0;

switch nargin
case 1
    ord=length(size(h2));
    if ord>2,single=ord;
        eval(['h' num2str(ord) '=h2; clear(''h2'')']);
    end
otherwise
    ord=nargin+1;
end

if (ord>=2) && ((single==0)||(single==2))
    R=size(h2,1);
    f2=NaNmat(R,R);
    for tau2=1:R
        for tau1=tau2:R
        f2(tau1,tau2)=h2(tau1,tau2);
        end
    end
    varargout{1}=f2;
    clear f2;
end

if (ord>=3) && ((single==0)||(single==3))
    R=size(h3,1);
    f3=NaNmat(R,R,R);
    for tau3=1:R
        for tau2=tau3:R
            for tau1=tau2:R
                   f3(tau1,tau2,tau3)=h3(tau1,tau2,tau3);
            end
        end
    end
    if single==3, varargout{1}=f3;else varargout{2}=f3; clear f3;end
end

if (ord>=4) && ((single==0)||(single==4))
    R=size(h4,1);
    f4=NaNmat(R,R,R,R);
    eqflag=0;
    for tau4=1:R
        for tau3=tau4:R
            for tau2=tau3:R
                for tau1=tau2:R
                    f4(tau1,tau2,tau3,tau4)=h4(tau1,tau2,tau3,tau4);
                end
            end
        end
    end
    if single==4, varargout{1}=f4;else varargout{3}=f4; clear f4;end
end

if (ord>=5) && ((single==0)||(single==5))
    R=size(h5,1);
    f5=NaNmat(R,R,R,R,R);
    eqflag=0;
    for tau5=1:R
        for tau4=tau5:R
            for tau3=tau4:R
                for tau2=tau3:R
                    for tau1=tau2:R
                        f5(tau1,tau2,tau3,tau4,tau5)=h5(tau1,tau2,tau3,tau4,tau5);
                    end
                end
            end
        end
    end
    if single==5, varargout{1}=f5;else varargout{4}=f5; clear f5;end
end
if (ord>=6) && ((single==0)||(single==6))
end
if (ord>=7) && ((single==0)||(single==7))
end
if (ord>=8) && ((single==0)||(single==8))
end
if (ord>=9) && ((single==0)||(single==9))
end
if (ord>=10) && ((single==0)||(single==10))
end

