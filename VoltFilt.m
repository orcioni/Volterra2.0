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


% function out=VoltFilt(xn, kernel, ord, delay)
%
% VoltFilt is an implementation of a Volterra kernel of order ord 
% (ord=10 at most), given the column input vector (xn) and a Volterra 
% kernel array.
%
% delay (default value 0) can be used to account for a possibly causal delay
% (formerly known or assessed) featured by the system we want to model with
% the Volterra filter. 
% If such a delay happens to be, the value of delay parameter avoids the
% explicit computation of the zero valued points of the kernels due to 
% causality delay, indeed often reducing considerably the effective 
% computations burden.
% delay can be also used as a means to trim the lower (in time lag sense) 
% part of kernels.
%
% out will be the total order Volterra filter output to the xn input 
% stimulus. 
% See also BasicVoltFilt (suggested)

function out=VoltFilt(xn,kernel,ord,delay)

R=length(kernel);
R_1=R-1;

L=length(xn);   
if nargin<5
    if nargin==4
        delay=0;
    else
        error('check number and type of function arguments')
    end
end

out=cell(ord+2,1);

if ord==0,out=kernel;end

if ord==1, out=filter([zeros(delay,1);kernel],1,xn);end

if ord==2
    temp=zeros(L,1);
    
    for J1=delay:R_1
        p1=[zeros(J1,1);xn(1:end-J1)];
        for J2=J1:R_1
            if J2==J1
                a=kernel(J1+1,J2+1);
            else
                a=2*kernel(J1+1,J2+1);
            end
            temp=temp+a*p1.*[zeros(J2,1);xn(1:end-J2)];
        end
        out=temp;
    end
end


if ord==3
    temp=zeros(L,1);
    
    for J1=delay:R_1
        p1=[zeros(J1,1);xn(1:end-J1)];
        for J2=J1:R_1
            p2=[zeros(J2,1);xn(1:end-J2)];
            for J3=J2:R_1
                temp=temp+consts(3,6,[J1 J2 J3])*kernel(J1+1,J2+1,J3+1)...
                    *p1...
                    .*p2...
                    .*[zeros(J3,1);xn(1:end-J3)];
            end
        end
        out=temp;
    end
end

if ord==4
    temp=zeros(L,1);
    
    for J1=delay:R_1
        p1=[zeros(J1,1);xn(1:end-J1)];
        for J2=J1:R_1
            p2=[zeros(J2,1);xn(1:end-J2)];
            for J3=J2:R_1
                p3=[zeros(J3,1);xn(1:end-J3)];
                for J4=J3:R_1
                    temp=temp+consts(4,24,[J1 J2 J3 J4])*kernel(J1+1,J2+1,J3+1,J4+1)...
                        *p1...
                        .*p2...
                        .*p3...
                        .*[zeros(J4,1);xn(1:end-J4)];
                end
            end
        end
    end
    out=temp;
end

if ord==5
    temp=zeros(L,1);
    
    for J1=delay:R_1
        p1=[zeros(J1,1);xn(1:end-J1)];
        for J2=J1:R_1
            p2=[zeros(J2,1);xn(1:end-J2)];
            for J3=J2:R_1
                p3=[zeros(J3,1);xn(1:end-J3)];
                for J4=J3:R_1
                    p4=[zeros(J4,1);xn(1:end-J4)];
                    for J5=J4:R_1
                        temp=temp+consts(5,120,[J1 J2 J3 J4 J5])*kernel(J1+1,J2+1,J3+1,J4+1,J5+1)...
                            *p1...
                            .*p2...
                            .*p3...
                            .*p4...
                            .*[zeros(J5,1);xn(1:end-J5)];
                    end
                end
            end
        end
    end
    out=temp;
end    

if ord==6
    temp=zeros(L,1);
    
    for J1=delay:R_1
        p1=[zeros(J1,1);xn(1:end-J1)];
        for J2=J1:R_1
            p2=[zeros(J2,1);xn(1:end-J2)];
            for J3=J2:R_1
                p3=[zeros(J3,1);xn(1:end-J3)];
                for J4=J3:R_1
                    p4=[zeros(J4,1);xn(1:end-J4)];
                    for J5=J4:R_1
                        [zeros(J5,1);xn(1:end-J5)];
                        for J6=J5:R_1
                            temp=temp+consts(6,720,[J1 J2 J3 J4 J5 J6])*kernel(J1+1,J2+1,J3+1,J4+1,J5+1,J6+1)...
                                *p1...
                                .*p2...
                                .*p3...
                                .*p4...
                                .*p5...
                                .*[zeros(J6,1);xn(1:end-J6)];
                        end
                    end
                end
            end
        end
    end
    out=temp;
end    

if ord==7
    temp=zeros(L,1);
    
    for J1=delay:R_1
        p1=[zeros(J1,1);xn(1:end-J1)];
        for J2=J1:R_1
            p2=[zeros(J2,1);xn(1:end-J2)];
            for J3=J2:R_1
                p3=[zeros(J3,1);xn(1:end-J3)];
                for J4=J3:R_1
                    p4=[zeros(J4,1);xn(1:end-J4)];
                    for J5=J4:R_1
                        p5=[zeros(J5,1);xn(1:end-J5)];
                        for J6=J5:R_1
                            p6=[zeros(J6,1);xn(1:end-J6)];
                            for J7=J6:R_1
                                temp=temp+consts(7,5040,[J1 J2 J3 J4 J5 J6 J7])*kernel(J1+1,J2+1,J3+1,J4+1,J5+1,J6+1,J7+1)...
                                    *p1...
                                    .*p2...
                                    .*p3...
                                    .*p4...
                                    .*p5...
                                    .*p6...
                                    .*[zeros(J7,1);xn(1:end-J7)];
                            end
                        end
                    end
                end
            end
        end
    end
    out=temp;
end 

if ord==8
    temp=zeros(L,1);
    
    for J1=delay:R_1
        p1=[zeros(J1,1);xn(1:end-J1)];
        for J2=J1:R_1
            p2=[zeros(J2,1);xn(1:end-J2)];
            for J3=J2:R_1
                p3=[zeros(J3,1);xn(1:end-J3)];
                for J4=J3:R_1
                    p4=[zeros(J4,1);xn(1:end-J4)];
                    for J5=J4:R_1
                        p5=[zeros(J5,1);xn(1:end-J5)];
                        for J6=J5:R_1
                            p6=[zeros(J6,1);xn(1:end-J6)];
                            for J7=J6:R_1
                                p7=[zeros(J7,1);xn(1:end-J7)];
                                for J8=J7:R_1
                                    temp=temp+consts(8,40320,[J1 J2 J3 J4 J5 J6 J7 J8])*kernel(J1+1,J2+1,J3+1,J4+1,J5+1,J6+1,J7+1,J8+1)...
                                        *p1...
                                        .*p2...
                                        .*p3...
                                        .*p4...
                                        .*p5...
                                        .*p6...
                                        .*p7...
                                        .*[zeros(J8,1);xn(1:end-J8)];
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    out=temp;
    
end 

if ord==9
    temp=zeros(L,1);
    
    for J1=delay:R_1
        p1=[zeros(J1,1);xn(1:end-J1)];
        for J2=J1:R_1
            p2=[zeros(J2,1);xn(1:end-J2)];
            for J3=J2:R_1
                p3=[zeros(J3,1);xn(1:end-J3)];
                for J4=J3:R_1
                    p4=[zeros(J4,1);xn(1:end-J4)];
                    for J5=J4:R_1
                        p5=[zeros(J5,1);xn(1:end-J5)];
                        for J6=J5:R_1
                            p6=[zeros(J6,1);xn(1:end-J6)];
                            for J7=J6:R_1
                                p7=[zeros(J7,1);xn(1:end-J7)];
                                for J8=J7:R_1
                                    p8=[zeros(J8,1);xn(1:end-J8)];
                                    for J9=J8:R_1
                                        temp=temp+consts(9,362880,[J1 J2 J3 J4 J5 J6 J7 J8 J9])*kernel(J1+1,J2+1,J3+1,J4+1,J5+1,J6+1,J7+1,J8+1,J9+1)...
                                            *p1...
                                            .*p2...
                                            .*p3...
                                            .*p4...
                                            .*p5...
                                            .*p6...
                                            .*p7...
                                            .*p8...
                                            .*[zeros(J9,1);xn(1:end-J9)];
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    out=temp;
    
end 

if ord==10
    temp=zeros(L,1);
    
    for J1=delay:R_1
        p1=[zeros(J1,1);xn(1:end-J1)];
        for J2=J1:R_1
            p2=[zeros(J2,1);xn(1:end-J2)];
            for J3=J2:R_1
                p3=[zeros(J3,1);xn(1:end-J3)];
                for J4=J3:R_1
                    p4=[zeros(J4,1);xn(1:end-J4)];
                    for J5=J4:R_1
                        p5=[zeros(J5,1);xn(1:end-J5)];
                        for J6=J5:R_1
                            p6=[zeros(J6,1);xn(1:end-J6)];
                            for J7=J6:R_1
                                p7=[zeros(J7,1);xn(1:end-J7)];
                                for J8=J7:R_1
                                    p8=[zeros(J8,1);xn(1:end-J8)];
                                    for J9=J8:R_1
                                        p9=[zeros(J9,1);xn(1:end-J9)];
                                        for J10=J9:R_1
                                            temp=temp+consts(10,3628800,[J1,J2,J3,J4,J5,J6,J7,J8,J9,J10])*kernel(J1+1,J2+1,J3+1,J4+1,J5+1,J6+1,J7+1,J8+1,J9+1,J10+1)...
                                                *p1...
                                                .*p2...
                                                .*p3...
                                                .*p4...
                                                .*p5...
                                                .*p6...
                                                .*p7...
                                                .*p8...
                                                .*p9...
                                                .*[zeros(J10,1);xn(1:end-J10)];
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    out=temp;
    
end 




function c=consts(ord,ordfact,v)
i=1;
while ~isempty(v)
    L=length(v);
    v=v(find(v(1)~=v));
    gruppi(i)=factorial(L-length(v));
    i=i+1;
end
c=ordfact/prod(gruppi);

