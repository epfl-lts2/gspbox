% sgwt_filter_design : Return list of scaled wavelet kernels and derivatives
%
% g{1} is scaling function kernel,
% g{2} ... g{Nscales+1} are wavelet kernels
%
% function [g,t]=sgwt_filter_design(lmax,Nscales,varargin)
%
% Inputs :
% lmax - upper bound on spectrum
% Nscales - number of wavelet scales
%
% selectable parameters :
% designtype
% lpfactor - default 20. lmin=lmax/lpfactor will be used to determine
%            scales, then scaling function kernel will be created to
%            fill the lowpass gap.
%
% Outputs :
% g - wavelet and scaling function kernel function handles
% t - values of scale parameters used for wavelet kernels

% This file is part of the SGWT toolbox (Spectral Graph Wavelet Transform toolbox)
% Copyright (C) 2010, David K. Hammond.
%
% The SGWT toolbox is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% The SGWT toolbox is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with the SGWT toolbox.  If not, see <http://www.gnu.org/licenses/>.

function [g,t]=sgwt_filter_design(lmax,Nscales,varargin)
control_params={'designtype','abspline3','lpfactor',20,...
    'a',2,...
    'b',2,...
    't1',1,...
    't2',2,...
    };
argselectAssign(control_params);
argselectCheck(control_params,varargin);
argselectAssign(varargin);

switch designtype
    case 'abspline3'
        lmin=lmax/lpfactor;
        t=sgwt_setscales(lmin,lmax,Nscales);
        gl = @(x) exp(-x.^4);
        glp = @(x) -4*x.^3 .*exp(-x.^4);
        gb= @(x) sgwt_kernel_abspline3(x,a,b,t1,t2);
        for j=1:Nscales
            g{j+1}=@(x) gb(t(j)*x);
        end
        % find maximum of g's ...
        % I could get this analytically as it is a cubic spline, but
        % this also works.
        f=@(x) -gb(x);
        xstar=fminbnd(f,1,2);
        gamma_l=-f(xstar);
        lminfac=.6*lmin;
        g{1}=@(x) gamma_l*gl(x/lminfac);
        
    case 'mexican_hat'
        lmin=lmax/lpfactor;
        t=sgwt_setscales(lmin,lmax,Nscales);
        gb=@(x) x.*exp(-x);
        gl = @(x) exp(-x.^4);
        for j=1:Nscales
            g{j+1}=@(x) gb(t(j)*x);
        end
        lminfac=.4*lmin;
        g{1}=@(x) 1.2*exp(-1)*gl(x/lminfac);
         
    case 'meyer'
        t=(4/(3*lmax)) * 2.^(Nscales-1:-1:0);
        g{1}= @(x) sgwt_kernel_meyer(t(1)*x,'sf');
        for j=1:Nscales
            g{j+1}= @(x) sgwt_kernel_meyer(t(j)*x,'wavelet');
        end
        
    case 'simple_tf'
        t=(1/(2*lmax)) * 2.^(Nscales-1:-1:0);
        g{1}= @(x) sgwt_kernel_simple_tf(t(1)*x,'sf');
        for j=1:Nscales
            g{j+1}= @(x) sgwt_kernel_simple_tf(t(j)*x,'wavelet');
        end

    otherwise        
        error('Unknown design type');
end
