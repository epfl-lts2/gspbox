function filter = gsp_jtv_filter_evaluate(g,filtertype, x, v, param)
%GSP_JTV_FILTER_EVALUATE Evaluate the time-vertex filters
%   Usage: filter = gsp_jtv_filter_evaluate(g,filtertype, x, v)
%          filter = gsp_jtv_filter_evaluate(g,filtertype, x, v, param)
%
%   Input parameters:
%       g          : Cell array of time-vertex filters
%       filtertype : Filter domain (ts,js,ts-array,js-array)
%       x          : Graph domain values (Eigenvalues)
%       v          : Time/Frequency domain values
%       param      : Optional parameters
%
%   Output parameters:
%       filter     : Response of the time-vertex filters (by default in the
%                    joint spectral domain)
%
%   This function apply all the filters in *g* to the data *x,v*. The dummy
%   variable *v* identifies time or frequency values for ts(-array) or
%   js(-array) filters respectively. The filter is defined in the
%   joint-spectral domain by default, i.e. filters $= g(\lambda,\omega)$.
%
%   Every time-vertex filter correspond to one matrix of the cube *filter*.
%
%   Additional parameters
%   ---------------------
%   * *param.domain*   : Specifies in which domain the function should
%                        return the filters: 'time-spectral' or
%                        'joint-spectral' (default param.domain='joint-spectral')
%   * *param.parfor*   : Parallelized computation of filters in g in case of
%                        ts,js-array filters
%

% Author: Francesco Grassi
% Date : September 2016
% Testing: test_jtv_filter


if ~gsp_check_filtertype(filtertype)
    error('Invalid filtertype');
end

if nargin<5
    param = struct;
end

if ~isfield(param,'domain'), param.domain = 'joint-spectral'; end
if ~isfield(param,'parfor'), param.parfor = 0; end

if ~iscell(g)
    g = {g};
end



if isnumeric(g{1})  %Case 1: The user provides numerical filters
    if size(g{1},2) == 1  %Case 2.1: ts/js-array filters (array on 1d filters)
        filter = zeros(size(g{1},1),size(g,2),size(g,1));
        for n=1:numel(g)
            filter(:,:,n)=cell2mat(g{n,:});
        end
        
    elseif size(g{1},2)>1  %Case 2.1: ts/js filters
        
        filter=zeros(size(g{1},1),size(g{1},2),numel(g));
        for n=1:numel(g)
            filter(:,:,n)=g{n};
        end
    else
        error('Invalid numerical joint time-vertex filter');
    end
    
else %Case 2: The user provides function handle filters
    
    N = length(x);
    Nt = length(v);
    
    switch filtertype
        
        case 'ts'
            
            Nf = numel(g);
            
            time    = v;
            [tt,xx] = meshgrid(time,x);
            filter  = zeros(N,Nt,Nf);
            
            for n=1:Nf
                filter(:,:,n) = g{n}(xx,tt);
            end
            if strcmpi(param.domain,'joint-spectral')
                filter = fft(filter,[],2)/sqrt(size(filter,2));
            end
        case 'js'
            
            Nf      = numel(g);
            
            omega   = v;
            [ww,xx] = meshgrid(omega,x);
            filter  = zeros(N,Nt,Nf);
            
            for n=1:Nf
                filter(:,:,n) = g{n}(xx,ww);
            end
            if strcmpi(param.domain,'time-spectral')
                filter = ifft(filter,[],2)*sqrt(size(filter,2));
            end
        case {'ts-array','js-array'}
            
            Nf = size(g,1);
            M = size(g,2);
            filter=zeros(N,M,Nf);
            
            %             for ii=1:Nf
            %                 fd(:,:,ii) = cell2mat(cellfun(@(h) h(x),g(ii,:),'uniformoutput',0));
            %             end
            
            if param.parfor
                parfor ii=1:Nf
                    for jj = 1:M
                        filter(:,jj,ii) = g{ii,jj}(x);
                    end
                end
            else
                for ii=1:Nf
                    for jj = 1:M
                        filter(:,jj,ii) = g{ii,jj}(x);
                    end
                end
            end
            if and(strcmpi(param.domain,'time-spectral'),strcmpi(filtertype,'js-array'))
                filter = ifft(filter,[],2)*sqrt(size(filter,2));
            elseif and(strcmpi(param.domain,'joint-spectral'),strcmpi(filtertype,'ts-array'))
                filter = fft(filter,[],2)/sqrt(size(filter,2));
            end
            
    end
    
end





end
