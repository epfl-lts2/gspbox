function fd = gsp_jtv_filter_evaluate_simple(G,g,x,param)

if nargin<4
    param = struct;
end

if ~iscell(g)
    g = {g};
end

if isnumeric(g{1})
    % for different filter, the user might just give fd
    fd=zeros(size(g{1},1),size(g{1},2),numel(g));
    for ii=1:numel(g)
        fd(:,:,ii)=g{ii};
    end
else
    
    if ~isfield(param,'filtertype')
        param.filtertype = 'time-spectral';
    end
    
    
    
    N = G.N;
    T = G.jtv.T;
    
    if nargin(g{1})==1
        Nf = size(g,1);
        M = size(g,2);
        fd=zeros(N,M,Nf);
        
        
        switch param.filtertype
            case 'ts'
                for ii=1:Nf
                    for jj = 1:M
                        fd(:,jj,ii)=g{ii,jj}(x(:));
                    end
                end
                fd = fft(fd,[],2)/sqrt(size(fd,2));
            case 'js'
                for ii=1:Nf
                    for jj = 1:M
                        fd(:,jj,ii)=g{ii,jj}(x(:));
                    end
                end
            otherwise
                error('Unknown filter type')
        end
    elseif nargin(g{1})==2
        Nf = numel(g);
        switch param.filtertype
            case 'ts'
                [xx,tt]=meshgrid(x,0:1/G.jtv.fs:(T-1)/G.jtv.fs);
                fd=zeros(N,T,Nf);
                
                for ii=1:Nf
                    fd(:,:,ii)=transpose(g{ii}(xx,tt));
                end
                fd = fft(fd,[],2)/sqrt(size(fd,2));
            case 'js'
                [xx,ff]=meshgrid(x,gsp_jtv_fa(G));
                fd=zeros(N,T,Nf);
                
                for ii=1:Nf
                    fd(:,:,ii)=transpose(g{ii}(xx,ff));
                end
            otherwise
                error('Unknown filter type')
        end
        
        
    else
        error('Two many variable in the filter')
    end
    
    
end



end