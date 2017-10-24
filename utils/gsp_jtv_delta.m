function jd = gsp_jtv_delta( G,vertex,time,param )
%GSP_JTV_DELTA Generate time-vertex signal with delta at given time-vertex locations
%      Usage:  [ jd ] = gsp_jtv_delta( G,vertex,time )
%                     = gsp_jtv_delta( G,vertex,time,param )
%
%   Input parameters:
%           G               : Time-Vertex Graph or Cell or Array of two element containing the sizes of the 2d delta signal
%           vertex          : locations on the graph
%           time            : locations in time
%           param           : Structure of optional parameters
%   Output parameters:
%         	jd              : joint delta
%
%   Additional parameters
%   ---------------------
%   * *param.lag*    : If 1 the size of the delta in time is $2T-1$ to take into account negative time location.
%   * *param.concat* : If 1 concatenates 2d-delta in a 3d matrix (N1,N2,#loc) instead of summing them. (default 0)

%   Author: Francesco Grassi
%   Date: July 2016

if nargin<4
    param=struct;
end

if ~isfield(param,'lag'), param.lag=0;end;
if ~isfield(param,'concat'), param.concat=0;end;

if isstruct(G)
      if or(~isfield(G.jtv,'T'),~isfield(G.jtv,'fs'));
        error('GSP_JTV_DELTA need time dimension. Use GSP_JTV_GRAPH.')
      end
        N = G.N;
        T = G.jtv.T;
elseif iscell(G)
    N = G{1};
    T = G{2};
elseif isvector(G)
    N = G(1);
    T = G(2);
end

if param.lag
    time=time+T;
    T=2*T-1;
end

if max(vertex)>N || max(time)>T
    error('Vertex or time index greater than time-vertex graph size.')
end

n = numel(vertex);
t = numel(time);

if n~=t
    error('vertex and time indexes must have same sizes')
end

if param.concat
    jd = zeros(N,T,n);
    
    for ii = 1:n;
        jd(vertex(ii),time(ii),ii) = 1;
    end
else
    jd = zeros(N,T);
    
    for ii = 1:n;
        jd(vertex(ii),time(ii)) = 1;
    end
end

end

