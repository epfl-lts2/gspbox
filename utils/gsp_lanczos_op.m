function c = gsp_lanczos_op(G,fi,s,param)
%GSP_LANCZOS_OP Perform the lanczos approximation of the signal s

if nargin < 4, param = struct; end
if ~isfield(param,'verbose'), param.verbose = 1; end;
if ~isfield(param,'order'), param.order = 30; end


Nf = numel(fi);
Nv = size(s,2);
c = zeros(G.N*Nf,Nv);

for jj = 1:Nv
    
    if sum(abs(s(:,jj)))>eps
        [V,H] = lanczos(G.L, param.order, s(:,jj));

        [Uh, Eh] = eig(H);


        Eh = diag(Eh);
        Eh(Eh<0) = 0;
        fie = gsp_filter_evaluate(fi,Eh);
        V = V*Uh;


        for ii=1:Nf
           c((1:G.N) + G.N*(ii-1),jj) = V * (fie(:, ii) .* (V'*s(:,jj)));
        end
    else
        c(:,jj) = 0;
    end
end

% [V,H] = lanczos(G.L, k, s);
% 
% for jj = 1:Nv
%     ind = (1:size(H,1))+size(H,1)*(jj-1);
%     [Uh, Eh] = eig(H(:,ind));
% 
% 
%     Eh = diag(Eh);
%     fie = gsp_filter_evaluate(fi,Eh);
%     V(:,ind) = V(:,ind)*Uh;
% 
% 
%     for ii=1:Nf
%        c((1:G.N) + G.N*(ii-1),jj) = V(:,ind) * fie(:, ii) .* (V(:,ind)'*s(:,jj));
%     end
% end

end




function [V,H,orth] = lanczos(A,order,x)

[N,M] = size(x);

% normalization
norm2vec = @(x) (sum(x.^2,1)).^0.5;
q = x./repmat(norm2vec(x),N,1);

% Initialization
hiv =0:order:(order*M-1); % helping indice vector

V = zeros(N,M*order);
V(:,1+hiv) = q;


H = zeros(order+1,order*M);

r = A*q;
H(1,1+hiv) = sum(q .* r, 1 );
r = r - repmat(H(1,1+hiv),N,1).*q; 
H(2,1+hiv) = norm2vec(r);

if (nargout > 2)
    orth = zeros(M,1);
    orth(1) = norm(V'*V - M);
end

for k = 2:order
    
    if (sum(abs(H(k,k-1+hiv))) <= eps)
        H = H(1:k-1,sum_ind(1:k-1,hiv));
        V = V(:,sum_ind(1:k-1,hiv));
        if (nargout > 2)
            orth = orth(1:k-1);
        end
        return;
    end
    
    H(k-1,hiv+k) = H(k,hiv+k-1);
    v = q;
    q = r./repmat(H(k-1,k+hiv),N,1);
    V(:,k+hiv) = q;
  
    r = A*q;
    r = r - repmat(H(k-1,k+hiv),N,1).*v;
    H(k,k+hiv) = sum(q .* r, 1 );
    
    r = r - repmat(H(k,k+hiv),N,1).*q;
    % The next line has to be checked
    r = r - V*(V'*r); % full reorthogonalization
    H(k+1,k+hiv) = norm2vec(r);
    
    if (nargout > 2)
        orth(k) = [orth, norm(V'*V - M)];
    end
end
   
H = H(1:order,1:order);


% H = zeros(order+1,order);
% 
% q = x/norm(x);
% V(:,1) = q;
%  
% r = A*q;
% 
% H(1,1) = q'*r;
% r = r - H(1,1)*q ; 
% H(2,1) = norm(r);
% 
% if (nargout > 2)
%  orth = norm(V'*V - 1);
% end
% 
% for k = 2:order
%     
%     if (abs(H(k,k-1)) <= 1e-15)
%         H = H(1:k-1,1:k-1);
%         return;
%     end
%     
%     H(k-1,k) = H(k,k-1);
%     v = q;
%     q = r/H(k-1,k);
%     V(:,k) = q;
%   
%     r = A*q;
%     r = r - H(k-1,k)*v;
%     H(k,k) = q'*r;
%     
%     r = r - H(k,k)*q;
%     r = r - V*(V'*r); % full reorthogonalization
%     H(k+1,k) = norm(r);
%     
%     if (nargout > 2)
%         orth = [orth, norm(V'*V - eye(k))];
%     end
% end
%    
% H = H(1:order,1:order);
% 
end

function ind = sum_ind(ind1,ind2)
    ind = repmat(ind1(:),1,numel(ind2)) + repmat((ind2(:))',numel(ind1),1);
    ind = (ind(:))';
end
