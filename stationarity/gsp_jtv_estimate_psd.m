function [psd, ft] = gsp_jtv_estimate_psd(G, X, param)
%GSP_JTV_ESTIMATE_PSD Estimate the PSD of a time vertex process
%   Usage:  psd = gsp_jtv_estimate_psd(G, X)
%           psd = gsp_jtv_estimate_psd(G, X, param)
%
%   Input parameters:
%         G          : Graph
%         X          : Signal(s) (matrix of NxTxR)
%         param      : Optional parameters
%   Output parameters:
%         psd        : PSD matrix
%         ft         : Filter type
%
%   This function estimates the PSD for vertex time processes. It is a
%   generalization of the Bartlett method.
%
%   Additional parameters
%   ---------------------
%
%   * *param.estimator*  : Method used for the estimation. By default 'TA'
%     - 'TA': Average over the time domain.
%     - 'FC': Fourier convolution a trick used to reduce the computation. 
%     - 'VA': Average over the vertex domain
%     - 'TVA': Average over the time and the vertex domain
%     - 'FA' : Average on multiple realization in the joint spectral domain
%   * *param.L*       : Length of the window for the 'GSTFT' method
%   * *param.kernel*  : Kernel for the convolution for 'FC' method.
%
%
%   References: perraudin2016stationary

% Author : Nathanael Perraudin
% Date: 6 January 2016
% Testing: test_gsp_estimate_vertex_time_psd

if nargin<3
    param = struct;
end
if ~isfield(param,'estimator')
    param.estimator = 'FA';
end

% number of realizations
[N, T, R] = size(X);

switch param.estimator
    
    % This method compute the Fourier transform of the signals and average
    % the modulus
    case 'FA'
        if ~isfield(param, 'L')    
            param.L = compute_default_L(X,T,R);
        end
        Xhat = gsp_jft(G,reshape(X,N,param.L,R*T/param.L));
        psd = mean(abs(Xhat).^2,3);
        ft = 'js';
       
    
    % This method splits time in windows, and takes an average of the psd
    % over each time window. The size of each window (L) should be 
    % equal to the maximum correlation distance in the process 
    case 'TA'
        
        if ~gsp_check_fourier(G)
            error('The Fourier basis is needed the Fourier basis')
            %G = gsp_compute_fourier_basis(G);
        end

        if ~isfield(param, 'L')    
            param.L = compute_default_L(X,T,R);
        end        
        if ~isfield(param, 'a'),       param.a = round(param.L/2); end
        if ~isfield(param, 'M'),       param.M = G.jtv.T; end
        if ~isfield(param,'win_type'), param.win_type = 'itersine'; end
        
        if round(T/param.L) - (T/param.L) ~= 0
            warning('Optimally the length of the signal should be divisible by param.L.')            
            while round(T/param.L) - (T/param.L) ~= 0
                param.L = param.L + 1;
            end
        end

        % Compute a Gabor window
        w = gabwin(param.win_type, param.a, param.L);

        % estimate the psd as the average over each realization
        psd = zeros(G.N, G.jtv.T);
        for r = 1:R
            
            % do a windowed JFT
            coeff = gsp_jtwgft(G, w, X(:,:,r), param);

            % compute the psd for each window
            psd_win = abs(coeff(:, 1:end-1, :)).^2;
            
            % average all windows to obtain an estimate of the psd
            psd_est = squeeze(mean(psd_win, 2));
            
            % average over realizations
            psd = psd + transpose(psd_est);
        end

        psd = psd / R;
        
        % normalize energy
        psd = psd * (norm(X(:), 'fro')^2/R / sum(psd(:)));
        
    % Go to the joint frequency domain and smooth out the PSD with a kernel.     
        ft = 'js';

    case 'FC'
        if ~isfield(param, 'kernel')
            param.kernel = exp(-(-20:20).^2/3);
        end
        h = param.kernel;
        
        if ~gsp_check_fourier(G)
            G = gsp_compute_fourier_basis(G);
        end
                        
        Xhat2 = abs(gsp_jft(G,X)).^2;
        

        psd = conv2(h', h', mean(Xhat2,3), 'same') / norm(h, 1)^2;
%         psd = conv2(mean(Xhat2,3), param.kernel, 'same') / norm(param.kernel, 1);


        % normalize energy
        psd = psd * (norm(X(:), 'fro')^2/R / sum(psd(:)));
        ft = 'js';

    % For scalability, 
    case 'VA'

        if ~isfield(G,'boundary'); param.boundary = 'periodic'; end

        switch param.boundary
        
            case 'periodic'
                Xhat = fft(X,[],2)/sqrt(T);
            
            case 'reflecting'
                error('Sorry. Not implemented yet...')
                
            otherwise
                error('Unknown boundary condition');
        end
        
        if ~isfield(G,'lmax')
            G = gsp_estimate_lmax(G);
            
            warning(['GSP_PSD_ESTIMATION: The variable lmax is not ',...
                'available. The function will compute it for you. ',...
                'However, if you apply many time this function, you ',...
                'should precompute it using the function: ',...
                'gsp_estimate_lmax']);
        end
        
        if ~isfield(param,'Nfilt'),   param.Nfilt = 50; end
        if ~isfield(param,'Nrandom'), param.Nrandom = max(10, R); end
        if ~isfield(param,'g0')
            sigma = sqrt(2*G.lmax/param.Nfilt^2 * (param.Nfilt + 1));
            param.g0 = @(x) exp(-x.^2/sigma^2);
        end
        
        % Design the frame
        [ g , mu ] = gsp_design_translates(G, param.g0, param.Nfilt);
        %[ g , mu ] = gsp_design_itersine(G,Nfilt );
        
        % Estimate the energy of the window
        if gsp_check_fourier(G)
            mu_y2 = sum(abs(gsp_filter_evaluate(g,G.e)).^2,1)';
        else
            w = randn(N, param.Nrandom);
            x_filt2    = gsp_vec2mat( gsp_filter(G,g,w,param), param.Nfilt);
            n2_x_filt2 = sum(abs(x_filt2).^2,1);
            mu_y2      = reshape(mean(n2_x_filt2,3),[],1);
        end
        
        psd = cell(1,T);
        
        for ii = 1:T
            % Perform the filtering
            x_filt = gsp_vec2mat( ...
                gsp_filter(G, g, squeeze(Xhat(:,ii,:)), param), ...
                param.Nfilt);
            % estimate the points
            n2_x_filt = sum(abs(x_filt).^2,1);
            mu_y = reshape(mean(n2_x_filt,3),[],1);
            
            % Interpolate to obtain a nice filter.
            psd{ii} = @(s) max(spline(mu,mu_y./mu_y2,s), 0);
        end
        ft = 'js-array';

        
   case 'TVA'

        if ~isfield(param, 'L') 
            
            threshold = 0.05;
            
            % compute the average correlation distance
            C_T = zeros(T); Toep = toeplitz(0:(T-1)); Tmax = round(0.5*T);

            for r = 1:R, C_T = C_T + abs(X(:,:,r)'*X(:,:,r)) / R; end
            cor = zeros(Tmax,1);
            for t = 1:Tmax, cor(t) = mean( vec(C_T(Toep == t-1)) ); end
            cor = normalize_data(cor, 0, 1);
            
            % select the window length conservatively
            param.L = min(4*find(cor>=threshold, 1, 'last' ), T);

            % visualize the correlation 
            % figure; plot(0:(Tmax-1), cor, '-o', [1 1]*param.L, [0 1],
            % 'r-'); xlabel('correlation distance'); 
            
            % param.L = round(T/4);
        end        
        
       % Make sure that the length of the signal should be divisible by param.L  
       if round(T/param.L) - (T/param.L) ~= 0
           % warning('Ideally the length of the signal should be divisible by param.L.')           
           while (round(T/param.L) - (T/param.L) ~= 0) || (round(param.L/2) - (param.L/2) ~= 0)
              param.L = param.L + 1;  
           end
       end
       
       if ~isfield(param, 'a'),       param.a = round(param.L/2);  end
       if ~isfield(param, 'M'),       param.M = G.jtv.T;           end
       if ~isfield(param,'win_type'), param.win_type = 'itersine'; end
       if ~isfield(G,'boundary');     param.boundary = 'periodic'; end
        
              
       % Compute a Gabor window
       w = gabwin(param.win_type, param.a, param.L);

        switch param.boundary
            
            case 'periodic'
                S = zeros(N, param.M, ceil(T/param.a), R);
                for r = 1:R
                    % do a discrete gabore transform
                    X_gabore = dgt(transpose(X(:,:,r)), w, param.a, param.M);
                    S(:,:,:,r) = permute( X_gabore, [3,1,2]);
                end                
            case 'reflecting'
                error('Sorry. Not implemented yet...')
            otherwise
                error('Unknown boundary condition');
        end
        
        if ~isfield(G,'lmax')
            G = gsp_estimate_lmax(G);
            
            warning(['GSP_PSD_ESTIMATION: The variable lmax is not ',...
                'available. The function will compute it for you. ',...
                'However, if you apply many time this function, you ',...
                'should precompute it using the function: ',...
                'gsp_estimate_lmax']);
        end
        
        if ~isfield(param,'Nfilt'),   param.Nfilt = 50;          end
        if ~isfield(param,'Nrandom'), param.Nrandom = max(10,R); end
        if ~isfield(param,'g0')
            sigma = sqrt(2*G.lmax/param.Nfilt^2 * (param.Nfilt + 1));
            param.g0 = @(x) exp(-x.^2/sigma^2);
        end
        
        % Design the frame
        [ g , mu ] = gsp_design_translates(G, param.g0, param.Nfilt);
        %[ g , mu ] = gsp_design_itersine(G,Nfilt );
        
        % Estimate the energy of the window
        if gsp_check_fourier(G)
            mu_y2 = sum(abs(gsp_filter_evaluate(g, G.e)).^2, 1)';
        else
            w = randn(N, param.Nrandom);
            x_filt2 = gsp_vec2mat( gsp_filter_analysis(G,g,w,param), param.Nfilt);
            n2_x_filt2 = sum(abs(x_filt2).^2, 1);
            mu_y2 = reshape(mean(n2_x_filt2,3), [], 1);
        end
        
        psd = cell(1, param.M);
        
%         if gsp_check_fourier(G) && param.use_fast
%             warning('This may use a lot of ram!')
%             [N, T, Nt, Ns ] = size(S);
%             Nf = numel(g);
%             ShatG = gsp_gft(G,S);
%             filt_eval = gsp_filter_evaluate(g, G.e);
%             filt_eval = repmat(reshape(filt_eval,[N,1,1,1,Nf]), 1,T,Nt,Ns ,1);
%             % Perform the filtering
%             S_filt = gsp_igft(G, repmat(ShatG,[1,1,1,1,Nf]) .* filt_eval);
%             n2_S_filt = squeeze(mean(mean(sum(abs(S_filt).^2,1),3),4));
%             for ii = 1:param.M           
%                 
%                 % Interpolate to obtain a nice filter.
%                 psd{ii} = @(s) max(spline(mu,reshape(n2_S_filt(ii,:),[],1)./mu_y2,s), 0);
%             end
%         else
        for ii = 1:param.M
            
            % Perform the filtering
            x_filt = gsp_vec2mat( ...
                gsp_filter(G, g, reshape(S(:,ii,:,:),N,[]), param), ...
                param.Nfilt);
            
            % estimate the points
            n2_x_filt = sum(abs(x_filt).^2,1);
            mu_y = reshape(mean(n2_x_filt,3),[],1);
            
            % Interpolate to obtain a nice filter.
            psd{ii} = @(s) max(spline(mu,mu_y./mu_y2,s), 0);
        end
%         end

        
        ft = 'js-array';

         

    otherwise
        error('Unknown method')
end



end


function [x] = normalize_data(x, xmin, xmax)
%NORMALIZE_VECTOR Normalize x between xmin and xmax
%   x might be a scalar, vector, or matrix.

%normalize to [0,1]
x = (x - min(min(x))) ./ (max(max(x)) - min(min(x)));

if exist('xmax', 'var') && exist('xmin', 'var')
   x = x .* (xmax - xmin) + xmin; 
end

end

function L =  compute_default_L(X,T,R)
% compute the correlation distance
C_T = zeros(T);
for r = 1:R
    C_T = C_T + abs(X(:,:,r)'*X(:,:,r)) / R;
end
Toep = toeplitz(0:(T-1));
Tmax = round(0.5*T);
cor = zeros(Tmax,1);
for t = 1:Tmax
    cor(t) = mean( vec(C_T(Toep == t-1)) ); 
end
cor = normalize_data(cor, 0, 1);
threshold = 0.05;
L = min(2*find(cor>=threshold, 1, 'last'), T);

%L = round(T/4);
%figure; plot(0:(Tmax-1), cor, '-o', [1 1]*param.L, [0 1], 'r-');

end