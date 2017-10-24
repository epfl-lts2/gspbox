function errors = test_filter()
%TEST_FILTERS This function test all the filters

errors = 0;

errors = errors + test_mexican_hat1();
errors = errors + test_mexican_hat2();
errors = errors + test_meyer();
errors = errors + test_abspline();
errors = errors + test_simple_tf();
errors = errors + test_itersine();
errors = errors + test_uniform_half_cosine();
errors = errors + test_design_warped_translates();
errors = errors + test_held();
errors = errors + test_simoncelli();
errors = errors + test_papadakis();
errors = errors + test_regular();


errors = errors + test_filter_lanczos();
errors = errors + test_filter_lanczos_multiple();

errors = errors + test_filter_lanczos_syn();
errors = errors + test_filter_lanczos_syn_multiple();


errors = errors + test_cheb_coeff();
errors = errors + test_cheb_op();
errors = errors + test_filter_analysis();
errors = errors + test_filter_synthesis();
errors = errors + test_filter_inverse();
errors = errors + test_localization();

errors = errors + test_multiple_dimensions();

errors = errors + test_reshape();
errors = errors + test_matrix_op();


try  %#ok<TRYNC>
    close(100)
end

end



function errors = test_regular()

errors = 0;
try
  
   figure(100);
   G = gsp_sensor(100);
   G = gsp_estimate_lmax(G);
   g = gsp_design_regular(G);
   gsp_plot_filter(G,g);
   close(100);


   fprintf('FILTER: regular kernel 1 ok\n');
catch
    errors = errors + 1;
    warning('FILTER: Error kernel regular 1 test')
end

try
  
   figure(100);
   G = gsp_sensor(100);
   G = gsp_estimate_lmax(G);
   param.d = 5;
   g = gsp_design_regular(G,param);
   gsp_plot_filter(G,g);
   close(100);


   fprintf('FILTER: regular kernel 2 ok\n');
catch
    errors = errors + 1;
    warning('FILTER: Error kernel regular 2 test')
end


try
  
   figure(100);
   G = gsp_sensor(100);
   G = gsp_estimate_lmax(G);
   param.d = 0;
   g = gsp_design_regular(G,param);
   gsp_plot_filter(G,g);
   close(100);


   fprintf('FILTER: regular kernel 3 ok\n');
catch
    errors = errors + 1;
    warning('FILTER: Error kernel regular 3 test')
end

end


function errors = test_papadakis()

errors = 0;
try
  
   figure(100);
   G = gsp_sensor(100);
   G = gsp_estimate_lmax(G);
   g = gsp_design_papadakis(G);
   gsp_plot_filter(G,g);
   close(100);


   fprintf('FILTER: papadakis kernel 1 ok\n');
catch
    errors = errors + 1;
    warning('FILTER: Error kernel papadakis 1 test')
end

try
  
   figure(100);
   G = gsp_sensor(100);
   G = gsp_estimate_lmax(G);
   param.a = 0.25;
   g = gsp_design_papadakis(G,param);
   gsp_plot_filter(G,g);
   close(100);


   fprintf('FILTER: papadakis kernel 2 ok\n');
catch
    errors = errors + 1;
    warning('FILTER: Error kernel papadakis 2 test')
end

end


function errors = test_simoncelli()

errors = 0;
try
  
   figure(100);
   G = gsp_sensor(100);
   G = gsp_estimate_lmax(G);
   g = gsp_design_simoncelli(G);
   gsp_plot_filter(G,g);
   close(100);


   fprintf('FILTER: simoncelli kernel 1 ok\n');
catch
    errors = errors + 1;
    warning('FILTER: Error kernel simoncelli 1 test')
end

try
  
   figure(100);
   G = gsp_sensor(100);
   G = gsp_estimate_lmax(G);
   param.a = 0.25;
   g = gsp_design_simoncelli(G,param);
   gsp_plot_filter(G,g);
   close(100);


   fprintf('FILTER: simoncelli kernel 2 ok\n');
catch
    errors = errors + 1;
    warning('FILTER: Error kernel simoncelli 2 test')
end

end


function errors = test_held()

errors = 0;
try
  
   figure(100);
   G = gsp_sensor(100);
   G = gsp_estimate_lmax(G);
   g = gsp_design_held(G);
   gsp_plot_filter(G,g);
   close(100);


   fprintf('FILTER: held kernel 1 ok\n');
catch
    errors = errors + 1;
    warning('FILTER: Error kernel held 1 test')
end

try
  
   figure(100);
   G = gsp_sensor(100);
   G = gsp_estimate_lmax(G);
   param.a = 0.25;
   g = gsp_design_held(G,param);
   gsp_plot_filter(G,g);
   close(100);


   fprintf('FILTER: held kernel 2 ok\n');
catch
    errors = errors + 1;
    warning('FILTER: Error kernel held 2 test')
end

end


function errors = test_design_warped_translates()

errors = 0;
try
  
    figure(100);
    Nf = 10;
    G = gsp_sensor(100);
    G = gsp_estimate_lmax(G);
    G = gsp_spectrum_cdf_approx(G);
    g = gsp_design_warped_translates(G, Nf);   
    gsp_plot_filter(G,g);
    close(100);
    
   fprintf('FILTER: test_design_warped_translates 1 ok\n');
catch
    errors = errors + 1;
    warning('FILTER: Error in test_design_warped_translates 1 test')
end

try
  
    figure(100);
    Nf = 10;
    G = gsp_sensor(100);
    G = gsp_estimate_lmax(G);
    param.log =1 ;
    param.warping_type = 'spectrum_approximation';
    G = gsp_spectrum_cdf_approx(G);
    g = gsp_design_warped_translates(G, Nf,param);   
    gsp_plot_filter(G,g);
     close(100);
    
   fprintf('FILTER: test_design_warped_translates 2 ok\n');
catch
    errors = errors + 1;
    warning('FILTER: Error in test_design_warped_translates 2 test')
end

try
  
    figure(100);
    Nf = 10;
    G = gsp_sensor(100);
    G = gsp_estimate_lmax(G);
    param.log =1 ;
    param.warping_type = 'custom';
    G = gsp_spectrum_cdf_approx(G);
    g = gsp_design_warped_translates(G, Nf,param);   
    gsp_plot_filter(G,g);
    close(100);
    
   fprintf('FILTER: test_design_warped_translates 3 ok\n');
catch
    errors = errors + 1;
    warning('FILTER: Error in test_design_warped_translates 3 test')
end

try
  
    figure(100);
    Nf = 10;
    G = gsp_sensor(100);
    G = gsp_estimate_lmax(G);
    param.log =0 ;
    param.warping_type = 'custom';
    G = gsp_spectrum_cdf_approx(G);
    g = gsp_design_warped_translates(G, Nf,param);   
    gsp_plot_filter(G,g);
    close(100);
    
   fprintf('FILTER: test_design_warped_translates 4 ok\n');
catch
    errors = errors + 1;
    warning('FILTER: Error in test_design_warped_translates 4 test')
end


try
  
    figure(100);
    Nf = 10;
    G = gsp_sensor(100);
    G = gsp_compute_fourier_basis(G);
    tol=1e-8;
    [unique_E,unique_E_inds]=unique(round(G.e*1/tol)*tol);
    param.approx_spectrum.x=unique_E;
    param.approx_spectrum.y=(unique_E_inds-1)/(G.N-1);  
    param.interpolation_type='monocubic';
    param.log =0 ;
    param.warping_type = 'spectrum_interpolation';
    g = gsp_design_warped_translates(G, Nf,param);   
    gsp_plot_filter(G,g);
    close(100);
    
   fprintf('FILTER: test_design_warped_translates 5 ok\n');
catch
    errors = errors + 1;
    warning('FILTER: Error in test_design_warped_translates 5 test')
end

try
  
    figure(100);
    Nf = 10;
    G = gsp_sensor(100);
    G = gsp_compute_fourier_basis(G);
    tol=1e-8;
    approx_spectrum_inds=1:11:100;
    N_bar=length(approx_spectrum_inds);
    [unique_subset_E, unique_subset_E_inds]=unique(round(G.e(approx_spectrum_inds)*1/tol)*tol);
    param.approx_spectrum.x=unique_subset_E;
    param.approx_spectrum.y=(unique_subset_E_inds-1)/(N_bar-1);  
    param.interpolation_type='pwl';
    param.log =1 ;
    param.warping_type = 'spectrum_interpolation';
    G = gsp_spectrum_cdf_approx(G);
    g = gsp_design_warped_translates(G, Nf,param);   
    gsp_plot_filter(G,g);
    close(100);
    
   fprintf('FILTER: test_design_warped_translates 6 ok\n');
catch
    errors = errors + 1;
    warning('FILTER: Error in test_design_warped_translates 6 test')
end

end





function errors = test_mexican_hat1()

errors = 0;
try
  
   G = gsp_sensor(100);
   gsp_design_mexican_hat(G, 5);



   fprintf('FILTER: mexican hat 1 ok\n');
catch
    errors = errors + 1;
    warning('FILTER: Error in the mexican hat 1 test')
end

end

function errors = test_mexican_hat2()

errors = 0;
try
  
    figure(100);
    Nf = 4;
    G = gsp_sensor(100);
    G = gsp_estimate_lmax(G);
    g = gsp_design_mexican_hat(G, Nf);   
    gsp_plot_filter(G,g);
    close(100);
    
   fprintf('FILTER: mexican hat 2 ok\n');
catch
    errors = 1;
    warning('FILTER: Error in the mexican hat 2 test')
end

end


function errors = test_uniform_half_cosine()

errors = 0;
try
  
    figure(100);
    Nf = 4;
    G = gsp_sensor(100);
    G = gsp_estimate_lmax(G);
    g = gsp_design_half_cosine(G, Nf);   
    gsp_plot_filter(G,g);
    close(100);
    
   fprintf('FILTER: uniform half cosine ok\n');
catch
    errors = 1;
    warning('FILTER: Error in uniform half cosine test')
end

end

function errors = test_abspline()

errors = 0;
try
  
    figure(100);
    Nf = 4;
    G = gsp_sensor(100);
    G = gsp_estimate_lmax(G);
    g = gsp_design_abspline(G, Nf);   
    gsp_plot_filter(G,g); 
    close(100);
    
   fprintf('FILTER: abspline ok\n');
catch
    errors = 1;
    warning('FILTER: Error abspline test')
end

end


function errors = test_meyer()

errors = 0;
try
  
    figure(100);
    Nf = 4;
    G = gsp_sensor(100);
    G = gsp_estimate_lmax(G);
    g = gsp_design_meyer(G, Nf);   
    gsp_plot_filter(G,g);
    close(100);
    
   fprintf('FILTER: meyer ok\n');
catch
    errors = 1;
    warning('FILTER: Error in meyer test')
end

end


function errors = test_simple_tf()

errors = 0;
try
  
    figure(100);
    Nf = 4;
    G = gsp_sensor(100);
    G = gsp_estimate_lmax(G);
    g = gsp_design_simple_tf(G, Nf);   
    gsp_plot_filter(G,g);
    close(100);
    
   fprintf('FILTER: simple_tf ok\n');
catch
    errors = 1;
    warning('FILTER: Error in simple_tf test')
end

end


function errors = test_itersine()

errors = 0;
try
  
        Nf = 20;
        G = gsp_sensor(100);
        G = gsp_estimate_lmax(G);
        g = gsp_design_itersine(G, Nf);   
        gsp_plot_filter(G,g);  
        [A,B] = gsp_filterbank_bounds(G,g);
    
   fprintf('FILTER: itersine ok\n');
catch
    errors = errors + 1;
    warning('FILTER: Error in itersine test')
end

if norm(A-B)<eps(1000)
    fprintf('FILTER: itersine 2 ok\n');
else
    errors = errors + 1;
    warning('FILTER: Error in itersine 2 test')
end

end


function errors = test_cheb_coeff()

 errors = 0;
 try
   
   Nf = 4;
    G = gsp_sensor(100);
    G = gsp_estimate_lmax(G);
    g = gsp_design_meyer(G, Nf);  
    c = gsp_cheby_coeff(G, g);

     
  fprintf('FILTER: cheb coeff ok\n');

 catch
	 errors = 1;
	 warning('FILTER: Error in cheb coeff test')
end
 
end


function errors = test_cheb_op()

errors = 0;
try
  
        Nf = 4;
        G = gsp_sensor(100);
        G = gsp_estimate_lmax(G);
        g = gsp_design_meyer(G, Nf);  
        c = gsp_cheby_coeff(G, g);
        f = rand(G.N,1);
        r = gsp_cheby_op(G, c, f);
    
   fprintf('FILTER: cheb op ok\n');
catch
    errors = 1;
    warning('FILTER: Error in cheb op test')
end

end



function errors = test_filter_analysis()

    Nf = 4;
    G = gsp_sensor(10);
    G = gsp_estimate_lmax(G);
    G = gsp_compute_fourier_basis(G);
    g = gsp_design_meyer(G, Nf);  
    f = rand(G.N,1);
    
    param.method = 'exact';
    param.order = 100;
    f1 = gsp_filter_analysis(G,g,f,param);
    param.method = 'cheby';
    f2 = gsp_filter_analysis(G,g,f,param);


    if norm(f1(:)-f2(:))/norm(f1) < 1e-4
        errors = 0;
       fprintf('FILTER: analysis ok\n');
    else
        errors = 1;
        warning('FILTER: Error in analysis test')
    end

end


function errors = test_filter_synthesis()

    Nf = 4;
    G = gsp_sensor(10);
    G = gsp_estimate_lmax(G);
    G = gsp_compute_fourier_basis(G);
    g = gsp_design_simple_tf(G, Nf);  
    f = rand(G.N,1);
    f = f/norm(f);
    
    param.method = 'cheby';
    param.order = 1000;
    f1 = gsp_filter_analysis(G,g,f,param);
    f2 = gsp_filter_synthesis(G,g,f1,param);
    param.method = 'exact';
    f3 = gsp_filter_synthesis(G,g,f1,param);
    f4 = f2/norm(f2);
    
    if norm(f-f4) < 1e-5
        errors = 0;
       fprintf('FILTER: synthesis ok\n');
    else
        errors = 1;
        norm(f-f4)
        warning('FILTER: Error in synthesis test')
    end

    if norm(f3-f2)/norm(f2) < 1e-5
       fprintf('FILTER: synthesis 2 ok\n');
    else
        errors = errors+1;
        warning('FILTER: Error in synthesis 2 test')
    end
end


function errors = test_filter_inverse()

    Nf = 4;
     G = gsp_sensor(10);
%     G = gsp_estimate_lmax(G);
    G = gsp_compute_fourier_basis(G);
    
    g = gsp_design_mexican_hat(G, Nf);  
    f = rand(G.N,1);
    
    param.method = 'exact';
    f1 = gsp_filter_analysis(G,g,f,param);
    f2 = gsp_filter_inverse(G,g,f1,param);
    


    if norm(f-f2) < 1e-10
        errors = 0;
       fprintf('FILTER: inverse ok\n');
    else
        errors = 1;
        warning('FILTER: Error in inverse test')
    end

end


function errors = test_localization()

   % Number of nodes
    N = 100;

    % Kernel
    tau = 0.1;

    g = @(x) exp(-tau*x);


    G1 = gsp_sensor(N);


    G1 = gsp_estimate_lmax(G1);
    s1 = sqrt(G1.N)*gsp_filter_analysis(G1,g,eye(N));

    % Do the test
    G1 = gsp_compute_fourier_basis(G1);
    shat = g(G1.e);
    s = gsp_igft(G1,shat);
    
    s2 = zeros(G1.N);
    for ii = 1:G1.N
        s2(:,ii) = gsp_translate_old(G1, s, ii);
    end



    if norm(s1-s2) < 1e-10
        errors = 0;
       fprintf('FILTER: localisation ok\n');
    else
        errors = 1;
        warning('FILTER: Error in localisation test')
    end

end


function errors = test_multiple_dimensions()

   % Number of nodes
    N = 100;


    G1 = gsp_sensor(N);
    
    Nf = 4;
    k = 6;

    G1 = gsp_estimate_lmax(G1);
    g = gsp_design_meyer(G1, Nf);  


    G1 = gsp_estimate_lmax(G1);
    
    f = rand(N,k);
    
    s1 = gsp_filter_analysis(G1,g,f);
    
    s2 = zeros(G1.N*Nf,k);
    for ii = 1:k
        s2(:,ii) =  gsp_filter_analysis(G1,g,f(:,ii));
    end



    if norm(s1(:)-s2(:)) < 1e-10
        errors = 0;
       fprintf('FILTER: multiple_dimensions ok\n');
    else
        errors = 1;
        warning('FILTER: Error in multiple_dimensions test')
    end
    
    f1 = gsp_filter_synthesis(G1,g,s1);
    
    f2 = zeros(G1.N,k);
    for ii = 1:k
        f2(:,ii) =  gsp_filter_synthesis(G1,g,s1(:,ii));
    end

     if norm(f1(:)-f2(:)) < 1e-10
       fprintf('FILTER: multiple_dimensions 2 ok\n');
    else
        errors = errors + 1;
        warning('FILTER: Error in multiple_dimensions 2 test')
     end

        f1 = gsp_filter_inverse(G1,g,s1);
    
    f2 = zeros(G1.N,k);
    for ii = 1:k
        f2(:,ii) =  gsp_filter_inverse(G1,g,s1(:,ii));
    end

    if norm(f1(:)-f2(:)) < 1e-10
        fprintf('FILTER: multiple_dimensions 3 ok\n');
    else
        errors = errors + 1;
        warning('FILTER: Error in multiple_dimensions 3 test')
    end
    
%     if norm(f1(:)-f(:))/norm(f(:)) < 1e-6
%         fprintf('FILTER: multiple_dimensions 4 ok\n');
%     else
%         errors = errors + 1;
%         norm(f1(:)-f(:))/norm(f(:))
%         warning('FILTER: Error in multiple_dimensions 4 test')
%     end
    
end


function errors = test_reshape()

   % Number of nodes
    N = 10;
    M = 4;
    
    f = rand(M*N,N);
    
    ft = gsp_vec2mat(f,M);
    [f1, M2] = gsp_mat2vec(ft);

    errors = 0;
    
    if norm(f1(:)-f(:))/norm(f(:)) < 1e-12
        fprintf('FILTER: reshape ok\n');
    else
        errors = errors + 1;
        warning('FILTER: erros reshape test')
    end
    
    if norm(M2-M) < 1e-12
        fprintf('FILTER: reshape 2 ok\n');
    else
        errors = errors + 1;
        warning('FILTER: erros reshape test 2')
    end
end



function errors = test_matrix_op()

    errors = 0;
    % Number of nodes
    N = 100;


    G1 = gsp_sensor(N);
    
    Nf = 4;

    G1 = gsp_estimate_lmax(G1);
    g = gsp_design_meyer(G1, Nf);  


    G1 = gsp_estimate_lmax(G1);
    
    x = rand(N,1);
    
    f = gsp_filter_analysis(G1,g,x);
    F = gsp_filterbank_matrix(G1,g);
    f1 = F'*x;
    
    if norm(f1(:)-f(:))/norm(f(:)) < 1e-12
        fprintf('FILTER: operator matrix ok\n');
    else
        errors = errors + 1;
        warning('FILTER: error in operator matrix test')
    end
    
   
end



function errors = test_filter_lanczos()

    Nf = 4;
    G = gsp_sensor(10);
    G = gsp_estimate_lmax(G);
    G = gsp_compute_fourier_basis(G);
    g = gsp_design_meyer(G, Nf);  
    f = rand(G.N,1);
    
    param.method = 'exact';
    param.order = 100;
    f1 = gsp_filter_analysis(G,g,f,param);
    param.method = 'lanczos';
    f2 = gsp_filter_analysis(G,g,f,param);


    if norm(f1(:)-f2(:))/norm(f1(:)) < 1e-10
        errors = 0;
       fprintf('FILTER: lanczos ok\n');
    else
        errors = 1;
        norm(f1(:)-f2(:))/norm(f1(:))
        warning('FILTER: Error in lanczos test')
    end

end

function errors = test_filter_lanczos_multiple()

    Nf = 4;
    G = gsp_sensor(10);
    G = gsp_estimate_lmax(G);
    G = gsp_compute_fourier_basis(G);
    g = gsp_design_meyer(G, Nf);  
    f = rand(G.N,5);
    
    param.method = 'exact';
    param.order = 100;
    f1 = gsp_filter_analysis(G,g,f,param);
    param.method = 'lanczos';
    f2 = gsp_filter_analysis(G,g,f,param);


    if norm(f1(:)-f2(:))/norm(f1(:)) < 1e-10
        errors = 0;
       fprintf('FILTER: lanczos multiple ok\n');
    else
        errors = 1;
        norm(f1(:)-f2(:))/norm(f1(:))
        warning('FILTER: Error in lanczos multiple test')
    end

end



function errors = test_filter_lanczos_syn()

    Nf = 4;
    G = gsp_sensor(10);
    G = gsp_estimate_lmax(G);
    G = gsp_compute_fourier_basis(G);
    g = gsp_design_meyer(G, Nf);  
    f = rand(G.N*Nf,1);
    
    param.method = 'exact';
    param.order = 100;
    f1 = gsp_filter_synthesis(G,g,f,param);
    param.method = 'lanczos';
    f2 = gsp_filter_synthesis(G,g,f,param);


    if norm(f1(:)-f2(:))/norm(f1(:)) < 1e-10
        errors = 0;
       fprintf('FILTER: lanczos synthesis ok\n');
    else
        errors = 1;
        norm(f1(:)-f2(:))/norm(f1(:))
        warning('FILTER: Error in lanczos synthesis test')
    end

end

function errors = test_filter_lanczos_syn_multiple()

    Nf = 4;
    G = gsp_sensor(10);
    G = gsp_estimate_lmax(G);
    G = gsp_compute_fourier_basis(G);
    g = gsp_design_meyer(G, Nf);  
    f = rand(G.N*Nf,5);
    
    param.method = 'exact';
    param.order = 100;
    f1 = gsp_filter_synthesis(G,g,f,param);
    param.method = 'lanczos';
    f2 = gsp_filter_synthesis(G,g,f,param);


    if norm(f1(:)-f2(:))/norm(f1(:)) < 1e-10
        errors = 0;
       fprintf('FILTER: lanczos synthesis multiple test ok\n');
    else
        errors = 1;
        norm(f1(:)-f2(:))/norm(f1(:))
        warning('FILTER: Error in lanczos synthesis multiple test')
    end

end
