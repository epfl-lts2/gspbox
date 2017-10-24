function [ error ] = test_gwft()
%TEST_GWFT test the fast grap windowed fourier transform

error=0;
tol=10e-8;

error = error+test1(tol);
error = error+test2(tol);

error = error+test3(tol);

end


function error=test1(tol)

    error = 0;
    N=16;

    % create a signal
    f=rand(N,1);
    % create a random graph
    G = gsp_random_regular(N,5);
    G=gsp_compute_fourier_basis(G);
    % create a random window;
    g=rand(N,1);

    c1=WGFT_wgft(g,G.U,f);

    c2=gsp_gwft( G,g,f);

    param.lowmemory=0;
    c3=gsp_gwft(G,g,f,param);

    % 
    % figure
    % subplot(121)
    % imagesc(c1)
    % subplot(122)
    % imagesc(c3)

    if norm(c1-c2)>tol
        error=error+1;
    end


    if norm(c1-c3)>tol
        error=error+1;
    end

    if error
       fprintf(' - Test WGFT: ERROR\n'); 
    else

       fprintf(' - Test WGFT: OK\n');
    end
end

function error=test2(tol)

    error = 0;
    N=16;

    % create a signal
    f=rand(N,1);
    % create a random graph
    G = gsp_random_regular(N,5);
    G=gsp_compute_fourier_basis(G);
    G = gsp_create_laplacian( G,'normalized' );
    param.verbose = 0;
    G=gsp_compute_fourier_basis(G,param);
    % create a random window;
    g=rand(N,1);

    c1=WGFT_wgft(g,G.U,f);

    c2=gsp_gwft(G,g,f);

    param.lowmemory=0;
    c3=gsp_gwft(G,g,f,param);

    % 
    % figure
    % subplot(121)
    % imagesc(c1)
    % subplot(122)
    % imagesc(c3)

    if norm(c1-c2)>tol
        error=error+1;
    end


    if norm(c1-c3)>tol
        error=error+1;
    end

    if error
       fprintf(' - Test WGFT Normalized Laplacian: ERROR\n'); 
    else

       fprintf(' - Test WGFT Normalized Laplacian: OK\n');
    end
end



function error=test3(tol)

    error = 0;
    N=6;

    % create a signal
    f=rand(N,1);
    % create a random graph
    G = gsp_random_regular(N,5);
    G=gsp_compute_fourier_basis(G);
    % create a random window;
    g=rand(N,1);
    f = rand(N,2);
    
    param.lowmemory=1;
    
    c1 = zeros(N,N,2);
    c1(:,:,1) = gsp_gwft( G,g,f(:,1),param);
    c1(:,:,2) = gsp_gwft( G,g,f(:,2),param);
    
    c2 = gsp_gwft( G,g,f,param);

    param.lowmemory=0;
    c3 = gsp_gwft(G,g,f,param);

    % 
    % figure
    % subplot(121)
    % imagesc(c1)
    % subplot(122)
    % imagesc(c3)

    if norm(c1(:)-c2(:))>tol
        error=error+1;
    end


    if norm(c1(:)-c3(:))>tol
        error=error+1;
    end

    if error
       fprintf(' - Test WGFT: ERROR\n'); 
    else

       fprintf(' - Test WGFT: OK\n');
    end
end
