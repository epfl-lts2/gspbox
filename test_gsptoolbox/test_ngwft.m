function [ error ] = test_ngwft()
%TEST_FAST_GWFT test the fast grap windowed fourier transform

error=0;
tol=10e-14;

error = error + test_ring(tol);
error = error + test_frame(tol);
error = error + test_ring2(tol);
error = error + test_frame2(tol);

end


function error=test_ring(tol)

    N=16;

    % Create the graph
    G=gsp_ring(N);
    G=gsp_compute_fourier_basis(G);

    % Create the window 
    fhat = rand(N,1);
    fhat = fhat/norm(fhat);

    f = gsp_igft(G,fhat);

    % Compute the ambiguity function
    Af = gsp_gwft(G,f,f);
    Afn = gsp_ngwft(G,f,f);
    
    if norm(Af(:)-Afn(:))>tol
        error = 1;
        fprintf(' - Test NWGFT ring: ERROR\n'); 

    else
        error = 0;
        fprintf(' - Test NWGFT ring: OK\n'); 

    end


   
end


function error=test_frame(tol)

    N=16;

    % Create the graph
    G=gsp_sensor(N);
    G=gsp_compute_fourier_basis(G);

    % Create the window 
    fhat = rand(N,1);
    fhat = fhat/norm(fhat);

    f = gsp_igft(G,fhat);

    % Compute the ambiguity function
    Af = gsp_ngwft(G,f,f);
    param.lowmemory=0;
    Afn = gsp_ngwft(G,f,f,param);
    
    if norm(Af(:)-Afn(:))>tol
        error = 1;
        fprintf(' - Test NWGFT frame: ERROR\n'); 
        norm(Af(:)-Afn(:))

    else
        error = 0;
        fprintf(' - Test NWGFT frame: OK\n'); 

    end


   
end

function error=test_ring2(tol)

    N=16;

    % Create the graph
    G=gsp_ring(N);
    G = gsp_create_laplacian( G,'normalized' );
    G=gsp_compute_fourier_basis(G);

    % Create the window 
    fhat = rand(N,1);
    fhat = fhat/norm(fhat);

    f = gsp_igft(G,fhat);

    % Compute the ambiguity function
    Af = gsp_gwft(G,f,f);
    Afn = gsp_ngwft(G,f,f);
    
    if norm(Af(:)-Afn(:))>tol
        error = 1;
        fprintf(' - Test NWGFT ring: ERROR\n'); 

    else
        error = 0;
        fprintf(' - Test NWGFT ring: OK\n'); 

    end


   
end


function error=test_frame2(tol)

    N=16;

    % Create the graph
    G=gsp_sensor(N);
    G = gsp_create_laplacian( G,'normalized' );
    G=gsp_compute_fourier_basis(G);

    % Create the window 
    fhat = rand(N,1);
    fhat = fhat/norm(fhat);

    f = gsp_igft(G,fhat);

    % Compute the ambiguity function
    Af = gsp_ngwft(G,f,f);
    param.lowmemory=0;
    Afn = gsp_ngwft(G,f,f,param);
    
    if norm(Af(:)-Afn(:))>tol
        error = 1;
        fprintf(' - Test NWGFT frame: ERROR\n'); 
        norm(Af(:)-Afn(:))

    else
        error = 0;
        fprintf(' - Test NWGFT frame: OK\n'); 

    end


   
end

