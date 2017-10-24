function [ errors ] = test_plotting(  )
errors = 0;

errors = errors +test_signal_spectral();
errors = errors +test_signal();
errors = errors +test_directed();

errors = errors +test_sgram();

end


function errors = test_signal_spectral()
    errors = 0;

    try
        figure(100)
          N = 32;
          G = gsp_path(N);
          G = gsp_compute_fourier_basis(G);
          f = sin((1:N)'*2*pi*4/N);
          fhat = gsp_gft(G,f);
          gsp_plot_signal_spectral(G,fhat);
          close(100)
        fprintf('PLOTTING: spectral ok\n');
    catch
            errors = errors + 1;
            warning('PLOTTING: Error in the spectral test')
    end

end


function errors = test_signal()
    errors = 0;

    try
        figure(100)
        G = gsp_ring(15);
        f = sin((1:15)'*2*pi/15);
        paramplot.show_edges = 1;
        gsp_plot_signal(G,f,paramplot)
        close(100)
        fprintf('PLOTTING: vertex 2d ok\n');
    catch
        errors = errors + 1;
        warning('PLOTTING: Error in the vertex 2d test')
    end
    
    try
        figure(100)
        N = 10;
        W = rand(N);
        W = W+W';
        W = W-diag(diag(W));
        G = gsp_graph(W,rand(N,3));
        f = rand(N,1);
        gsp_plot_signal(G,f,paramplot)
        close(100)
        fprintf('PLOTTING: vertex 3d ok\n');
    catch
        errors = errors + 1;
        warning('PLOTTING: Error in the vertex 3d test')
    end
    
     try
        figure(100)
        G = gsp_ring(15);
        f = sin((1:15)'*2*pi/15);
        paramplot.bar = 1;
        gsp_plot_signal(G,f,paramplot)
        close(100)
        fprintf('PLOTTING: vertex 2d bar ok\n');
    catch
        errors = errors + 1;
        warning('PLOTTING: Error in the vertex 2d bar test')
    end

end

function errors = test_directed()
    errors = 0;

    try
        figure(100)
        G = gsp_ring(15);
        f = sin((1:15)'*2*pi/15);
        paramplot.show_edges = 1;
        G.directed = 1;
        gsp_plot_signal(G,f,paramplot)
        close(100)
        fprintf('PLOTTING: directed 2d ok\n');
    catch
        errors = errors + 1;
        warning('PLOTTING: Error in the directed 2d test')
    end
    
    try
        figure(100)
        N = 10;
        W = rand(N);
        W = W-diag(diag(W));
        G = gsp_graph(W,rand(N,3));
        G.directed = 1;
        f = rand(N,1);
        gsp_plot_signal(G,f,paramplot)
        close(100)
        fprintf('PLOTTING: directed 3d ok\n');
    catch
        errors = errors + 1;
        warning('PLOTTING: Error in the directed 3d test')
    end
    
     try
        figure(100)
        G = gsp_ring(15);
        G.directed = 1;
        f = sin((1:15)'*2*pi/15);
        paramplot.bar = 1;
        gsp_plot_signal(G,f,paramplot)
        close(100)
        fprintf('PLOTTING: directed 2d bar ok\n');
    catch
        errors = errors + 1;
        warning('PLOTTING: Error in the directed 2d bar test')
    end

end


function errors = test_sgram()
    errors = 0;

    try
       figure(100)
      N = 15;
      G = gsp_ring(2*N);
      G = gsp_compute_fourier_basis(G);
      x = [0:N,(N-1):-1:1]';
      s = 3;
      g = exp(-x.^2/s^2);
      f = gsp_modulate(G,gsp_translate(G,g,N),N);
      c = gsp_gwft(G,f,g);
      gsp_plot_sgram(G,c);
      close(100)
      fprintf('PLOTTING: sgram ok\n');
    catch
        errors = errors + 1;
        warning('PLOTTING: Error in the sgram test')
    end

end
