function error = test_resistance_distance()
    error = 0;
    
    W = rand(10);
    W = W-diag(diag(W));
    W = (W + W')/2;
    
    G = gsp_graph(W);
    
    % Method 1
    rd1 = gsp_resistance_distance(G);
    
    % Method 2
    pseudo=pinv(full(G.L));
    rd2=zeros(size(G.L));
    N=size(G.L,1);
    for i=1:N
        for j=1:N
            rd2(i,j)=pseudo(i,i)+pseudo(j,j)-2*pseudo(i,j);
        end
    end
        
    
    if norm(rd1-rd2) < 1e-10
         fprintf('RESISTANCE DISTANCE: test 1 ok\n');
    else
        error = error + 1;
        warning('RESISTANCE DISANCE: error in test 1');
    end
    
    
    % The weigh 
    dist = [0, 3, 1;...
            3, 0, 2;...
            1, 2, 0];
    % The weight is the inverse of the distance...
    W = dist.^(-1);
    % Fix the diagonal
    W([1,5,9])=0;
    
    G = gsp_graph(W);

    
    % Change laplacian    
    rd3 = gsp_resistance_distance(G);

    
    % Method 4 Hand
    
    rd4 = [0, 3/2, 5/6;...
           3/2, 0, 4/3;...
           5/6, 4/3, 0];
       
    if norm(rd3-rd4) < 1e-10
         fprintf('RESISTANCE DISTANCE: test 2 ok\n');
    else
        error = error + 1;
        warning('RESISTANCE DISANCE: error in test 2');
    end
    
    
    G = gsp_sensor(20);
    rd3 = gsp_resistance_distance(G);
    rd4 = compute_resistance_distances(G.L);
    if sum(abs(rd3(:)-rd4(:))) < 1e-10
         fprintf('RESISTANCE DISTANCE: test 3 ok\n');
    else
        error = error + 1;
        warning('RESISTANCE DISANCE: error in test 3');
    end
    
    rd5 = gsp_compute_resistance_distances_old(G.L);
    
    if sum(abs(rd3(:)-rd5(:))) < 1e-10
         fprintf('RESISTANCE DISTANCE: test 4 ok\n');
    else
        error = error + 1;
        warning('RESISTANCE DISANCE: error in test 4');
    end
    
end
