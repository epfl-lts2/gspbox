function errors = test_gsp_hope_distanz()
    
errors = 0;

G = gsp_path(10);

d = gsp_hop_distanz(G,1,5);

if d==4
    fprintf(' GSP_HOP_DISTANZ 1 OK\n')
else
    warning(' GSP_HOP_DISTANZ 1 pas OK')
    errors = errors + 1;
end

d = gsp_hop_distanz(G,1,8);

if d==7
    fprintf(' GSP_HOP_DISTANZ 2 OK\n')
else
    warning(' GSP_HOP_DISTANZ 2 pas OK')
    errors = errors + 1;
end


end