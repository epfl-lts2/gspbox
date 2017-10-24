function d = gsp_diameter(G)
    
d = 0;
v = bfs(G.A,1);
if any(v<0)
    warning('Disconnected graph!')
end
for ii = 1:G.N
    d = max(d,max(bfs(G.A,ii)));
end

end