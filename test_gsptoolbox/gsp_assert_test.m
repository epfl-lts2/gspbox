function e = gsp_assert_test(X1,X2,tol, name)

if X1==0
    if norm(X2(:),'fro')<tol
        fprintf(['Test ',name, ': OK\n']);
        e = 0;
    else
        warning(['Test ',name, ': Error, badness ',...
            num2str(norm(X2(:),'fro')),...
            ', tol ',num2str(tol)]);
        e = 1;
    end
else
    if norm(X1(:)-X2(:),'fro')/norm(X2(:),'fro')<tol
        fprintf(['Test ',name, ': OK\n']);
        e = 0;
    else
        warning(['Test ',name, ': Error, badness ',...
            num2str(norm(X1(:)-X2(:),'fro')/norm(X2(:),'fro')), ...
            ', tol ',num2str(tol)]);
        e = 1;
    end
end

end