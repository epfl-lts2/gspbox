function errors = test_gsp_remove_mean()


X = rand(100,10);

X1 = gsp_remove_mean(X);

errors = gsp_assert_test(0,sum(X1),eps(1000), 'remove mean dim 1');

X2 = gsp_remove_mean(X,2);

errors = errors + gsp_assert_test(0,sum(X2,2),eps(1000), 'remove mean dim 2');

X = rand(100,10,3);

X3 = gsp_remove_mean(X,3);

errors = errors + gsp_assert_test(0,sum(X3,3),eps(1000), 'remove mean dim 3');

end