lambda_list = [.1 1 10 100];
lambda_list = [.01 .1 1];
for ilambda = 1:length(lambda_list)
    lambda = lambda_list(ilambda);
    test_radon;
end