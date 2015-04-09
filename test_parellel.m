V = zeros(10000,1000);

% matlabpool(4)

tic
parfor i = 1:10000
    for j = 1:1000
        V(i,j) = randn(1);
    end
end
toc

% matlabpool close