using CUDA

A = CUDA.rand(100,100,100,100,100);
AA = findmax(A; dims=1);
AA1 = dropdims(AA[1]; dims=1);
AA2 = getindex.(dropdims(AA[2]; dims=1), 1);