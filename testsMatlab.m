lambda = 10
N = 10
seed = 1
rand('seed',seed)
for k=9:lambda,
  display(k)
  xmean = rand(N,1); % objective variables initial point
  display(xmean)
end  