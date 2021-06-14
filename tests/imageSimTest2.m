
idx = linspace(0, 10*pi, 1001);
idx = idx(1:end-1);

D = 100;

a = sign(sin(idx)) + 1;
a = a/max(a(:));

a(a==0) = .3*exp(2*pi*i * 0.5);
NA = 3.01*6.7500e-04;

plot(a)

plot(abs(imageSim(a, D, .0135, NA)))


%%


idx = linspace(0, 10*pi, 201);
idx = idx(1:end-1);

D = [100, 100];

a = sign(sin(idx)) + 1;


A = ones(size(idx))'* a;

a = a/max(a(:));

NA = 3.01*6.7500e-04;


imagesc(abs(imageSim(A, D, .0135, NA)))
