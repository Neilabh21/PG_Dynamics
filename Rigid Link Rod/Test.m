function x_all = main()
%Some code
d1 = 1;
d3 = 1;
k1 = 2;
k3 = 3;
g = 9.81;
w3 = 5;
k1in= 5;
k3in= 4;
alpha1in=10;
alpha3in=10;

x0 = zeros(1,4);

% Values for w1
w1_pool = 1:1:10;

% Preallocate result matrix
x_all = zeros(numel(w1_pool),numel(x0))

% switch dispay off
options = optimoptions('fsolve','Display','off');

% Call fsolve in a loop
for k = 1 : numel(w1_pool)
    w1 = w1_pool(k);
    x_all(k,:) = fsolve(@obj_fun,x0,options);
end

function F =  obj_fun(x)
F(1:4) = [x(1) - 2; x(2) - 5; x(3) - 3; x(4) - 7];
% F(2)=x(2)*(-1i*(w1-d1)+k1)-conj(g)*x(4)*x(1)-sqrt(2*k1in)*conj(alpha1in);
% F(3)=x(3)*(1i*(w3-d3)+k3)+(g/2)*x(1)^2-sqrt(2*k3in)*alpha3in;
% F(4)=x(4)*(1i*(w3-d3)+k3)+(g/2)*x(2)^2-sqrt(2*k3in)*alpha3in;
end
end