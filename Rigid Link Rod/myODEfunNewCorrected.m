function dYdt = myODEfunNewCorrected(t,Y)
%dYdt = myODEfun(t,Y,M,m,n,g,k,l)
%n_size = n-1;               % No units
n_size=length(Y)/2;
n=n_size+1;
k = 100;                   % Nm/radian
M = 0.1;                   % kg  
mbeam=0.1;                 % kg
m=mbeam/n;                 % kg
g = 9.81;                  % m/s^2
l = 1 / n;                % Metres

dYdt = zeros((2 * (n_size)), 1);
%l = 1 / n_size+1;  
A = zeros(n_size, n_size);
for i = 1:(n_size)
    for j = 1:(n_size)
        if j<=i
            A(i, j) = ((M*l^2)+(n-i-(1/2))*m*l^2)*cos(Y(i) - Y(j));
        elseif j>i
            A(i, j) = ((M*l^2)+(n-j-(1/2))*m*l^2)*cos(Y(i) - Y(j));
        end
    end
end

%A_ = zeros(n_size, n_size);
A_= A-((m*l^2)/6)*eye(n_size);
t
DetA_=det(A_)

B = zeros(n_size, n_size);
for i = 1:(n_size)
    for j = 1:(n_size)            
        if j<=i
            B(i, j) = ((M*l^2)+(n-i-(1/2))*m*l^2)*sin(Y(i) - Y(j));
        elseif j>i
            B(i, j) = ((M*l^2)+(n-j-(1/2))*m*l^2)*sin(Y(i) - Y(j));
        end
    end
end


C=zeros(n_size,1);
for i=1:n_size
    (M+m*(n-(i+1)+(1/2)))*g*l*sin(Y(i));
end

D=zeros(n_size, n_size);
for i=1:n_size
    for j=1:n_size
        if i==1
            D(i,1)=1;
            %D(i,2)=-1;
  
        elseif i==n_size
            D(i,n_size-1)=-1;
            D(i,n_size)=1;
            
        else
            if j==i
                D(i,j)=2;
            elseif j==((i-1)||(i+1))
                D(i,j)=-1;
            end
        end
    end
end
% for i=n_size
%     for j=1:n_size
%         if j==n_size-1
%             D(i,j)=-1;
%         elseif j==n_size
%             D(i,j)=1;
%         else
%             D(i,j)=0;
%         end
%     end
% end

% Alphavec Alphasquarevec
alphavec=zeros(n_size,1);
alphadsqvec=zeros(n_size,1);
alphadvec=zeros(n_size,1);
for i=1:n_size
    alphavec(i)=Y(i);
    alphadvec(i)=Y(n_size+i);
    alphadsqvec(i)=(Y(n_size+i)^2);
end

b=zeros(n_size, n_size);
b=-(B*alphadsqvec)-k*D*alphavec+C;

for i = 1:n_size
    dYdt(i) = Y(n_size + i);
end

dYdt((n_size + 1):end) = A_\b; % contains 2 to n alpha_dots and alpha_doubledots
end
