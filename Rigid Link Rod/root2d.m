function F = root2d(x)

F(1) = k * (x(1) - x(2)) - g * l * sin(x(1)) * (M + m * (n_links - 0.5));

if (n_links > 2)
    for i_1 = 2:(n_links-1)
        F(i_1) = k * (-x(i_1-1) + 2 * x(i_1) - x(i_1+1)) - g * l * sin(x(i_1)) * (M + m * (n_links - i_1 + 0.5));
    end
end

F(n_links) = k * (x(n_links) - x(n_links - 1)) - g * l * sin(x(n_links)) * (M + m * 0.5);

% F(1) = k * (abs(x(1)) - abs(x(2))) - g * l * sin(abs(x(1))) * (M + m * (n_links - 0.5));
% 
% for i_1 = 2:(n_links-1)
%     F(i_1) = k * (-abs(x(i_1-1)) + 2 * abs(x(i_1)) - abs(x(i_1+1))) - g * l * sin(abs(x(i_1))) * (M + m * (n_links - i_1 + 0.5));
% end
% 
% F(n_links) = k * (abs(x(n_links)) - abs(x(n_links - 1))) - g * l * sin(abs(x(n_links))) * (M + m * 0.5);

end