figure();
f = @(x,y,z) x+z + sqrt((x-z).^2 + y.^2) - 10;
g = @(x,y,z) x+z - sqrt((x-z).^2 + y.^2) - 10;
h = @(x,y,z) x+z + sqrt((x-z).^2 + y.^2) - 100;
i = @(x,y,z) x+z - sqrt((x-z).^2 + y.^2) - 100;
fimplicit3(f, [0, 60], 'EdgeColor', 'none', 'FaceAlpha', .5, "FaceColor", 'b')
hold on;
fimplicit3(g, [0, 60], 'EdgeColor', 'none', 'FaceAlpha', .5, "FaceColor", 'b')
fimplicit3(h, [0, 60], 'EdgeColor', 'none', 'FaceAlpha', .5, "FaceColor", 'r')
fimplicit3(i, [0, 60], 'EdgeColor', 'none', 'FaceAlpha', .5, "FaceColor", 'r')
% , 'EdgeColor', 'none', 'FaceAlpha', .5
title('$2 \times 10  \leq (a + c) \pm \sqrt{(a-c)^2 + b^2} \leq 2 \times 100$', 'Interpreter', 'latex')
xlabel('a')
ylabel('b')
zlabel('c')
print(gcf,'Cones_3.png','-dpng','-r600'); 
