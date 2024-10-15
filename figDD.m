% Step (i): Construct a coarse triangulation TH
figure;
subplot(1, 2, 1);
hold on;

% Define the coarse grid (TH)
for i = 0:5
    for j = 0:5
        % Shade specific squares
        if (i == 1 && j == 2) || (i == 3 && j == 1)
            fill([i, i+1, i+1, i], [j, j, j+1, j+1], [0.6, 0.6, 0.6]);
        end
        plot([i, i+1], [j, j], 'k'); % Horizontal line
        plot([i, i], [j, j+1], 'k'); % Vertical line
        plot([i, i+1], [j+1, j], 'k'); % Diagonal line

        
    end
end
axis equal;
axis off;
axis([0,4,0,4])
title('Coarse meshes');
hold off;

% Step (ii): Define subdomains aligned with TH
% Subdomains are shaded regions as shown

% Step (iii): Subdivide TH to obtain Th
subplot(1, 2, 2);
hold on;
fill([0.5, 2.5, 2.5, 0.5], [1.5, 1.5, 3.5, 3.5], [0.6, 0.6, 0.6], 'EdgeColor', 'none');
fill([2.5, 4, 4, 2.5], [0.5, 0.5, 2.5, 2.5], [0.6, 0.6, 0.6], 'EdgeColor', 'none');

% Plot the same coarse grid
for i = 0:5
    for j = 0:5
        plot([i, i+1], [j, j], 'k'); % Horizontal line
        plot([i, i], [j, j+1], 'k'); % Vertical line
        plot([i, i+1], [j+1, j], 'k'); % Diagonal line
    end
end


% Step (iv): Refine specific regions (Th)
% Define the refined grids
i = 1;
j = 2;
for k = -0.5:0.25:1.25
    for l = 2.5:0.25:4.25
        plot([k, k+0.25] + i, [l, l] - 3 + j, 'b'); % Horizontal line
        plot([k, k] + i, [l, l+0.25] - 3 + j, 'b'); % Vertical line
        plot([k, k+0.25] + i, [l+0.25, l] - 3 + j, 'b'); % Diagonal line
    end
end
plot([2.5,2.5],[0.5,3.5],'b')
plot([0.5,2.5],[3.5,3.5],'b')

i = 3;  
j = 1;
for k = -0.5:0.25:1.25
    for l = 0.5:0.25:2.25
        plot([k, k+0.25] + i, [l, l] - 1 + j, 'b'); % Horizontal line
        plot([k, k] + i, [l, l+0.25] - 1 + j, 'b'); % Vertical line
        plot([k, k+0.25] + i, [l+0.25, l] - 1 + j, 'b'); % Diagonal line
    end
end
plot([2.5,4],[2.5,2.5],'b')
% plot([2.5,4],[0.5,0.5],'b')
% Step (v): Highlight enlarged subdomains

axis equal;
axis off;
axis([0,4,0,4])
title('Fine meshes in two subdomains');
hold off;

set(gcf, 'Color', 'w');
export_fig ('figDD.eps')