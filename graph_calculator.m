% Input slope and elevation map
[intensity,R] = readgeoraster("Site001PSR.tif");
elevation_map = double(intensity);

% Get dimensions of the elevation map
[rows, cols] = size(elevation_map);
    
% Initialize lists to store source, target, and weights of edges
W_el = 0; % Elevation weight
W_d = 1;  % Distance weight

s = []; % Start node
t = []; % End node
weights = [];   % Edge weight
    
% Define the 8-connected neighbors for each pixel
offsets = [-1, -1; -1, 0; -1, 1; 0, -1; 0, 1; 1, -1; 1, 0; 1, 1];
    
% Iterate over each pixel in the elevation map
for i = 1:rows
    for j = 1:cols
        % Compute index of current pixel
        current_index = sub2ind([rows, cols], i, j);
            
        % Iterate over each neighbor of the current pixel
        for k = 1:size(offsets, 1)
             % Compute indices of neighbor pixel
             neighbor_row = i + offsets(k, 1);
             neighbor_col = j + offsets(k, 2);
                
             % Check if neighbor is within bounds
             if neighbor_row >= 1 && neighbor_row <= rows && ...
                neighbor_col >= 1 && neighbor_col <= cols
                % Compute index of neighbor pixel
                neighbor_index = sub2ind([rows, cols], neighbor_row, neighbor_col);
                
                % Compute weight (elevation difference) between current pixel and neighbor
                elev = abs(elevation_map(i, j) - elevation_map(neighbor_row, neighbor_col));
                dist = 5.0;
                % If neighbor is diagonal the path is longer
                if k == 1 || k == 3 || k == 6 || k == 8
                    dist = sqrt((7.07^2) + (elev)^2);
                else
                    dist = sqrt((5^2) + (elev)^2);
                end
                slope = asind(elev/dist);
                % Normalize the distance and elevation based on max 20 deg
                max_slope = 20;
                max_dist = sqrt((7.07^2) + (7.07*sind(max_slope))^2);
                max_elev = 7.07*sind(max_slope);
                dist = dist / max_dist;
                elev = elev / max_elev;
                % Ignore options that 
                if slope >= max_slope
                    weight = inf;
                else 
                    weight = W_el*elev + W_d*dist; 
                end
                % Add edge from current pixel to neighbor with weight
                s(end+1) = current_index;
                t(end+1) = neighbor_index;
                weights(end+1) = weight;
            end
        end
    end
end

% Create graph object
G = graph(s, t, weights);
% Store the graph, nodes, and edges
fileName = "Site001PSR_inside.mat";
save(fileName,"s","t","weights","G")
