clc;
format longG
[elev_intensity,R] = readgeoraster("Site001PSR.tif");    % Elevation map
elev_intensity = double(elev_intensity); % Intensity(1) is y, Intensity(2) is x

y_tot = size(elev_intensity,1);
x_tot = size(elev_intensity,2);

[slope_intensity,Q] = readgeoraster("Site001PSRSlope.tif");
slope_intensity = double(slope_intensity);

% Load in the nodal map with weight of slope
load("Site001PSR_SlopeInside.mat", "s","t","weights","G")

% Coordinates of crater entrance
X01a_en = [844+1 (1+158)];
X01b_en = [712+1 (1+184)];
M01g_en = [1189+1 (1+89  )];
M01f_en = [780+1 (1+558)];
M01c_en = [436+1 (1+571)];
M01b_en = [118+1 (1+403)];
% M01b_en = [113+1 (1+497)];
entrance = [M01g_en; X01b_en; X01a_en; M01c_en; M01b_en; M01f_en];

% Coordinates of crater centers
X01a_c = [856+1 (1+144)];
X01b_c = [718+1 (1+184)];
M01g_c = [1230+1 (1+105)];
M01f_c = [735+1 (1+601)];
M01c_c = [406+1 (1+602)];
M01b_c = [140+1 (1+435)];
% M01b_c = [159+1 (1+452)];
center = [M01g_c; X01b_c; X01a_c; M01c_c; M01b_c; M01f_c];

% Coordinates of crater exits
X01a_ex = [871+1 (1+129)];
X01b_ex = [724+1 (1+183)];
M01g_ex = [1266+1 (1+106)];
M01f_ex = [683+1 (1+632)];
M01c_ex = [374+1 (1+639)];
M01b_ex = [202+1 (1+390)];
% M01b_ex = [102+1 (1+407)];
exit = [M01g_ex; X01b_ex; X01a_ex; M01c_ex; M01b_ex; M01f_ex];

% Plot the image
image = imread("Site001PSRSlope.tif");
figure(1), clf
imshow(image, [])
hold on
% text(Site(1),(y_tot - Site(2))+30,"LZ")

% Cycle through each of the locations from entrance to center to exit
for site = 1:6
    str = '';
    if site == 3
        str = "X01a";
    elseif site == 2
        str = "X01b";
    elseif site == 1
        str = "M01g";
    elseif site == 6
        str = "M01f";
    elseif site == 4
        str = "M01c";
    elseif site == 5
        str = "M01b";
    end
    fprintf("From Site to <strong>%s</strong>\n", str)
    for j = 1:2
        % Define starting position and goal
        if j == 1
            start = entrance(site,:);
            goal = center(site,:);
        elseif j == 2
            start = center(site,:);
            goal = exit(site,:);
        end 
        main(start,goal,elev_intensity,slope_intensity,G)
    end
end

function [path, total_cost] = A_star_path_planning_shortestpath(map, start_pixel, target_pixel, G)
     % Get dimensions of the elevation map
    [rows, cols] = size(map);
    
    % Compute x and y coordinates for each node
    [x, y] = meshgrid(1:cols, 1:rows);
    
    % Find the shortest path using built-in shortestpath function
    [path_nodes, total_cost] = shortestpath(G, start_pixel, target_pixel);
    
    % Convert node indices to coordinates
    path_x = zeros(1,numel(path_nodes));
    path_y = zeros(1,numel(path_nodes));
    path = cell(size(path_nodes));
    for i = 1:numel(path_nodes)
        [row, col] = ind2sub(size(map), path_nodes(i));
        path{i} = [col, rows - row];
        path_x(i) = col;
        path_y(i) = rows-row;
    end

    plot(path_x, rows - path_y,'LineWidth',2)
end

function main(start,goal,elev_intensity,slope_intensity,G)  
    fprintf("Coordinates start: [%i, %i]\n", start(1), start(2))
    fprintf("Coordinates end: [%i, %i]\n", goal(1), goal(2))
        
    % Convert start and stop pixel coordinates to nodes
    y_tot = size(elev_intensity,1);
    x_tot = size(elev_intensity,2);
    start_node = sub2ind([y_tot, x_tot], start(2), start(1));
    goal_node = sub2ind([y_tot, x_tot], goal(2), goal(1));
    
    [path] = A_star_path_planning_shortestpath(slope_intensity, start_node, goal_node, G);
    
    [max_slope,avg_slope,tot_dist] = max_of_path(path,elev_intensity,slope_intensity,y_tot);
    
    fprintf("   <strong>Max slope encountered</strong>: %f degrees\n", max_slope)
    fprintf("   Total Distance driven: %f meters\n", tot_dist)
    fprintf("   Average slope: %f degrees\n", avg_slope)
end

function [max_slope,avg_slope,tot_dist] = max_of_path(path, elev_intensity, slope_intensity,y_tot)
    % Calculate total distance, max slope, and average slope
    tot_dist = 0.0;
    max_slope = 0.0;
    avg_slope = 0.0;
    for i = 1:(size(path,2) - 1)
        curr_xy = path{i};
        next_xy = path{i + 1};
%         curr_node = (1 + y_tot*(curr_xy(1)-1) + (curr_xy(2) - 1))
%         next_node = (1 + y_tot*(next_xy(1)-1) + (next_xy(2) - 1))
        del_elev = abs(elev_intensity(curr_xy(2), curr_xy(1)) ...
            - elev_intensity(next_xy(2), next_xy(1)));
        next_slope = slope_intensity(y_tot - next_xy(2),next_xy(1));
        % Determine if the movement is diagonal or straight
        dist = 0.0;
        diag = (next_xy(1) - curr_xy(1))^2 + (next_xy(2) - curr_xy(2))^2;
        if diag == 2
            dist = sqrt((7.07^2) + (del_elev)^2);
            tot_dist = tot_dist + dist;
        elseif diag == 1
            dist = sqrt((5^2) + (del_elev)^2);
            tot_dist = tot_dist + dist;
        end
        % Check if slope is a new max
        avg_slope = avg_slope + next_slope;
        if next_slope > max_slope
            max_slope = next_slope;
%             disp(next_xy(2))
%             disp(next_xy(1))
        end
    end
    avg_slope = avg_slope / size(path,2);
end