clc;
format longG
[elev_intensity,R] = readgeoraster("Site001PSR.tif");    % Elevation map
elev_intensity = double(elev_intensity); % Intensity(1) is y, Intensity(2) is x

y_tot = size(elev_intensity,1);
x_tot = size(elev_intensity,2);

[slope_intensity,Q] = readgeoraster("Site001PSRSlope.tif");
slope_intensity = double(slope_intensity);

% Plot the image
image = imread("Site001PSRSlope.tif");
figure(1), clf
imshow(image, [])
hold on
% text(Site(1),(y_tot - Site(2))+30,"LZ")

tot_dist_m01g = [];
max_slope_m01g = [];
avg_slope_m01g = [];
tot_dist_x01b = [];
max_slope_x01b = [];
avg_slope_x01b = [];
tot_dist_x01a = [];
max_slope_x01a = [];
avg_slope_x01a = [];
tot_dist_m01b = [];
max_slope_m01b = [];
avg_slope_m01b = [];
tot_dist_m01c = [];
max_slope_m01c = [];
avg_slope_m01c = [];
tot_dist_m01f = [];
max_slope_m01f = [];
avg_slope_m01f = [];

% figure(2),clf
% hold on
% xlabel("Total Distance")
% ylabel("Max Slope")
% for W_sl = 0.9:0.02:1
    
    % Load in the nodal map with weight of slope
%     filename = strcat('Site001PSR_Slope',num2str(W_sl), '.mat');
%     load(filename, "s","t","weights","G")
    load("Site001PSR_Slope.mat","s","t","weights","G")
    
    % Different coordinates for different PSRs
    Site = [742+1 (1+299)];
    X01a = [862+1 (1+167)];
    X01b = [724+1 (1+196)];
    M01g = [1219+1 (1+143)];
    M01f = [742+1 (1+516)];
    M01c = [436+1 (1+571)];
    M01b = [231+1 (1+439)];
    locations = [M01g; X01b; X01a; M01c; M01b; M01f];
    
    % Stop along the way
    M01b_stop1 = [591 (y_tot - 373)];
    M01b_stop2 = [425 (y_tot - 423)];
    stop_plot = 5;  % which index in locations stop is associated with
    stops = [M01b_stop1;M01b_stop2];
    stops = [];
    j = 1;
    
    % Plot the image
%     image = imread("Site001PSRSlope.tif");
%     figure(1), clf
%     imshow(image, [])
%     hold on
%     text(Site(1),(y_tot - Site(2))+30,"LZ")
    
    % Cycle through each of the locations from site to PSR
    for site = 1:size(locations,1)
        site;
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
        % Define starting position and goal
        if site == stop_plot && size(stops,1) ~= 0
            for i = 1:size(stops)+1
               if i == 1
                   start = Site;
                   goal = stops(i,:);
                   [max_slope, tot_dist, avg_slope] = main(start,goal,str,elev_intensity,slope_intensity,G);
               elseif i ~= size(stops)+1
                   start = stops(i-1,:);
                   goal = stops(i,:);
                   [max_slope, tot_dist, avg_slope] = main(start,goal,str,elev_intensity,slope_intensity,G);
               else 
                   start = stops(end,:);
                   goal = locations(site,:);
                   [max_slope, tot_dist, avg_slope] = main(start,goal,str,elev_intensity,slope_intensity,G);
               end
            end
        else
            start = Site;
            goal = locations(site,:);
            [max_slope, tot_dist, avg_slope] = main(start,goal,str,elev_intensity,slope_intensity,G);
            if site == 3
                str = "X01a";
                tot_dist_x01a(end+1) = tot_dist;
                max_slope_x01a(end+1) = max_slope;
                avg_slope_x01a(end+1) = avg_slope;
            elseif site == 2
                str = "X01b";
                tot_dist_x01b(end+1) = tot_dist;
                max_slope_x01b(end+1) = max_slope;
                avg_slope_x01b(end+1) = avg_slope;
            elseif site == 1
                str = "M01g";
                tot_dist_m01g(end+1) = tot_dist;
                max_slope_m01g(end+1) = max_slope;
                avg_slope_m01g(end+1) = avg_slope;
            elseif site == 6
                str = "M01f";
                tot_dist_m01f(end+1) = tot_dist;
                max_slope_m01f(end+1) = max_slope;
                avg_slope_m01f(end+1) = avg_slope;
            elseif site == 4
                str = "M01c";
                tot_dist_m01c(end+1) = tot_dist;
                max_slope_m01c(end+1) = max_slope;
                avg_slope_m01c(end+1) = avg_slope;
            elseif site == 5
                str = "M01b";
                tot_dist_m01b(end+1) = tot_dist;
                max_slope_m01b(end+1) = max_slope;
                avg_slope_m01b(end+1) = avg_slope;
            end
        end
    end
% end
data = [tot_dist_m01b;tot_dist_m01g;tot_dist_m01f;tot_dist_m01c;tot_dist_x01b;tot_dist_x01a];
max_slope_data = [max_slope_m01b;max_slope_m01g;max_slope_m01f;max_slope_m01c;max_slope_x01b;max_slope_x01a];
avg_slope_data = [avg_slope_m01b;avg_slope_m01g;avg_slope_m01f;avg_slope_m01c;avg_slope_x01b;avg_slope_x01a];
titles = {'Path to M01b', 'Path to M01g', 'Path to M01f', 'Path to M01c', 'Path to X01b', 'Path to X01a'};
% legend_entries = ['Wsl=0,Wd=1';'Wsl=0.1,Wd=0.9';'Wsl=0.2,Wd=0.8';'Wsl=0.3,Wd=0.7';'Wsl=0.4,Wd=0.6'; ...
% 'Wsl=0.5,Wd=0.5';'Wsl=0.6,Wd=0.4';'Wsl=0.7,Wd=0.3';'Wsl=0.8,Wd=0.2';'Wsl=0.9,Wd=0.1';'Wsl=1,Wd=0'];
legend_entries = {'Wsl=0.90,Wd=0.10', 'Wsl=0.92,Wd=0.08', 'Wsl=0.94,Wd=0.06', 'Wsl=0.96,Wd=0.04', 'Wsl=0.98,Wd=0.02', ...
'Wsl=1.0,Wd=0.0'};
% legend_entries = {'Wsl=0,Wd=1', 'Wsl=0.1,Wd=0.9', 'Wsl=0.2,Wd=0.8', 'Wsl=0.3,Wd=0.7', 'Wsl=0.4,Wd=0.6', ...
% 'Wsl=0.5,Wd=0.5', 'Wsl=0.6,Wd=0.4', 'Wsl=0.7,Wd=0.3', 'Wsl=0.8,Wd=0.2', 'Wsl=0.9,Wd=0.1', 'Wsl=1,Wd=0'};

% Reconsctruction of the path
% for i = 1:6
%     figure(i),clf
%     filename = strcat("WeightComp",num2str(i),'.png');
%     t = tiledlayout(1,2,"TileSpacing","compact");
%     nexttile
%     hold on
%     for j = 1:size(data,2)
%         plot(data(i,j),max_slope_data(i,j),'o','LineWidth',2)
%     end
%     plot(data(i,:),max_slope_data(i,:),'LineWidth',2)
%     legend(legend_entries,'FontSize',13)
%     title(titles{i},'FontSize',13)
%     xlabel("Distance (meters)",'FontSize',13)
%     ylabel("Max Slope (degs)",'FontSize',13)
%     nexttile
%     hold on
%     for j = 1:size(data,2)
%         plot(data(i,j),avg_slope_data(i,j),'o','LineWidth',2)
%     end
%     plot(data(i,:),avg_slope_data(i,:),'LineWidth',2)
%     legend(legend_entries,'FontSize',13)
%     title(titles{i},'FontSize',13)
%     xlabel("Distance (meters)",'FontSize',13)
%     ylabel("Average Slope (degs)",'FontSize',13)
%     hold off
%     exportgraphics(t,filename,'BackgroundColor','none','ContentType','vector')
%     k = k + 1;
% end

function [path, total_cost] = A_star_path_planning_shortestpath(map, start_pixel, target_pixel, G, str)
     % Get dimensions of the elevation map
    [rows, cols] = size(map);
    
    % Compute x and y coordinates for each node
    [x, y] = meshgrid(1:cols, 1:rows);
    
    % Find the shortest path using built-in shortestpath function
    [path_nodes, total_cost] = shortestpath(G, start_pixel, target_pixel,'Method','positive');
    
    % Convert node indices to coordinates
    path_x = zeros(1,numel(path_nodes));
    path_y = zeros(1,numel(path_nodes));
    path = cell(size(path_nodes));
    for i = 1:numel(path_nodes)
        [row, col] = ind2sub(size(map), path_nodes(i));
        path{i} = [col, rows-row];
        path_x(i) = col;
        path_y(i) = rows-row;
    end

    plot(path_x, rows - path_y,'DisplayName',str,'LineWidth',2)
%     text(path_x(end)+10,path_y(end)+10,str,'Color','r','FontSize',12)
%     xlabel('x')
%     ylabel('y')
%     legend()
end

function [max_slope, tot_dist, avg_slope] = main(start,goal,str,elev_intensity,slope_intensity,G)  
    fprintf("From Site to <strong>%s</strong>\n", str)
    fprintf("Coordinates start: [%i, %i]\n", start(1), start(2))
    fprintf("Coordinates end: [%i, %i]\n", goal(1), goal(2))
        
    % Convert start and stop pixel coordinates to nodes
    y_tot = size(elev_intensity,1);
    x_tot = size(elev_intensity,2);
    start_node = sub2ind([y_tot, x_tot], start(2), start(1));
    goal_node = sub2ind([y_tot, x_tot], goal(2), goal(1));
    
    [path, total_cost] = A_star_path_planning_shortestpath(slope_intensity, start_node, goal_node, G, str);
    
    [max_slope,avg_slope,tot_dist] = max_of_path(path,elev_intensity,slope_intensity);
    
    fprintf("   Max slope encountered: %f degrees\n", max_slope)
    fprintf("   Total Distance driven: %f meters\n", tot_dist)
    fprintf("   Average slope: %f degrees\n", avg_slope)
end

function neighbors = get_neighbors(row, col, map_size)
%     disp("In neighbors")
    % Get 8-connected neighbors of a pixel
    rows = [row-1, row-1, row-1, row, row, row+1, row+1, row+1];
    cols = [col-1, col, col+1, col-1, col+1, col-1, col, col+1];
    valid_rows = rows >= 1 & rows <= map_size(1);
    valid_cols = cols >= 1 & cols <= map_size(2);
    valid_indices = valid_rows & valid_cols;
    neighbors = [rows(valid_indices)', cols(valid_indices)'];
end

function total_cost = compute_total_cost(path_nodes, elevation_map)
    disp("In total cost")
    % Compute total cost (elevation change) of the path
    total_cost = 0;
    for i = 1:numel(path_nodes)-1
        [row1, col1] = ind2sub(size(elevation_map), path_nodes(i));
        [row2, col2] = ind2sub(size(elevation_map), path_nodes(i+1));
        total_cost = total_cost + abs(elevation_map(row1, col1) - elevation_map(row2, col2));
    end
end

function G_adjacency_list = extract_graph_as_adjacency_list(G, elevation_map)
    % Convert adjacency matrix G to adjacency list representation
    [m, n] = size(elevation_map);
    G_adjacency_list = cell(m * n, 1);
    
    for i = 1:m * n
        neighbors_indices = find(G(i, :) ~= 0);
        neighbors_costs = G(i, neighbors_indices);
        neighbors_cells = num2cell(neighbors_indices);
        G_adjacency_list{i} = [neighbors_cells; num2cell(neighbors_costs)];
    end
end

function [max_slope,avg_slope,tot_dist] = max_of_path(path, elev_intensity, slope_intensity)
    % Calculate total distance, max slope, and average slope
    y_tot = size(elev_intensity,1);
    x_tot = size(elev_intensity,2);
    tot_dist = 0.0;
    max_slope = 0.0;
    avg_slope = 0.0;
    for i = 1:(size(path,2) - 1)
        curr_xy = path{i};
        next_xy = path{i + 1};
        start_row = y_tot - curr_xy(2) + 1;
        start_col = curr_xy(1);
        goal_row = y_tot - next_xy(2) + 1;
        goal_col = next_xy(1);
        start_node = sub2ind([y_tot, x_tot], start_row, start_col);
        goal_node = sub2ind([y_tot, x_tot], goal_row, goal_col);
        del_elev = abs(elev_intensity(y_tot - curr_xy(2), curr_xy(1)) ...
            - elev_intensity(y_tot - next_xy(2), next_xy(1)));
        curr_slope = slope_intensity(y_tot - curr_xy(2),curr_xy(1));
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
%             disp(next_xy(2));
%             disp(next_xy(1));
        end
    end
    avg_slope = avg_slope / size(path,2);
end