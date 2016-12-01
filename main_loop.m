clearvars
addpath(genpath('E:\Research/StochasticActivationClustering'))
yalmip('clear')
path_config
%yalmiptest(sdpsettings('solver','mosek-geometric'))

rng('default');
rng(1);


%% Parameters
no_accesspoints = 4;
no_accessareas = 16;
inter_ap_distance = 20;
inter_area_distance = 10;
no_time_slots = 6;
epsilon = 0.1;
sigma_sq = 0.1;
%% create list of access points and areas objects

for i =1:sqrt(no_accesspoints)
    for j =1:sqrt(no_accesspoints)
        AP_array(i,j) = accesspoint();
    end
end


for i =1:sqrt(no_accessareas)
    for j =1:sqrt(no_accessareas)
        AArea_array(i,j) = access_area();
    end
end

%% Add the coordinates
for i =1:sqrt(no_accesspoints)
    for j =1:sqrt(no_accesspoints)
        AP_array(i,j).x_pos = inter_ap_distance/2 + (i-1)*inter_ap_distance;
        AP_array(i,j).y_pos = inter_ap_distance/2 + (j-1)*inter_ap_distance;
    end
end
AP_flat_array = reshape(AP_array,[no_accesspoints,1]);

for i =1:sqrt(no_accessareas)
    for j =1:sqrt(no_accessareas)
        AArea_array(i,j).x_pos = inter_area_distance/2 + (i-1)*inter_area_distance;
        AArea_array(i,j).y_pos = inter_area_distance/2 + (j-1)*inter_area_distance;
        AArea_array(i,j).queue_load = rand(no_time_slots,1);
    end
end
AArea_flat_array = reshape(AArea_array,[no_accessareas,1]);

%% Calculate Power Matrix
p_dBm = 80;
p_mat = zeros(no_accessareas,no_accesspoints);
for j = 1 : no_accessareas
    for i = 1 : no_accesspoints
        p_mat(j,i) = get_path_loss(AP_flat_array(i),AArea_flat_array(j), p_dBm);
    end
end
noise_power_dBm = -174;
noise_power_watt  = 10^((noise_power_dBm-30)/10);


%% Set up Traffic Qeue
A_e = zeros(no_time_slots,no_accessareas);
for j = 1 : no_accessareas
        A_e(:,j) = AArea_flat_array(j).queue_load;
end


%% Create the traffic
traffic_load = zeros(no_time_slots,1);
for t = 1:no_time_slots
    traffic_load(t) = t;
end

%% Call the optimization
[x_new, u_new, q_new] = generic_cluster_opt(no_accesspoints, no_accessareas, no_time_slots, epsilon, p_mat, sigma_sq, A_e)
%% Save the data