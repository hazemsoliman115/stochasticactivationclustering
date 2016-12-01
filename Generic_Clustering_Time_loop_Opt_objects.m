clearvars
addpath(genpath('E:\Research/StochasticActivationClustering'))
yalmip('clear')
path_config
%yalmiptest(sdpsettings('solver','mosek-geometric'))

rng('default');
rng(1);

%% Parameters
inter_ap_distance = 15;
inter_area_distance = 7.5;
epsilon = 0.1;
sigma_sq = 0.1;
No_APs = 4;
No_Areas = 16;
No_time_slots = 10;

%% create list of access points and areas objects

for i =1:sqrt(No_APs)
    for j =1:sqrt(No_APs)
        AP_array(i,j) = accesspoint();
    end
end


for i =1:sqrt(No_Areas)
    for j =1:sqrt(No_Areas)
        AArea_array(i,j) = access_area();
    end
end

%% Add the coordinates
for i =1:sqrt(No_APs)
    for j =1:sqrt(No_APs)
        AP_array(i,j).x_pos = inter_ap_distance/2 + (i-1)*inter_ap_distance;
        AP_array(i,j).y_pos = inter_ap_distance/2 + (j-1)*inter_ap_distance;
    end
end
AP_flat_array = reshape(AP_array,[No_APs,1]);

for i =1:sqrt(No_Areas)
    for j =1:sqrt(No_Areas)
        AArea_array(i,j).x_pos = inter_area_distance/2 + (i-1)*inter_area_distance;
        AArea_array(i,j).y_pos = inter_area_distance/2 + (j-1)*inter_area_distance;
        AArea_array(i,j).queue_load = rand(No_time_slots,1);
    end
end
AArea_flat_array = reshape(AArea_array,[No_Areas,1]);

%% Calculate Power Matrix
p_dBm = 80;
p_mat = zeros(No_Areas,No_APs);
for j = 1 : No_Areas
    for i = 1 : No_APs
        p_mat(j,i) = get_path_loss(AP_flat_array(i),AArea_flat_array(j), p_dBm);
    end
end
noise_power_dBm = -174;
noise_power_watt  = 10^((noise_power_dBm-30)/10);


%% Set up Traffic Qeue
A_e = zeros(No_time_slots,No_Areas);
for j = 1 : No_Areas
        A_e(:,j) = AArea_flat_array(j).queue_load;
end
% A_e =1*[    0.9935    0.5396    0.1428    0.8125    0.3258    0.6704    0.5121 0.0509    0.0433    0.2358    0.8146    0.0644    0.2821    0.6022 0.1926    0.3210
%     0.8345    0.8953    0.0941    0.2838    0.8898    0.4305    0.6175 0.2351    0.4864    0.9408    0.0541    0.8970    0.8838    0.9636 0.7529    0.1230
%     0.6996    0.4466    0.8702    0.5278    0.7517    0.7678    0.4324 0.0633    0.2394    0.6842    0.1305    0.2034    0.5677    0.3457 0.0073    0.7213
%     0.9183    0.8770    0.2369    0.3394    0.7626    0.5360    0.8477 0.4217    0.9525    0.0649    0.8424    0.8262    0.1151    0.5956 0.3283    0.4403
%     0.0397    0.2536    0.3860    0.5547    0.4695    0.0399    0.4541 0.8638    0.9439    0.8704    0.6183    0.8818    0.2270    0.5990 0.9176    0.1267
%     0.0703    0.2738    0.5715    0.9744    0.2108    0.1348    0.0154 0.0816    0.6139    0.7014    0.5313    0.4868    0.5960    0.6157 0.5884    0.5898
%     0.4740    0.3284    0.5258    0.3117    0.0415    0.1934    0.8731 0.4731    0.9735    0.6049    0.2483    0.5985    0.2394    0.0592 0.8552    0.0361
%     0.3492    0.5476    0.0760    0.6688    0.3218    0.3357    0.6562 0.1255    0.3449    0.7324    0.2951    0.5273    0.1314    0.7503 0.6047    0.2002
%     0.9373    0.2201    0.8741    0.3260    0.0371    0.0523    0.8230 0.7729    0.8979    0.2534    0.8727    0.6248    0.1618    0.9482 0.8226    0.7883
%     0.4896    0.6714    0.9511    0.7745    0.6939    0.6051    0.9518 0.8414    0.4346    0.6005    0.4217    0.8550    0.8449    0.5347 0.8795    0.0121];

A_e = A_e + 1;

%% Create the traffic
traffic_load = zeros(No_time_slots,1);
for t = 1:No_time_slots
    traffic_load(t) = t;
end

x_old = zeros(No_time_slots,No_APs);
x_new = zeros(No_time_slots,No_APs);
for t=1:No_time_slots
    x_new(t,:) = [0.5 0.9 0.4 0.1];
end
u_old = zeros(No_time_slots,No_Areas,No_APs);
u_new = (1-0).*rand(No_time_slots,No_Areas,No_APs) + 0;
q_old = (1-0).*zeros(No_time_slots,No_APs,No_APs) + 0;
q_new = (1-0).*rand(No_time_slots,No_APs,No_APs) + 0;

P = zeros(No_Areas,No_APs);

for i=1:No_Areas
    for j =1:No_APs
        if (1<=i)&&(i<=4)
            if j==1
                P(i,j) = 1;
            elseif (j==2) || (j==3)
                P(i,j) = 0.5;
            else
                P(i,j) = 0.25;
            end
        elseif (5<=i)&&(i<=8)
            if j==2
                P(i,j) = 1;
            elseif (j==1) || (j==3)
                P(i,j) = 0.5;
            else
                P(i,j) = 0.25;
            end
        elseif (9<=i)&&(i<=12)
            if j==3
                P(i,j) = 1;
            elseif (j==1) || (j==4)
                P(i,j) = 0.5;
            else
                P(i,j) = 0.25;
            end
        elseif (13<=i)&&(i<=16)
            if j==4
                P(i,j) = 1;
            elseif (j==2) || (j==3)
                P(i,j) = 0.5;
            else
                P(i,j) = 0.25;
            end
        end
    end
end

P = p_mat;

sigma_sq = 0.1;



Q_0 = 0.2*ones(1,No_Areas);
Q_e_old = zeros(No_time_slots,No_Areas);
Q_e_old(1,:) = Q_0;
Q_e_new = rand(No_time_slots,No_Areas);
Q_e_new(1,:) = Q_0;
%A_e = rand(No_time_slots,No_Areas);
gam = 2*ones(1,No_Areas);

[x_new, u_new, q_new, Q_e_new] = Generic_Clustering_Time_loop_Opt_func_edited(No_time_slots,No_APs,No_Areas,P, sigma_sq, A_e, gam);

filename = date;
save(filename,'x_new', 'u_new', 'q_new', 'Q_e_new', 'A_e', 'P','inter_ap_distance', 'inter_area_distance', 'No_APs', 'No_Areas', 'No_time_slots')



% alpha_matrix = zeros(No_time_slots,No_Areas,No_APs);
% alpha_Q_nom_Pow = zeros(No_time_slots,No_Areas,No_APs);
% alpha_cluster = zeros(No_time_slots,No_Areas, No_APs, No_APs);
% alpha_Q_den = zeros(No_time_slots,No_Areas, No_APs, No_APs);
% alpha_Q_nom_Intf = zeros(No_time_slots,No_Areas, No_APs, No_APs);
% 
% 
% while sum(sum(abs(x_new-x_old)))/sum(sum(abs(x_old))) > epsilon
%     sum(sum(abs(x_new-x_old)))/sum(sum(abs(x_old)))
%     x_old = x_new;
%     u_old = u_new;
%     q_old = q_new;
%     Q_e_old = Q_e_new;
%     
%     %% Alpha Generation
%     tic
%     disp = 'Calculating Alphas'
%     for t=1:No_time_slots
%         den =  zeros(1,No_Areas);
%         for i=1:No_Areas
%             den(i) = 0;
%             for j =1:No_APs
%                 for k=1:No_APs
%                     if k ~= j
%                         den(i) = den(i) + u_old(t,i,j)*P(i,k)*x_old(t,k);
%                     end
%                 end
%             end
%             den(i) = den(i) + sum(reshape(u_old(t,i,:),size(x_old(t,:))).*P(i,:).*x_old(t,:));
%             for j =1:No_APs
%                 alpha_matrix(t,i,j) = ( u_old(t,i,j)*P(i,j)*x_old(t,j) )/( den(i) );
%                 for k=1:No_APs
%                     if k ~= j
%                         alpha_cluster(t,i,j,k) = ( gam(i)*Q_e_old(t,i)*u_old(t,i,j)*q_old(t,j,k)*P(i,k)*x_old(t,k) )/( den(i) );
%                     end
%                 end
%             end
%         end
%     end
%     for t=2:No_time_slots
%         den_Queue_den =  zeros(1,No_Areas);
%         for i=1:No_Areas
%             den_Queue_den(i) = 0;
%             for j =1:No_APs
%                 for k=1:No_APs
%                     if k ~= j
%                         den_Queue_den(i) = den_Queue_den(i) + u_old(t-1,i,j)*P(i,k)*x_old(t-1,k);
%                     end
%                 end
%             end
%             for j =1:No_APs
%                 for k=1:No_APs
%                     if k ~= j
%                         alpha_Q_den(t-1,i,j,k) = ( u_old(t-1,i,j)*P(i,k)*x_old(t-1,k) )/( den_Queue_den(i) );
%                     end
%                 end
%             end
%         end
%         den_Queue_nom =  zeros(1,No_Areas);
%         for i=1:No_Areas
%             den_Queue_nom(i) = 0;
%             for j =1:No_APs
%                 for k=1:No_APs
%                     if k ~= j
%                         den_Queue_nom(i) = den_Queue_nom(i) + u_old(t-1,i,j)*q_old(t-1,j,k)*P(i,k)*x_old(t-1,k);
%                     end
%                 end
%             end
%             den_Queue_nom(i) = den_Queue_nom(i) + Q_e_old(t,i)*(Q_e_old(t-1,i))^(-1)/(A_e(t-1,i))*sum(reshape(u_old(t-1,i,:),size(x_old(t-1,:))).*P(i,:).*x_old(t-1,:));
%             for j =1:No_APs
%                 alpha_Q_nom_Pow(t-1,i,j) = ( Q_e_old(t,i)*(Q_e_old(t-1,i))^(-1)/(A_e(t-1,i))*u_old(t-1,i,j)*P(i,j)*x_old(t-1,j) )/( den_Queue_nom(i) );
%                 for k=1:No_APs
%                     if k ~= j
%                         alpha_Q_nom_Intf(t-1,i,j,k) = ( u_old(t-1,i,j)*q_old(t-1,j,k)*P(i,k)*x_old(t-1,k) )/( den_Queue_nom(i) );
%                     end
%                 end
%             end
%         end
%     end
%     toc
%     
%     %%
%     options = sdpsettings('solver','mosek-geometric');
%     x = sdpvar(No_time_slots,No_APs);
%     s = sdpvar(No_time_slots,No_Areas);
%     Q_e = sdpvar(No_time_slots,No_Areas);
%     u = sdpvar(No_time_slots,No_Areas,No_APs);
%     q = sdpvar(No_time_slots,No_APs,No_APs,'full');
%     z = sdpvar(No_time_slots,No_Areas);
%     z_Q_den = sdpvar(No_time_slots,No_Areas);
%     z_Q_nom = sdpvar(No_time_slots,No_Areas);
%     nom = sdpvar(No_time_slots,No_Areas);
%     nom_Q_den = sdpvar(No_time_slots,No_Areas);
%     nom_Q_nom = sdpvar(No_time_slots,No_Areas);
%     
%     
%     %minimize( sum(sum(x)) + sum(sum(s)) )
%     obj = sum(sum(x)) + sum(sum(s)) ;
%     
%     tic
%     disp = 'Constraint 1'
%     Constr_1 = sdpvar(No_time_slots,No_Areas);
%     for t=1:No_time_slots
%         for i=1:No_Areas
%             z(t,i) = 0;
%             nom(t,i) = 1;
%             for j=1:No_APs
%                 for k=1:No_APs
%                     if k ~= j
%                         z(t,i) = z(t,i) + u(t,i,j)*P(i,k)*x(t,k);
%                         nom(t,i) = nom(t,i)*(u(t,i,j)*q(t,j,k)*P(i,k)*x(t,k)/alpha_cluster(t,i,j,k))^(-alpha_cluster(t,i,j,k));
%                     end
%                 end
%             end
%             %(gam(i)*Q_e(t,i))*z(t,i)*nom(t,i)*prod( (reshape(u(t,i,:),size(x(t,:))).*P(i,:).*x(t,:)./reshape(alpha_matrix(t,i,:),size(x(t,:)))).^(-reshape(alpha_matrix(t,i,:),size(x(t,:)))) ) <= 1;
%             %Constr = [ Constr , (gam(i)*Q_e(t,i))*z(t,i)*nom(t,i)*prod( (reshape(u(t,i,:),size(x(t,:))).*P(i,:).*x(t,:)./reshape(alpha_matrix(t,i,:),size(x(t,:)))).^(-reshape(alpha_matrix(t,i,:),size(x(t,:)))) ) <= 1];
%             Constr_1(t,i) = (gam(i)*Q_e(t,i))*z(t,i)*nom(t,i)*prod( (reshape(u(t,i,:),size(x(t,:))).*P(i,:).*x(t,:)./reshape(alpha_matrix(t,i,:),size(x(t,:)))).^(-reshape(alpha_matrix(t,i,:),size(x(t,:)))) );
%         end
%     end
%     toc
%     
%     tic
%     disp = 'Constraint 2'
%     Constr_2 = sdpvar(No_time_slots-1,No_Areas);
%     for t=2:No_time_slots
%         for i=1:No_Areas
%             z_Q_den(t-1,i) = 0;
%             nom_Q_den(t-1,i) = 1;
%             for j=1:No_APs
%                 for k=1:No_APs
%                     if k ~= j
%                         z_Q_den(t-1,i) = z_Q_den(t-1,i) + u(t-1,i,j)*q(t-1,j,k)*P(i,k)*x(t-1,k);
%                         nom_Q_den(t-1,i) = nom_Q_den(t-1,i)*(u(t-1,i,j)*P(i,k)*x(t-1,k)/alpha_Q_den(t-1,i,j,k))^(-alpha_Q_den(t-1,i,j,k));
%                     end
%                 end
%             end
%             %( Q_e(t,i)*(Q_e(t-1,i))^(-1)/(A_e(t-1,i))*sum(reshape(u(t-1,i,:),size(x(t-1,:))).*P(i,:).*x(t-1,:)) + z_Q_den(t-1,i) )*nom_Q_den(t-1,i)  <= 1;
%             %Constr = [ Constr , ( Q_e(t,i)*(Q_e(t-1,i))^(-1)/(A_e(t-1,i))*sum(reshape(u(t-1,i,:),size(x(t-1,:))).*P(i,:).*x(t-1,:)) + z_Q_den(t-1,i) )*nom_Q_den(t-1,i)  <= 1];
%             Constr_2(t-1,i) = ( Q_e(t,i)*(Q_e(t-1,i))^(-1)/(A_e(t-1,i))*sum(reshape(u(t-1,i,:),size(x(t-1,:))).*P(i,:).*x(t-1,:)) + z_Q_den(t-1,i) )*nom_Q_den(t-1,i);
%         end
%     end
%     toc
%     
%     tic
%     disp = 'Constraint 3'
%     Constr_3 = sdpvar(No_time_slots-1,No_Areas);
%     for t=2:No_time_slots
%         for i=1:No_Areas
%             z_Q_nom(t-1,i) = 0;
%             nom_Q_nom(t-1,i) = 1;
%             for j=1:No_APs
%                 for k=1:No_APs
%                     if k ~= j
%                         z_Q_nom(t-1,i) = z_Q_nom(t-1,i) + u(t-1,i,j)*P(i,k)*x(t-1,k);
%                         nom_Q_nom(t-1,i) = nom_Q_nom(t-1,i)*(u(t-1,i,j)*q(t-1,j,k)*P(i,k)*x(t-1,k)/alpha_Q_nom_Intf(t-1,i,j,k))^(-alpha_Q_nom_Intf(t-1,i,j,k));
%                     end
%                 end
%             end
%             %s(t,i)^(-1)*z_Q_nom(t-1,i)*nom_Q_nom(t-1,i)*prod( (Q_e(t,i)*(Q_e(t-1,i))^(-1)/(A_e(t-1,i))*reshape(u(t-1,i,:),size(x(t-1,:))).*P(i,:).*x(t-1,:)./reshape(alpha_Q_nom_Pow(t-1,i,:),size(x(t-1,:)))).^(-reshape(alpha_Q_nom_Pow(t-1,i,:),size(x(t-1,:)))) )*nom_Q_nom(t-1,i) <= 1;
%             %Constr = [ Constr , s(t,i)^(-1)*z_Q_nom(t-1,i)*nom_Q_nom(t-1,i)*prod( (Q_e(t,i)*(Q_e(t-1,i))^(-1)/(A_e(t-1,i))*reshape(u(t-1,i,:),size(x(t-1,:))).*P(i,:).*x(t-1,:)./reshape(alpha_Q_nom_Pow(t-1,i,:),size(x(t-1,:)))).^(-reshape(alpha_Q_nom_Pow(t-1,i,:),size(x(t-1,:)))) )*nom_Q_nom(t-1,i) <= 1];
%             Constr_3(t-1,i) = s(t,i)^(-1)*z_Q_nom(t-1,i)*nom_Q_nom(t-1,i)*prod( (Q_e(t,i)*(Q_e(t-1,i))^(-1)/(A_e(t-1,i))*reshape(u(t-1,i,:),size(x(t-1,:))).*P(i,:).*x(t-1,:)./reshape(alpha_Q_nom_Pow(t-1,i,:),size(x(t-1,:)))).^(-reshape(alpha_Q_nom_Pow(t-1,i,:),size(x(t-1,:)))) )*nom_Q_nom(t-1,i);
%         end
%     end
%     toc
%     
% %     Constr = [ Constr , 0.1 <= x <= 1];
% %     Constr = [ Constr , 1 <= s ];
% %     Constr = [ Constr , 0.1 <= u <= 1];
% %     Constr = [ Constr , 0.1 <= q <= 1];
% %     Constr = [ Constr , 0.1 <= Q_e];
%     tic
%     disp = 'Constraint combine'
%     Constr = [Constr_1 <= 1, Constr_2 <= 1, Constr_3 <= 1, 0.1 <= x <= 1, 1 <= s, 0.1 <= u <= 1, 0.1 <= q <= 1, 0.1 <= Q_e];
%     toc
%      for t=1:No_time_slots
%          for i=1:No_Areas
%              Constr = [ Constr , sum(u(t,i,:)) <= 1];
%          end
%      end
%     tic
%     sol = optimize(Constr,obj,options);
%     toc
%     
%     x_new = value(x)
%     u_new = value(u)
%     q_new = value(q)
%     s_new = value(s)
%     Q_e_new(2,:) = value(Q_e(2,:))
%     yalmip('clear')
% end
% 
