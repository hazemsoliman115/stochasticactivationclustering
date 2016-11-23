clearvars
yalmip('clear')

rng('default');
rng(1);

No_APs = 4;
No_Areas = 16;
No_time_slots = 10;

x_old = zeros(No_time_slots,No_APs);
x_new = zeros(No_time_slots,No_APs);
for t=1:No_time_slots
    x_new(t,:) = [0.5 0.9 0.4 0.1];
end
u_old = zeros(No_time_slots,No_Areas,No_APs);
u_new = (1-0).*rand(No_time_slots,No_Areas,No_APs) + 0;
q_old = (1-0).*zeros(No_time_slots,No_APs,No_APs) + 0;
q_new = (1-0).*rand(No_time_slots,No_APs,No_APs) + 0;

epsilon = 1e-2;



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

sigma_sq = 0.1;

Q_0 = 0.2*ones(1,No_Areas);
Q_e_old = zeros(No_time_slots,No_Areas);
Q_e_old(1,:) = Q_0;
Q_e_new = rand(No_time_slots,No_Areas);
Q_e_new(1,:) = Q_0;
A_e = rand(No_time_slots,No_Areas);
gam = 0.2*ones(1,No_Areas);

alpha_matrix = zeros(No_time_slots,No_Areas,No_APs);
alpha_Q_nom_Pow = zeros(No_time_slots,No_Areas,No_APs);
alpha_cluster = zeros(No_time_slots,No_Areas, No_APs, No_APs);
alpha_Q_den = zeros(No_time_slots,No_Areas, No_APs, No_APs);
alpha_Q_nom_Intf = zeros(No_time_slots,No_Areas, No_APs, No_APs);


while sum(sum(abs(x_new-x_old)))/sum(sum(abs(x_old))) > epsilon
    sum(sum(abs(x_new-x_old)))/sum(sum(abs(x_old)))
    x_old = x_new;
    u_old = u_new;
    q_old = q_new;
    Q_e_old = Q_e_new;
    
    %% Alpha Generation
    tic
    disp = 'Calculating Alphas'
    for t=1:No_time_slots
        den =  zeros(1,No_Areas);
        for i=1:No_Areas
            den(i) = 0;
            for j =1:No_APs
                for k=1:No_APs
                    if k ~= j
                        den(i) = den(i) + u_old(t,i,j)*P(i,k)*x_old(t,k);
                    end
                end
            end
            den(i) = den(i) + sum(reshape(u_old(t,i,:),size(x_old(t,:))).*P(i,:).*x_old(t,:));
            for j =1:No_APs
                alpha_matrix(t,i,j) = ( u_old(t,i,j)*P(i,j)*x_old(t,j) )/( den(i) );
                for k=1:No_APs
                    if k ~= j
                        alpha_cluster(t,i,j,k) = ( gam(i)*Q_e_old(t,i)*u_old(t,i,j)*q_old(t,j,k)*P(i,k)*x_old(t,k) )/( den(i) );
                    end
                end
            end
        end
    end
    for t=2:No_time_slots
        den_Queue_den =  zeros(1,No_Areas);
        for i=1:No_Areas
            den_Queue_den(i) = 0;
            for j =1:No_APs
                for k=1:No_APs
                    if k ~= j
                        den_Queue_den(i) = den_Queue_den(i) + u_old(t-1,i,j)*P(i,k)*x_old(t-1,k);
                    end
                end
            end
            for j =1:No_APs
                for k=1:No_APs
                    if k ~= j
                        alpha_Q_den(t-1,i,j,k) = ( u_old(t-1,i,j)*P(i,k)*x_old(t-1,k) )/( den_Queue_den(i) );
                    end
                end
            end
        end
        den_Queue_nom =  zeros(1,No_Areas);
        for i=1:No_Areas
            den_Queue_nom(i) = 0;
            for j =1:No_APs
                for k=1:No_APs
                    if k ~= j
                        den_Queue_nom(i) = den_Queue_nom(i) + u_old(t-1,i,j)*q_old(t-1,j,k)*P(i,k)*x_old(t-1,k);
                    end
                end
            end
            den_Queue_nom(i) = den_Queue_nom(i) + Q_e_old(t,i)*(Q_e_old(t-1,i))^(-1)/(A_e(t-1,i))*sum(reshape(u_old(t-1,i,:),size(x_old(t-1,:))).*P(i,:).*x_old(t-1,:));
            for j =1:No_APs
                alpha_Q_nom_Pow(t-1,i,j) = ( Q_e_old(t,i)*(Q_e_old(t-1,i))^(-1)/(A_e(t-1,i))*u_old(t-1,i,j)*P(i,j)*x_old(t-1,j) )/( den_Queue_nom(i) );
                for k=1:No_APs
                    if k ~= j
                        alpha_Q_nom_Intf(t-1,i,j,k) = ( u_old(t-1,i,j)*q_old(t-1,j,k)*P(i,k)*x_old(t-1,k) )/( den_Queue_nom(i) );
                    end
                end
            end
        end
    end
    toc
    
    %%
    options = sdpsettings('solver','mosek-geometric');
    x = sdpvar(No_time_slots,No_APs);
    s = sdpvar(No_time_slots,No_Areas);
    Q_e = sdpvar(No_time_slots,No_Areas);
    u = sdpvar(No_time_slots,No_Areas,No_APs);
    q = sdpvar(No_time_slots,No_APs,No_APs,'full');
    z = sdpvar(No_time_slots,No_Areas);
    z_Q_den = sdpvar(No_time_slots,No_Areas);
    z_Q_nom = sdpvar(No_time_slots,No_Areas);
    nom = sdpvar(No_time_slots,No_Areas);
    nom_Q_den = sdpvar(No_time_slots,No_Areas);
    nom_Q_nom = sdpvar(No_time_slots,No_Areas);
    
    
    %minimize( sum(sum(x)) + sum(sum(s)) )
    obj = sum(sum(x)) + sum(sum(s)) ;
    
    tic
    disp = 'Constraint 1'
    Constr_1 = sdpvar(No_time_slots,No_Areas);
    for t=1:No_time_slots
        for i=1:No_Areas
            z(t,i) = 0;
            nom(t,i) = 1;
            for j=1:No_APs
                for k=1:No_APs
                    if k ~= j
                        z(t,i) = z(t,i) + u(t,i,j)*P(i,k)*x(t,k);
                        nom(t,i) = nom(t,i)*(u(t,i,j)*q(t,j,k)*P(i,k)*x(t,k)/alpha_cluster(t,i,j,k))^(-alpha_cluster(t,i,j,k));
                    end
                end
            end
            %(gam(i)*Q_e(t,i))*z(t,i)*nom(t,i)*prod( (reshape(u(t,i,:),size(x(t,:))).*P(i,:).*x(t,:)./reshape(alpha_matrix(t,i,:),size(x(t,:)))).^(-reshape(alpha_matrix(t,i,:),size(x(t,:)))) ) <= 1;
            %Constr = [ Constr , (gam(i)*Q_e(t,i))*z(t,i)*nom(t,i)*prod( (reshape(u(t,i,:),size(x(t,:))).*P(i,:).*x(t,:)./reshape(alpha_matrix(t,i,:),size(x(t,:)))).^(-reshape(alpha_matrix(t,i,:),size(x(t,:)))) ) <= 1];
            Constr_1(t,i) = (gam(i)*Q_e(t,i))*z(t,i)*nom(t,i)*prod( (reshape(u(t,i,:),size(x(t,:))).*P(i,:).*x(t,:)./reshape(alpha_matrix(t,i,:),size(x(t,:)))).^(-reshape(alpha_matrix(t,i,:),size(x(t,:)))) );
        end
    end
    toc
    
    tic
    disp = 'Constraint 2'
    Constr_2 = sdpvar(No_time_slots-1,No_Areas);
    for t=2:No_time_slots
        for i=1:No_Areas
            z_Q_den(t-1,i) = 0;
            nom_Q_den(t-1,i) = 1;
            for j=1:No_APs
                for k=1:No_APs
                    if k ~= j
                        z_Q_den(t-1,i) = z_Q_den(t-1,i) + u(t-1,i,j)*q(t-1,j,k)*P(i,k)*x(t-1,k);
                        nom_Q_den(t-1,i) = nom_Q_den(t-1,i)*(u(t-1,i,j)*P(i,k)*x(t-1,k)/alpha_Q_den(t-1,i,j,k))^(-alpha_Q_den(t-1,i,j,k));
                    end
                end
            end
            %( Q_e(t,i)*(Q_e(t-1,i))^(-1)/(A_e(t-1,i))*sum(reshape(u(t-1,i,:),size(x(t-1,:))).*P(i,:).*x(t-1,:)) + z_Q_den(t-1,i) )*nom_Q_den(t-1,i)  <= 1;
            %Constr = [ Constr , ( Q_e(t,i)*(Q_e(t-1,i))^(-1)/(A_e(t-1,i))*sum(reshape(u(t-1,i,:),size(x(t-1,:))).*P(i,:).*x(t-1,:)) + z_Q_den(t-1,i) )*nom_Q_den(t-1,i)  <= 1];
            Constr_2(t-1,i) = ( Q_e(t,i)*(Q_e(t-1,i))^(-1)/(A_e(t-1,i))*sum(reshape(u(t-1,i,:),size(x(t-1,:))).*P(i,:).*x(t-1,:)) + z_Q_den(t-1,i) )*nom_Q_den(t-1,i);
        end
    end
    toc
    
    tic
    disp = 'Constraint 3'
    Constr_3 = sdpvar(No_time_slots-1,No_Areas);
    for t=2:No_time_slots
        for i=1:No_Areas
            z_Q_nom(t-1,i) = 0;
            nom_Q_nom(t-1,i) = 1;
            for j=1:No_APs
                for k=1:No_APs
                    if k ~= j
                        z_Q_nom(t-1,i) = z_Q_nom(t-1,i) + u(t-1,i,j)*P(i,k)*x(t-1,k);
                        nom_Q_nom(t-1,i) = nom_Q_nom(t-1,i)*(u(t-1,i,j)*q(t-1,j,k)*P(i,k)*x(t-1,k)/alpha_Q_nom_Intf(t-1,i,j,k))^(-alpha_Q_nom_Intf(t-1,i,j,k));
                    end
                end
            end
            %s(t,i)^(-1)*z_Q_nom(t-1,i)*nom_Q_nom(t-1,i)*prod( (Q_e(t,i)*(Q_e(t-1,i))^(-1)/(A_e(t-1,i))*reshape(u(t-1,i,:),size(x(t-1,:))).*P(i,:).*x(t-1,:)./reshape(alpha_Q_nom_Pow(t-1,i,:),size(x(t-1,:)))).^(-reshape(alpha_Q_nom_Pow(t-1,i,:),size(x(t-1,:)))) )*nom_Q_nom(t-1,i) <= 1;
            %Constr = [ Constr , s(t,i)^(-1)*z_Q_nom(t-1,i)*nom_Q_nom(t-1,i)*prod( (Q_e(t,i)*(Q_e(t-1,i))^(-1)/(A_e(t-1,i))*reshape(u(t-1,i,:),size(x(t-1,:))).*P(i,:).*x(t-1,:)./reshape(alpha_Q_nom_Pow(t-1,i,:),size(x(t-1,:)))).^(-reshape(alpha_Q_nom_Pow(t-1,i,:),size(x(t-1,:)))) )*nom_Q_nom(t-1,i) <= 1];
            Constr_3(t-1,i) = s(t,i)^(-1)*z_Q_nom(t-1,i)*nom_Q_nom(t-1,i)*prod( (Q_e(t,i)*(Q_e(t-1,i))^(-1)/(A_e(t-1,i))*reshape(u(t-1,i,:),size(x(t-1,:))).*P(i,:).*x(t-1,:)./reshape(alpha_Q_nom_Pow(t-1,i,:),size(x(t-1,:)))).^(-reshape(alpha_Q_nom_Pow(t-1,i,:),size(x(t-1,:)))) )*nom_Q_nom(t-1,i);
        end
    end
    toc
    
%     Constr = [ Constr , 0.1 <= x <= 1];
%     Constr = [ Constr , 1 <= s ];
%     Constr = [ Constr , 0.1 <= u <= 1];
%     Constr = [ Constr , 0.1 <= q <= 1];
%     Constr = [ Constr , 0.1 <= Q_e];
    tic
    disp = 'Constraint combine'
    Constr = [Constr_1 <= 1, Constr_2 <= 1, Constr_3 <= 1, 0.1 <= x <= 1, 1 <= s, 0.1 <= u <= 1, 0.1 <= q <= 1, 0.1 <= Q_e];
    toc
     for t=1:No_time_slots
         for i=1:No_Areas
             Constr = [ Constr , sum(u(t,i,:)) <= 1];
         end
     end
    tic
    sol = optimize(Constr,obj,options);
    toc
    
    x_new = value(x)
    u_new = value(u)
    q_new = value(q)
    s_new = value(s)
    Q_e_new(2,:) = value(Q_e(2,:))
    yalmip('clear')
end

