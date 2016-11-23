clearvars

rng('default');
rng(1);

No_APs = 4;
No_Areas = 16;

x_old = zeros(1,No_APs);
x_new = zeros(1,No_APs);
x_new(1,:) = [0.5 0.9 0.4 0.1];
u_old = zeros(16,4);
u_new = (1-0).*rand(16,4) + 0;

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

Q = 0.2*ones(1,No_Areas);
gam = 0.2*ones(1,No_Areas);

alpha_matrix = zeros(No_Areas,No_APs);
iter_no = 0;
while sum(abs(x_new-x_old)) > epsilon
    iter_no = iter_no + 1;
    x_old = x_new;
    u_old = u_new;
    
    %alpha_1a = ( u_old(1,1)*P_11*x_old(1) )/( u_old(1,1)*P_11*x_old(1) + u_old(1,2)*P_12*x_old(2) + u_old(1,3)*P_13*x_old(3) + u_old(1,4)*P_14*x_old(4) );
    for i=1:No_Areas
        for j =1:No_APs
            alpha_matrix(i,j) = ( u_old(i,j)*P(i,j)*x_old(j) )/(sum(u_old(i,:).*P(i,:).*x_old(1,:)));
        end
    end
    
    options = sdpsettings('solver','mosek-geometric');
    x = sdpvar(1,No_APs);
    u = sdpvar(No_Areas,No_APs);
    %expression z(No_Areas);
    z = sdpvar(No_Areas);
    %minimize( sum(x) )
    obj = sum(x);
    C = [];
    
    for i=1:No_Areas
        z(i) = 0;
        for j=1:No_APs
            for k=1:No_APs
                if k ~= j
                    z(i) = z(i) + u(i,j)*P(i,k)*x(1,k);
                end
            end
        end
        %(gam(i)*Q(i))*z(i)*prod( (u(i,:).*P(i,:).*x(1,:)./alpha_matrix(i,:)).^(-alpha_matrix(i,:)) ) <= 1;
        C = [C, (gam(i)*Q(i))*z(i)*prod( (u(i,:).*P(i,:).*x(1,:)./alpha_matrix(i,:)).^(-alpha_matrix(i,:)) ) <= 1];
        
    end
    
    
    C = [C, 0.1*ones(1,4) <= x <= ones(1,4)];
    C = [C, 0.1*ones(16,4) <= u <= ones(16,4)];
    for i=1:16
        C = [C, sum(u(i,:)) <= 1];
    end
    sol = optimize(C,obj,options);
    x_new = value(x);
    u_new = value(u);

end

