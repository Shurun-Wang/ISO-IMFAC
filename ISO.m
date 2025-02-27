function [Xfood, fval, gbest_t] = ISO(N,T,fobj, dim,lb,ub)
%initial 
vec_flag=[1,-1];  % 常量C取值的随机正负
Threshold=0.2;
Thresold2=0.5;
C1=0.5;
C2=0.05;
C3=2;
dif = ub-lb;

% 生成了40个个体，每个个体包含变量的维度为2。
% dim为变量维度。ub,lb为求解问题的边界值

%生成随机初始化种群
X=lb+rand(N,dim).*(ub-lb);
%生成反向学习种群
X_obl=rand(N,dim).*(ub+lb)-X;
%合并种群
X_total=[X;X_obl];
%根据目标函数值选取最终初始种群
for L=1:size(X_total,1)
    fitness(L)=feval(fobj,X_total(L,:));
end
%对目标函数值进行升序
[B,I] = sort(fitness);
%获取前N个适应度函数值最优的种群为反向学习后的最终初始种群
Best_x=I(1,1:N);
fitness=B(1,1:N);
X=X_total(Best_x,:);
% 注意：randperm的输入是向量的长度  
numRows = size(X, 1);
indices = randperm(numRows);  
% 使用索引向量来重新排列两个向量  
fitness = fitness(1, indices);  
X = X(indices,:);  


% 得到30个个体的适应度值中最小的值及其位置，第一次操作，不需要对比
[GYbest, gbest] = min(fitness); 
% Xfood是食物的位置，也是目前所搜索到最优的位置，包括30个变量的具体值
Xfood = X(gbest,:);
% 种群分离，一半为雌性，一半为雄性
Nm=round(N/2); %eq.(2&3)
Nf=N-Nm;
% 雄性的变量值15*2
Xm=X(1:Nm,:);
% 雌性的变量值15*2
Xf=X(Nm+1:N,:);
% 雄性的适应度15
fitness_m=fitness(1:Nm);
% 雌性的适应度15
fitness_f=fitness(Nm+1:N);
% 雄性中最好的变量值
[fitnessBest_m, gbest1] = min(fitness_m);
Xbest_m = Xm(gbest1,:);
% 雌性中最好的变量值
[fitnessBest_f, gbest2] = min(fitness_f);
Xbest_f = Xf(gbest2,:);
for t = 1:T
    Temp=exp(-((t)/T));     %eq.(4)
    Q=C1*exp(((t-T)/(T)));  %eq.(5)
    if Q>1        
        Q=1;    
    end
    % Exploration Phase (no Food)
    if Q<Threshold
        for i=1:Nm
            % 所有的雄性个体都运动
            for j=1:1:dim
                % 每个个体中的所有变量依次更新位置值
                rand_leader_index = floor(Nm*rand()+1);  % 从1到15随机取个体号
                X_randm = Xm(rand_leader_index, :);
                flag_index = floor(2*rand()+1);
                Flag=vec_flag(flag_index);
                Am=exp(-fitness_m(rand_leader_index)/(fitness_m(i)+eps)); %eq.(7)
                Xnewm(i,j)=X_randm(j)+Flag*C2*Am*(dif(j)*rand+lb(j));       %eq.(6)
            end
        end
        for i=1:Nf
            for j=1:1:dim
                rand_leader_index = floor(Nf*rand()+1);
                X_randf = Xf(rand_leader_index, :);
                flag_index = floor(2*rand()+1);
                Flag=vec_flag(flag_index);
                Af=exp(-fitness_f(rand_leader_index)/(fitness_f(i)+eps));%eq.(9)
                Xnewf(i,j)=X_randf(j)+Flag*C2*Af*(dif(j)*rand+lb(j));%eq.(8)
            end
        end
    else %Exploitation Phase (Food Exists)
        if Temp>Thresold2  %hot
            % Xfood为当前全局最优，所有的个体朝着这里靠近
            for i=1:Nm
                flag_index = floor(2*rand()+1);
                Flag=vec_flag(flag_index);
                for j=1:1:dim
                    Xnewm(i,j)=Xfood(j)+C3*Flag*Temp*rand*(Xfood(j)-Xm(i,j));%eq.(10)
                end
            end
            for i=1:Nf
                flag_index = floor(2*rand()+1);
                Flag=vec_flag(flag_index);
                for j=1:1:dim
                    Xnewf(i,j)=Xfood(j)+Flag*C3*Temp*rand*(Xfood(j)-Xf(i,j));%eq.(10)
                end
            end
        else %cold
            freen = exp(4.*(t/T).^2);
            if rand>0.6 %fight
                % 雄性个体朝着雄性群体中的最优点靠近
                for i=1:Nm
                    for j=1:1:dim
                        FM=exp(-(fitnessBest_f)/(fitness_m(i)+eps));%eq.(13)
                        Xnewm(i,j)=Xm(i,j)+freen*Xm(i,j) +C3*FM*rand*(Q*Xbest_f(j)-Xm(i,j));%eq.(11)                       
                    end
                end
                % 雌性个体朝着雌性群体中的最优点靠近
                for i=1:Nf
                    for j=1:1:dim
                        FF=exp(-(fitnessBest_m)/(fitness_f(i)+eps));%eq.(14)
                        Xnewf(i,j)=Xf(i,j)+freen*Xf(i,j)+C3*FF*rand*(Q*Xbest_m(j)-Xf(i,j));%eq.(12)
                    end
                end
            else%mating
                % 雌性个体朝着雄性个体靠近
                for i=1:Nm
                    for j=1:1:dim
                        Mm=exp(-fitness_f(i)/(fitness_m(i)+eps));%eq.(17)
                        Xnewm(i,j)=Xm(i,j)+freen*Xm(i,j) +C3*rand*Mm*(Q*Xf(i,j)-Xm(i,j));%eq.(15
                    end
                end
                % 雄性个体朝着雌性个体靠近
                for i=1:Nf
                    for j=1:1:dim
                        Mf=exp(-fitness_m(i)/(fitness_f(i)+eps));%eq.(18)
                        Xnewf(i,j)=Xf(i,j)+freen*Xf(i,j) +C3*rand*Mf*(Q*Xm(i,j)-Xf(i,j));%eq.(16)
                    end
                end
                % 是否产卵，如果产卵，则用随机变量替代之前雄性和雌性种群中最差的个体
%                 flag_index = floor(2*rand()+1);
%                 egg=vec_flag(flag_index);
%                 if egg==1;
%                     [GYworst, gworst] = max(fitness_m);
%                     Xnewm(gworst,:)=lb+rand(1,dim).*(ub-lb);%eq.(19)
%                     [GYworst, gworst] = max(fitness_f);
%                     Xnewf(gworst,:)=lb+rand(1,dim).*(ub-lb);%eq.(20)
%                 end
                [~, index]=sort(fitness_m);
                [~, index1]= sort(fitness_f);%排序
                Xnewm(index(end-3),:)=lb+rand*(ub-lb);
                Xnewm(index(end-2),:)=lb+rand*(ub-lb);
                Xnewm(index(end-1),:)=lb+rand*(ub-lb);
                Xnewm(index(end),:)=lb+rand*(ub-lb);
                Xnewf(index1(end-3),:)=lb+rand*(ub-lb);
                Xnewf(index1(end-1),:)=lb+rand*(ub-lb);
                Xnewf(index1(end-2),:)=lb+rand*(ub-lb);
                Xnewf(index1(end),:)=lb+rand*(ub-lb);
            end

        end
    end
    % *************** 最后的更新处理  ************
    for j=1:Nm
        % 越界值处理
        Flag4ub=Xnewm(j,:)>ub;
        Flag4lb=Xnewm(j,:)<lb;
        Xnewm(j,:)=(Xnewm(j,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
        % 计算更新完的个体适应度
        y = feval(fobj,Xnewm(j,:));
        if y<fitness_m(j)
            fitness_m(j)=y;
            Xm(j,:)= Xnewm(j,:);
        end
    end
    % 得到雄性中最好个体的适应度以及个体索引位置
    [Ybest1,gbest1] = min(fitness_m);
    
    for j=1:Nf
         Flag4ub=Xnewf(j,:)>ub;
         Flag4lb=Xnewf(j,:)<lb;
        Xnewf(j,:)=(Xnewf(j,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
        y = feval(fobj,Xnewf(j,:));
        if y<fitness_f(j)
            fitness_f(j)=y;
            Xf(j,:)= Xnewf(j,:);
        end
    end
    % 得到雌性中最好个体的适应度以及个体索引位置    
    [Ybest2,gbest2] = min(fitness_f);
    
    % Ybest1为局部雄性最好的适应度值，fitnessBest_m为全局雄性最好的适应度值
    if Ybest1<fitnessBest_m
        Xbest_m = Xm(gbest1,:);
        fitnessBest_m=Ybest1;
    end
    if Ybest2<fitnessBest_f
        Xbest_f = Xf(gbest2,:);
        fitnessBest_f=Ybest2;
        
    end
    % 当前雄性和雌性最好的适应度值的比较
    if Ybest1<Ybest2
        gbest_t(t)=min(Ybest1);
    else
        gbest_t(t)=min(Ybest2);
        
    end
    % Xfood为当前全局最优
    if fitnessBest_m<fitnessBest_f
        GYbest=fitnessBest_m;
        Xfood=Xbest_m;
    else
        GYbest=fitnessBest_f;
        Xfood=Xbest_f;
    end
    
end
fval = GYbest;
end





