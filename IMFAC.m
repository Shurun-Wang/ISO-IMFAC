function objval = IMFAC(x)%[-5.12,5.12]
exp = 'E1';
%控制器参数
N=150;
eita = x(1);%1;
rou=x(2);%1;
miu = x(3);%0.4;
init_fai = x(4);%0.4;
fai(1:2) = init_fai;
lamda = x(5);%1;

y_list = [];
fai_list = [];
% *************************************************************
%初值
y(1:3)=0;
u(1:3)=0;
du(1:3) = 0;

%控制器伪偏导数初值
if strcmp(exp,'E1')
    [d1, d2, d3] = deal(0.5, 1, 2);
elseif strcmp(exp,'E2')
    [d1, d2, d3] = deal(0.5, 2, 1);
elseif strcmp(exp,'E3')
    [d1, d2, d3] = deal(1, 0.5, 2);
elseif strcmp(exp,'E4')
    [d1, d2, d3] = deal(1, 2, 0.5);
elseif strcmp(exp,'E5')
    [d1, d2, d3] = deal(2, 0.5, 1);
elseif strcmp(exp,'E6')
    [d1, d2, d3] = deal(2, 1, 0.5);
end

for k=1:N+1
    if k< 50
       refs(k)=d1;
    elseif k< 100
       refs(k)=d2;
    else
       refs(k)=d3;
    end
end

for k=3:N-1
    fai(k)=fai(k-1)+eita*(y(k)-y(k-1)-fai(k-1)*du(k-1))*du(k-1)/(miu+du(k-1)*du(k-1));
    fai_list = [fai_list fai(k)];
    %% 防止分母太小 公式4.7，4.8
    if (abs(fai(k))<10^(-5)) || ((du(k-1)*du(k-1))^0.5<10^(-5))
        fai(k)=init_fai;
    end
    %% 输入信号u的变化，也就是系统的控制方案，公式4.4  lamda限制了控制输出的变化，非常有用
    u(k) = u(k-1)+rou*fai(k)*(x(6)*(refs(k+1)-y(k)-refs(k)+y(k-1))+ x(7)*(refs(k+1)-y(k))+x(8)*(refs(k+1)-y(k)-2*refs(k)+2*y(k-1)+refs(k-1)-y(k-2)))/(lamda+fai(k).^2); % 只有一个fai     

    y(k+1)=0.6*y(k)-0.1*y(k-1)+1.8*u(k)-1.8*u(k)*u(k)+0.6*u(k)*u(k)*u(k)-0.15*u(k-1)+0.15*u(k-1)*u(k-1)-0.05*u(k-1)*u(k-1)*u(k-1); 
    y_list = [y_list y(k+1)];
    du(k)=u(k)-u(k-1);
end

y_list = [0 0 0 y_list];
refs = refs(1:end-1);
k = 1:1:N;
if isnan(sum(abs(y_list-refs)))
    objval = 1e8;
else
    objval_1 = sum(abs(y_list-refs)); % IAE
    objval_2 = sum(k.*abs(y_list-refs))/N; % ITAE
    objval_3 = max(abs(y_list-refs));
    objval_4 = min(abs(y_list-refs));
    objval_5 = mean(abs(y_list-refs));
    objval_6 = std(abs(y_list-refs));
    objval = objval_1 + objval_2 + objval_3 + objval_4 + objval_5 + objval_6;
end
% 
% % step=1;
% figure(1)
% % plot(0,'-k','MarkerSize',mark,'LineWidth',2);hold on;
% % plot(0,'-.bs','MarkerSize',mark,'LineWidth',2);hold on;
% set(gca,'LineWidth',1,'fontsize',28);
% plot(refs,'k','LineWidth',1);hold on;
% plot(y);hold on;grid on;
% % plot(1:step:N,y(1:step:N),'bs','MarkerSize',mark,'LineWidth',2);hold on;
% grid on;xlabel('时刻');ylabel('跟踪性能');legend({'y^{*}(k)','y(k)'},'Interpreter','tex');
