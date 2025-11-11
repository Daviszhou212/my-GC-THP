function [Sv F r_propose] = CalS(idx, Kv,M,NRF,Tx_CB,Rcron,PT,miu,weight,rankRcron)
% 给定weight向量，最大化WSR
% 输出：
% 1、组调度向量Sv（每个用户组调度几个用户，在这个函数里求得）
% 2、选取的模拟波束F（在calrateS函数里求得）
% 3、WSR优化后的各个用户速率向量（速率DE，在CalrateG函数里求得）
group_nums=length(Kv);
Sv = zeros(group_nums,1);
add_flag=1;
Sv_star=Sv;
[R_star F r_propose]=calrateS(idx, Sv_star,Kv,M,NRF,Tx_CB,Rcron,PT,miu,weight);
% 从Sv = 0开始，然后在一个while循环中，一次增加一个数据流，检查最大化WSR的分配
while(add_flag==1 & sum(Sv)<NRF)
    add_flag = 0;
    for kj=1:group_nums
        Sv_temp=Sv;
        Sv_temp(kj)=Sv_temp(kj)+1;
        [R_temp F_temp r_temp]=calrateS(idx, Sv_temp,Kv,M,NRF,Tx_CB,Rcron,PT,miu,weight);
        if Sv_temp(kj)<=rankRcron(kj) & R_temp > R_star
            Sv_star=Sv_temp;
            R_star = R_temp;
            F = F_temp;
            r_propose = r_temp;
            add_flag = 1;
        end;
    end;
   Sv=Sv_star;
end;
