function [Sv F r_propose] = CalS(idx, Kv,M,NRF,Tx_CB,Rcron,PT,miu,weight,rankRcron)
group_nums=length(Kv);
Sv = zeros(group_nums,1);
add_flag=1;
Sv_star=Sv;
[R_star F r_propose]=calrateS(idx, Sv_star,Kv,M,NRF,Tx_CB,Rcron,PT,miu,weight);
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
