cvx_begin quiet % 'quiet' 模式可以抑制求解过程的输出
cvx_solver Mosek_5
    variable qo(kL) % 定义 q 是一个 kL x 1 的变量
    
    % CVX 的 log 函数是凹的，q'*rate 是线性的
    % rate_per_user 是一个 1 x nUE 的行向量
    rate_per_user = qo.' * r_realUser; 
    
    % 最大化 比例公平(PF) 效用
    maximize( geo_mean(rate_per_user + 10^-6) ) 
    
    % 添加约束
    subject to
        sum(qo) == 1;
        qo >= 0;
cvx_end