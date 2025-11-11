function y=CalUq(q,rate)
    y=-sum(log(q'*rate+10^-6));
end