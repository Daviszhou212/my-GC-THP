function r_propose = CalrateG(M,Sv,PT,miu,R,F)
% 求解不动点方程，来得到速率DE
S=sum(Sv);
J = length(Sv);
% \tilde{R}
for kj=1:J
    R_effct(:,:,kj) = F * (F') * R(:,:,kj) * F * (F');
end
eps = zeros(J,1);
epspre=zeros(J,1);
left = 1;
% 迭代20次
while left<=20
    Tinv=eye(M);
    for kj=1:J
        Tinv=Tinv+Sv(kj)*R_effct(:,:,kj)/(miu+eps(kj))/M;
    end
    T=Tinv^-1;
    for kj=1:J
        eps(kj)= real(1/M * trace(R_effct(:,:,kj)*T));
    end
    epstlt(left)=eps(1);
    left = left+1;
%     if norm(epspre-eps)<10^-6
%         break;
%     end;
    epspre=eps;
end
% plot(epstlt);
for kj=1:J
    r_propose(kj)= log2(1+PT*M*real(eps(kj))/S);
end

% % % r_propose(r_propose==max(r_propose)) = min(r_propose);
