function [Rw F r_propose] = calrateS(idx, Sv,Kv,M,NRF,Tx_CB,Rcron,PT,miu,weight)
% 根据Sv向量，从Tx_CB中为每个组贪婪选择波束
% 直到填满RF链，从而确定模拟预编码器F
if sum(Sv)==0
    F=zeros(M,NRF);
    Rw=0;
    r_propose=zeros(1,sum(Kv));
else 
    F=[];
    Tx_CBtemp=Tx_CB;
    for kj=1:length(Sv)
        [M,Nc]=size(Tx_CBtemp);
        if Sv(kj)~=0
            Ep=zeros(Nc,1);
            for knc=1:Nc
                Ep(knc)=real(Tx_CBtemp(:,knc)'*Rcron(:,:,kj)*Tx_CBtemp(:,knc));
            end
            [Eval Eindx]=sort(Ep,'descend');
            F = [F Tx_CBtemp(:,Eindx(1:Sv(kj)))];
            Tx_CBtemp=Tx_CBtemp(:,Eindx(Sv(kj)+1:Nc));
        end
    end
    % 射频链还有空余，再加
    if sum(Sv) < NRF
        [M,Nc]=size(Tx_CBtemp);
        Rcronsum=zeros(M,M);
        for kj=1:length(Sv)
            Rcronsum=Rcronsum+Sv(kj)*Rcron(:,:,kj);
        end;
        Ep=zeros(Nc,1);
        for knc=1:Nc
            Ep(knc)=real(Tx_CBtemp(:,knc)'*Rcronsum*Tx_CBtemp(:,knc));
        end;
        [Eval Eindx]=sort(Ep,'descend');
        F = [F Tx_CBtemp(:,Eindx(1:NRF-sum(Sv)))];
    end;
    r_group = CalrateG(M,Sv,PT,miu,Rcron,F);
    r_propose=zeros(1,sum(Kv));
    for kj=1:length(Sv)
        if kj==1
            r_propose(1:Kv(kj)) = ones(1,Kv(kj))*Sv(kj)*r_group(kj)/Kv(kj);
        else
            r_propose(sum(Kv(1:kj-1))+1:sum(Kv(1:kj)))=ones(1,Kv(kj))*Sv(kj)*r_group(kj)/Kv(kj);
        end
    end
    Rw=r_propose*weight;
end;

