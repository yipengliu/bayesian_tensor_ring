function [G_hat,GSigma,EGGT_hat]=update_factor_Full(Y,G,EGGT,gammas,tau)% 100% observation
        [rn, I1, r1] = size(G{1});
        [r1, I2, r2] = size(G{2});
        % compute E(G_{\n}^{T} G_{\n})
        ENGGT=eye(r1*r1,r2*r2);
        for i=2:length(G)
            ENGGT=ENGGT*EGGT{i};
        end
        ENGGT=reshape(permute(reshape(ENGGT',[rn rn r1 r1]),[1,3,2,4]),[rn*r1 rn*r1]);
        Aw= kron(diag(gammas{end}),diag(gammas{1}));
        
          % compute E(G_{\n})
        FslashT =tens2mat(Y, 1,2:length(G))*reshape(permute(TCP(G(2:end)),[2,3,1]),[numel(Y)/I1,rn*r1]);%% I1*Rn-1Rn
        GSigma=(tau*ENGGT+Aw)^(-1);%%Rn-1Rn*Rn-1Rn
        G_hat =permute(reshape(tau *(FslashT)*GSigma,[I1,rn,r1]),[2,1,3])  ;%
        EGGT_hat=reshape(permute(reshape(unfold(G_hat,2)'*unfold(G_hat,2) +GSigma,[rn r1 rn r1]),[1,3,2,4]),[rn*rn r1*r1]);
end





