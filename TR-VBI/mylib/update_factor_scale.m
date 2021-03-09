function [G_hat,GSigma,Vn,EGGTn]=update_factor_scale(Y,O,G,EGGT,gammas,tau)% update core data
        [rn, I1, r1] = size(G{1});
        N=length(size(O));
        Aw= kron(diag(gammas{end}),diag(gammas{1}));
        B=reshape(permute(TCP(G(2:end)),[2,3,1]),[numel(Y)/I1,rn*r1]);
        O = tens2mat(O,   1, 2:N);  % I_1 * (I_2...I_n)
        Y = tens2mat(Y,    1, 2:N);  % I_1 * (I_2...I_n)
        VV=reshape(permute(TCP(EGGT(2:end)),[2,3,1]),[numel(Y)/I1,rn*r1*rn*r1]);
        VV=reshape(permute(reshape(VV'*O',[rn rn r1 r1 I1]),[1 3 2 4 5]),[rn*r1 rn*r1 I1]);
  
        % the size of FlashT is R^(n-1)R(n)*I(n)
%         FslashT =tens2mat(Y.*O, 1,2:N)*reshape(permute(TCP(G(2:end)),[2,3,1]),[numel(Y)/I1,rn*r1]);
        for i=1:I1
            idx = find(O(i,:)==1);
            GSigma(:,:,i) = (tau * ZZ(:,:,i) + Aw )^(-1);    
            G_hat(:,i,:) =reshape(tau *(Y(i,idx)*B(idx,:))*GSigma(:,:,i),[rn,1,r1])  ;%i(n)*R(n-1)*R(n)
            temp(:,:,i)=GSigma(:,:,i)+vec(G_hat(:,i,:))'*vec(G_hat(:,i,:));
        end
       Vn=diagten(GSigma,r1);
       EGGTn=reshape(permute(reshape(tempt,[rn,r1,rn,r1,I1]),[1,3,5,2,4]),[rn*rn I1 r1*r1])
end





