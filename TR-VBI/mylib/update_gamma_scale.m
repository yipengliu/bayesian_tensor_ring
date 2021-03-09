function [gammas1,c_gammaN,d_gammaN]=update_gamma_scale(c_gamma0,d_gamma0,G,gammas,Vn,beta)
                 [rn, I1, r1] = size(G{1});
                 [r1, I2, r2] = size(G{2});
   
          c_gammaN = (0.25*(I1*rn+I2*r2)+ c_gamma0-1)*ones(r1,1);
             for r=1:r1
                 d_gammaN(r)=(beta{end}*gammas{end})'*(diag(G{1}(:,:,r)*(G{1}(:,:,r))') + squeeze(sum(Vn{1}(:,:,r),2)))+((beta{2}*gammas{2})')*(diag((squeeze(G{2}(r,:,:)))'*squeeze(G{2}(r,:,:))) + squeeze(sum(squeeze(Vn{2}(r,:,:)),1))');
                 d_gammaN(r)=d_gamma0+0.125*d_gammaN(r);
             end
            gammas1 = c_gammaN./d_gammaN';
end
  