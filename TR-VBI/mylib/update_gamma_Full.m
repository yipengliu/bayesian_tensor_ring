function [gammas1,c_gammaN,d_gammaN]=update_gamma_Full(c_gamma0,d_gamma0,G,gammas,Gsigma)
                 [rn, I1, r1] = size(G{1});
                 [r1, I2, r2] = size(G{2});
                 Vn{1}=reshape(Gsigma{1},[rn r1 rn r1]);
                 Vn{2}=reshape(Gsigma{2},[r1 r2 r1 r2]);
          c_gammaN = (0.5*(I1*rn+I2*r2)+ c_gamma0-1)*ones(r1,1);
             for r=1:r1
                 d_gammaN(r)=(gammas{end})'*(diag(G{1}(:,:,r)*(G{1}(:,:,r))') + diag(squeeze(Vn{1}(:,r,:,r))))+(gammas{2}')*(diag((squeeze(G{2}(r,:,:)))'*squeeze(G{2}(r,:,:))) + diag(squeeze(Vn{2}(r,:,r,:))));
                 d_gammaN(r)=d_gamma0+0.25*d_gammaN(r);
             end
            gammas1 = c_gammaN./d_gammaN';
end
  