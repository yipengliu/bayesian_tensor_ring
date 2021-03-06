function [gammas1]=update_gamma_v2(c_gamma,G,gammas,Gsigma)
              [rn, I1, r1] = size(G{1});
              [r1, I2, r2] = size(G{2});
              Gsigma{1}=reshape(Gsigma{1},[rn,r1,rn,r1,I1]);
              Gsigma{2}=reshape(Gsigma{2},[r1,r2,r1,r2,I2]);
            B = (I1*rn+I2*r2);
             for r=1:r1
                 A=(gammas{end})'*(diag(G{1}(:,:,r)*(G{1}(:,:,r))' + sum(squeeze(Gsigma{1}(:,r,:,r,:)),3)))+(gammas{2}')*(diag((squeeze(G{2}(r,:,:)))'*squeeze(G{2}(r,:,:))) + diag(sum(squeeze(Gsigma{2}(r,:,r,:,:)),3)));
                 d_gammaN(r)=(B+sqrt(B^2+2*c_gamma*A))/A;
             end
      gammas1=d_gammaN';
end
  

  