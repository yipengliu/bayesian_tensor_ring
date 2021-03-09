function [G1,G2,V1,V2,R1,gammas1]=reduce_dimension_v1(G,V,gammas,deta)
           [rn, I1, r1] = size(G{1});
           [r1, I2, r2] = size(G{2});
           G_hat=reshape(TCP(G(1:2)),[rn*I1 r2*I2]);
           tau=deta*(norm(G_hat,'fro'));
           V_hat=reshape(TCP(V(1:2)),[rn*I1 r2*I2]);
           [G1,G2,R1,V1,V2]=Thr_SVD(G_hat,V_hat,tau);
           G1=reshape(G1,[rn,I1,R1]);
           G2=reshape(G2,[R1,I2,r2]);
           V1=reshape(V1,[rn,I1,R1]);
           V2=reshape(V2,[R1,I2,r2]); 
           
           gammas1=mink(gammas{1},R1);
            
            
%             N=length(G);
%            [rn, I1, r1] = size(G{1});
% %             ep=ep/sqrt(N-1);
%             conV{1}=diag(Unfold(G{1},3)*Unfold(G{1},3)');%obtain R{1}\times R_{1} conVar
%             comTol{1}=tau*norm(Unfold(G{1},3),'fro');
%             R1=sum(conV{1}>comTol{1});
%             if it>=2&& R1~=r1&&R1>1
%             indices = conV{1} > comTol{1};
%             gammas1 = gammas{1}(indices);
%        
%             G1 = G{1}(:,:,indices);
%             G2 =G{2}(indices,:,:);
%             V1 = V{1}(:,:,indices);
%             V2 =V{2}(indices,:,:);
%             else
%                 G1=G{1};
%                 G2=G{2};
%                 V1=V{1};
%                 V2=V{2};
%                 gammas1=gammas{1};
%        
%             end

            end

 
