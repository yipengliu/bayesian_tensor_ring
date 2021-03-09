function [G1,G2,EGGT1,EGGT2,gammas1,R1]=reduce_dimension_Full(G,EGGT,gammas,it,deta)
            N=length(G);
           [rn, I1, r1] = size(G{1});
           [r1, I2, r2] = size(G{2});
%             ep=ep/sqrt(N-1);
            C{1}=Unfold(G{1},3);
            C{2}=Unfold(G{2},1);
            A=cell2mat(C);
            conV{1}=diag(A*A');%obtain R{1}\times R_{1} conVar
%             comTol{1}=find_tau(conV{1})
            comTol{1}=eps*(norm(A,'fro'));
            R1=sum(conV{1}>comTol{1});
            if it>=1&& R1~=r1&&R1>1
            indices = conV{1} > comTol{1};
            gammas1 = gammas{1}(indices);
            temp = ones(r1,r1);
            temp(indices,indices) = 0;
            temp = temp(:);
            G1 = G{1}(:,:,indices);
            G2 =G{2}(indices,:,:);
            EGGT1 = EGGT{1}(:,temp==0);
            EGGT2 = EGGT{2}(temp==0,:);
            else
                G1=G{1};
                G2=G{2};
                EGGT1=EGGT{1};
                EGGT2=EGGT{2};
                gammas1=gammas{1};
       
            end

            end

 
