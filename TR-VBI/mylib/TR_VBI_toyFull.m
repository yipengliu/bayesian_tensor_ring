function [model] = TR_VBI_toyFull(Y, init, maxRank , maxiters, tol)  % tensor ring with Bayesian variational inference 100% observation data
%  INPUTS
%     Y              - Input tensor
%     'O'          - Binary (0-1) tensor indicating missing entries
%                      (0: missing; 1: observed)
%     'init'         - Initialization method
%                     - 'ml'  : SVD initilization (default)
%                     - 'rand': Random matrices
%     'maxRank'      - The initialization of rank (larger than true rank)
%     'maxiters'     - max number of iterations (default: 100)
%     'tol'          - lower band change tolerance for convergence dection
%                      (default: 1e-4)
%   OUTPUTS
%      B         - recovered tensor
%      LB        - the value of ELBO
%%
dimY = size(Y);
N = ndims(Y);
R=maxRank;
noise='on';

%% Initialization
nObs = prod(dimY);
oR = nObs/prod(dimY);
 scrsz = get(0,'ScreenSize');
 h1 = figure('Position',[scrsz(3)*0.2 scrsz(4)*0.3 scrsz(3)*0.6 scrsz(4)*0.4]); 
 figure(h1);
for i=1:N
c_gamma0{i}     = 1e-6;
d_gamma0{i}     = 1e-6;
end
if  strcmp(noise,'on')
    a_tau0     = 1e-6;
    b_tau0      = 1e-6;
else
    a_tau0     = 1e-1;
    b_tau0      = 1e-6;
end


for i=1:N
gammas{i} = ones(R(i),1);
end
dscale = std(Y(:))/N;
Y = Y./dscale;
tau = 1e5;


%% facror prior initialization

switch init,
    case 'ml'    
    [G,GSigma] = TR_Initialization_ml(Y,  R);
    case 'rand' 
    [G,GSigma] = TR_Initialization_rand(Y,  R);
end
% --------- E(aa') = cov(a,a) + E(a)E(a')----------------
EGGT = cell(N,1);
    EGGT{1} = reshape(permute(reshape(unfold(G{1},2)'*unfold(G{1},2) +sum(GSigma{1},3),[R(N) R(1) R(N) R(1)]),[1,3,2,4]),[R(N)*R(N) R(1)*R(1)]);%% R(n-1)R(n)*R(n-1)R(n)

for n=2:N
    EGGT{n} = reshape(permute(reshape(unfold(G{n},2)'*unfold(G{n},2) +sum(GSigma{n},3),[R(n-1) R(n) R(n-1) R(n)]),[1,3,2,4]),[R(n-1)*R(n-1) R(n)*R(n)]);%% R(n-1)R(n)*R(n-1)R(n)
end
LB = 0;%lower bound

%% Model learning
for it=1:maxiters
    %% Update factor cores
     for n=1:N
     [G{n},GSigma{n},EGGT{n}]=update_factor_Full(TensPermute(Y, n), G([n:N, 1:n-1]),EGGT([n:N, 1:n-1]),gammas([n:N, 1:n-1]),tau);
     end
         %% smooth contraint
%     for n=1:N-1
%         [G{n}]=back_process( G{n});
%     end   
    %% Update hyperparameters gamma
    for n=1:N
        [gammas{n},c_gammaN{n},d_gammaN{n}]=update_gamma_Full(c_gamma0{n},d_gamma0{n},G([n:N, 1:n-1]),gammas([n:N, 1:n-1]),GSigma([n:N, 1:n-1]));
    end

    %% update noise tau
      B=Ui2U(G);
       EX2=eye(R(N)*R(N),R(1)*R(1));
        for i=2:length(G)
            EX2=EX2*EGGT{i};
        end
          EX2 = trace(EX2);
   err(it) = abs(Y(:)'*Y(:) - 2*Y(:)'*B(:) + EX2);
   sqrt(err(it))/(norm(Y(:),'fro')^2)
   
    if  strcmp(noise,'on')
    a_tauN = a_tau0 + 0.5*nObs;
    b_tauN = b_tau0 + 0.5*err(it);
    else
        a_tauN = a_tau0;
        b_tauN = b_tau0;
    end


    tau = a_tauN/b_tauN;
     
  %% Lower bound
%     temp1 = 0.5*nObs*(psi(a_tauN)-safelog(b_tauN)) -0.5*N*nObs*safelog(2*pi)  - 0.5*(a_tauN/b_tauN)*err(it);
%     temp2 =0;
%     for n=1:N
%         temp2=temp2+caculate_tempt2(c_gammaN{n},d_gammaN{n},c_gamma0{n},d_gamma0{n},G([n:N, 1:n-1]),gammas([n:N, 1:n-1]),V([n:N, 1:n-1]));
%     end
%     temp3 = -safelog(gamma(a_tau0)) + a_tau0*safelog(b_tau0) + (a_tau0-1)*(psi(a_tauN)-safelog(b_tauN)) - b_tau0*(a_tauN/b_tauN);
%     temp4=0;
%     for n=1:N
%     temp4=temp4+caculate_tempt4(G([n:N, 1:n-1]),GSigma([n:N, 1:n-1]));
%     end
%     temp5=0;
%     for n=1:N
%     temp5=temp5+ 2*sum(safelog(gamma(c_gammaN{n})) - (c_gammaN{n}-1).*psi(c_gammaN{n}) -safelog(d_gammaN{n}') + c_gammaN{n});
%     end
%     temp6 = safelog(gamma(a_tauN)) - (a_tauN-1)*psi(a_tauN) -safelog(b_tauN) + a_tauN;
% 
%     LB(it) = temp1 + temp2 + temp3 + temp4 + temp5 + temp6 ;
    %% reduce irrelevant dimensions
   deta=(sqrt(err(it))/norm(Y(:),'fro'))/sqrt(N);
    for n=1:N-1
      [G{n},G{n+1},EGGT{n},EGGT{n+1},gammas{n},r]=reduce_dimension_Full(G([n:N, 1:n-1]),EGGT([n:N, 1:n-1]),gammas([n:N, 1:n-1]),it,deta);
        R(n)=r;
    end
     [G{N},G{1},EGGT{N},EGGT{1},gammas{N},r]=reduce_dimension_Full(G([N, 1:N-1]),EGGT([N, 1:N-1]),gammas([N, 1:N-1]),it,deta);
         R(N)=r;
    %% plot 
%      % Notice:the number of gammas is corresponding to  the data dimension
        set(0,'CurrentFigure',h1);
        subplot(3,4,1); bar(gammas{1}); title('Posterior mean of \lambda_{1}'); xlabel('Latent components'); ylabel(''); axis tight;
        subplot(3,4,2); bar(gammas{2}); title('Posterior mean of \lambda_{2}'); xlabel('Latent components'); ylabel(''); axis tight;
        subplot(3,4,3); bar(gammas{3}); title('Posterior mean of \lambda_{3}'); xlabel('Latent components'); ylabel(''); axis tight;
        subplot(3,4,4); bar(gammas{4}); title('Posterior mean of \lambda_{4}'); xlabel('Latent components'); ylabel(''); axis tight;
%         subplot(3,4,5); bar(gammas{5}); title('Posterior mean of \lambda_{5}'); xlabel('Latent components'); ylabel(''); axis tight;
% %         subplot(3,4,6); bar(gammas{6}); title('Posterior mean of \lambda_{6}'); xlabel('Latent components'); ylabel(''); axis tight;
% %         subplot(3,4,7); bar(gammas{7}); title('Posterior mean of \lambda_{7}'); xlabel('Latent components'); ylabel(''); axis tight;
% %         subplot(3,4,8); bar(gammas{8}); title('Posterior mean of \lambda_{8}'); xlabel('Latent components'); ylabel(''); axis tight;
% %         subplot(3,4,9); bar(gammas{9}); title('Posterior mean of \lambda_{9}'); xlabel('Latent components'); ylabel(''); axis tight;
% %         subplot(3,4,10); plot(LB, '-r.','LineWidth',1.5,'MarkerSize',10 ); title('Lower bound'); xlabel('Iteration');  grid on;
        subplot(3,4,11); plot(err, '-b.','LineWidth',1.5,'MarkerSize',10 ); title('err'); xlabel('Iteration');  grid on;
% %         subplot(3,4,12);imshow(uint8(reshape(B,[256,256,3])));  xlabel('Iteration');grid on;
        set(findall(h1,'type','text'),'fontSize',12);
        drawnow;
        if err<tol
            break;
        end 

end
        %output
        SNR = 10*log10(var(B(:))*tau);
        model.X = B;
model.SNR = SNR;
model.TrueRank = R;


