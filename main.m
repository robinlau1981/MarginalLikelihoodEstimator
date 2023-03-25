clear
close all

%%---- Algorithm Initialization ----%%
dim=7; % dimension of the target likelihood function
N=2e4; % sample size
M=10;  % Number of mixing components in the initial proposal
df=5;  % degree of freedom of student's t distributions involved in all proposals
W_m=ones(1,M)*1/M; % Initial weight of each mixing component
Mu=unifrnd(-10,10,M,dim); % Initial mean of mixing components
Sigma=zeros(dim,dim,M); % Initial covariance of mixing components
for i=1:M
    Sigma(:,:,i)=1e3*eye(dim);
end

M0=M;
Mu0=Mu;
Sigma0=Sigma;
W_m0=W_m;

cor_thr=.8; % Threshold for merging
q_thr=1e1/N; % Threshold for comparison with biggest IS weight
gama=[.0001 .001 .01 .05:.05:1]; % the annealing schedule
thr_del=1e-6; % Threshold used for deleting mixing components

%%-------------------------------------%%
for j=1:length(gama)
    for i=1:N  % draw N random samples from the initial Student's t-mixture proposal
        r=W_m(1,1);
        rand_num=rand;
        for ii=1:M % Component index of the Gaussian mixture
            if rand_num<=r
                tao=gamrnd(df/2,1/(df/2)); % a random number drawn from a gamma distribution
                X(i,:)=mvnrnd(Mu(ii,:),Sigma(:,:,ii)/tao);
                break;
            else
                r=r+W_m(1,ii+1);
            end
        end

    end

    proposal=t_mixture_pdf(X,M,W_m,Mu,Sigma,df); % Proposal density values of these N samples;
    target=Target_likelihood_function(X).^gama(j).*t_mixture_pdf(X,M0,W_m0,Mu0,Sigma0,df).^(1-gama(j)); % Target density values of these N samples at this (jth) iteration;
    q(:,j)=target./proposal; % importance weights
    q(:,j)=q(:,j)/sum(q(:,j)); % normalized importance weights
    %---------------------------%
    W_m_minus=zeros(1,M);
    Mu_minus=zeros(M,dim);
    Sigma_minus=zeros(dim,dim,M);
    for ii=1:M
        rou(ii,:)=W_m(ii)*((t_pdf(X,Mu(ii,:),Sigma(:,:,ii),df))./proposal)';
        W_m_minus(ii)=rou(ii,:)*q(:,j);
        u(ii,:)=(df+dim)*ones(1,N)./(df*ones(1,N)+MahalanobisDist(X,Mu(ii,:),Sigma(:,:,ii)));
        Mu_minus(ii,:)=rou(ii,:).*u(ii,:)*(repmat(q(:,j),1,dim).*X)/(rou(ii,:).*u(ii,:)*q(:,j));
        for d=1:dim
            Sigma_minus(d,d,ii)=(q(:,j).*(X(:,d)-Mu_minus(ii,d)))'*((X(:,d)-Mu_minus(ii,d)).*(rou(ii,:).*u(ii,:))')/W_m_minus(ii);
        end
    end

    W_m_minus=W_m_minus/sum(W_m_minus);
    w_thr=1e-6/M; % used to delete component with too small component weight

    Cmp_ind=find(W_m_minus>w_thr);
    M=length(Cmp_ind);
    Mu=Mu_minus(Cmp_ind,:);
    Sigma=Sigma_minus(:,:,Cmp_ind);
    W_m=W_m_minus(Cmp_ind);
    W_m=W_m/sum(W_m);


    for i=1:N  % Particle index
        r=W_m(1,1);
        rand_num=rand;
        for ii=1:M % Component index of the Gaussian mixture
            if rand_num<=r
                tao=gamrnd(df/2,1/(df/2));% gamrnd(a,b)=gengamma(a,1/b);
                X(i,:)=mvnrnd(Mu(ii,:),Sigma(:,:,ii)/tao);

                break;
            else
                r=r+W_m(1,ii+1);
            end
        end
    end
    proposal=t_mixture_pdf(X,M,W_m,Mu,Sigma,df);%Gaussian_Mixture_pdf(X(i,:),M,W_m,Mu,Sigma);
    q(:,j)=Target_likelihood_function(X).^gama(j).*t_mixture_pdf(X,M0,W_m0,Mu0,Sigma0,df).^(1-gama(j))./proposal; %Target_pdf_v3(X,data).^gama(j)./proposal;
    Sigma_old2=Sigma;
    q(:,j)=q(:,j)/sum(q(:,j));

    %-----Component Merging-----%
    M2=M;
    Mu2=Mu;
    Sigma2=Sigma;
    W_m2=W_m;
    %tic;
    z=zeros(N,M2);
    for k=1:M
        z(:,k)=q(:,j).*W_m(k).*(t_pdf(X,Mu(k,:),Sigma(:,:,k),df))./t_mixture_pdf(X(i,:),M,W_m,Mu,Sigma,df); %z(:,k)=q(:,j).*mvnpdf(X,Mu2(k,:),Sigma2(:,:,k));
    end
    z_mean=mean(z);
    cor=zeros(M2,M2);
    for r=1:M2
        for c=r+1:M2
            cor(r,c)=(z(:,r)-z_mean(r)*ones(N,1))'*(z(:,c)-z_mean(c)*ones(N,1))/(norm(z(:,r)-z_mean(r)*ones(N,1))*norm(z(:,c)-z_mean(c)*ones(N,1)));
        end
    end

    m_ind=ones(1,M2);
    for r=1:M2
        if m_ind(r)==1
            for c=r+1:M2
                if m_ind(c)==1
                    if cor(r,c)>cor_thr
                        m_ind(c)=0;
                        Mu2(r,:)=W_m2([r c])*Mu2([r c],:)/sum(W_m2([r c]));
                        Sigma2(:,:,r)=(W_m2(r)*Sigma2(:,:,r)+W_m2(c)*Sigma2(:,:,c))/sum(W_m2([r c]));
                        W_m2(r)=sum(W_m2([r c]));
                    end
                end
            end
        end
    end

    M_remain_ind=find(m_ind==1);
    M=length(M_remain_ind);
    Mu=Mu2(M_remain_ind,:);
    Sigma=Sigma2(:,:,M_remain_ind);
    W_m=W_m2(M_remain_ind);

    %----------add component---------%
    j2=1;
    q_try=zeros(N,1);
    q_max=max(q(:,j))
    while (j2<10)&&(q_max>q_thr)
        if j2>1
            X=X_try;
            M=M_try;
            W_m=W_m_try;
            Sigma=Sigma_try;
            Mu=Mu_try;
            Sigma_old2=Sigma_try;
            q(:,j)=q_try;
        end

        [maxq(j),max_ind]=max(q(:,j));
        M_try=M+1;
        Mu_try=Mu;
        Sigma_try=Sigma;
        W_m_try=W_m;
        Mu_try(M_try,:)=X(max_ind(1),:);            % where? location of the maximum importance weight sample

        Sigma_try(:,:,M_try)=diag(ones(1,dim))*1;%1/9*diag((Mu_try(M_try,:)-Mu_min_dis).^2);

        W_m_try(M_try)=1/M_try;

        for i=1:N  % Particle index
            r=W_m_try(1,1);
            rand_num=rand;
            for ii=1:M_try % Component index of the Gaussian mixture
                if rand_num<=r
                    tao=gamrnd(df/2,1/(df/2));% gamrnd(a,b)=gengamma(a,1/b);
                    X_try(i,:)=mvnrnd(Mu_try(ii,:),Sigma_try(:,:,ii)/tao);
                    break;
                else
                    r=r+W_m_try(1,ii+1);
                end
            end
        end
        %X_try=X;
        proposal_try=t_mixture_pdf(X_try,M_try,W_m_try,Mu_try,Sigma_try,df);%Gaussian_Mixture_pdf(X(i,:),M,W_m,Mu,Sigma);
        q_try=Target_likelihood_function(X_try).^gama(j).*t_mixture_pdf(X_try,M0,W_m0,Mu0,Sigma0,df).^(1-gama(j))./proposal_try;% Target_pdf_v3(X,data).^gama(j)./proposal_try;
        q_try=q_try/sum(q_try);

        W_m_minus=zeros(1,M_try);
        Mu_minus=zeros(M_try,dim);
        Sigma_minus=zeros(dim,dim,M_try);
        for ii=1:M_try
            rou_try(ii,:)=W_m_try(ii)*((t_pdf(X_try,Mu_try(ii,:),Sigma_try(:,:,ii),df))./proposal_try)';
            W_m_minus(ii)=rou_try(ii,:)*q_try;
            u_try(ii,:)=(df+dim)*ones(1,N)./(df*ones(1,N)+MahalanobisDist(X_try,Mu_try(ii,:),Sigma_try(:,:,ii)));
            Mu_minus(ii,:)=rou_try(ii,:).*u_try(ii,:)*(repmat(q_try,1,dim).*X_try)/(rou_try(ii,:).*u_try(ii,:)*q_try);%W_m_minus(ii);
            %             for i=1:N
            %                 Sigma_minus(:,:,ii)=Sigma_minus(:,:,ii)+q(i,j)*(X(i,:)-Mu_minus(ii,:))'*(X(i,:)-Mu_minus(ii,:))*rou_try(ii,i)*u_try(ii,i)/W_m_minus(ii);
            %             end

            for d=1:dim
                Sigma_minus(d,d,ii)=(q_try.*(X_try(:,d)-Mu_minus(ii,d)))'*((X_try(:,d)-Mu_minus(ii,d)).*(rou_try(ii,:).*u_try(ii,:))')/W_m_minus(ii);
            end
        end

        W_m_minus=W_m_minus/sum(W_m_minus);
        w_thr=thr_del/M_try; % used to delete component with too small component weight
        Cmp_ind=find(W_m_minus>w_thr);
        M_try=length(Cmp_ind);
        Mu_try=Mu_minus(Cmp_ind,:);
        Sigma_try=Sigma_minus(:,:,Cmp_ind);
        W_m_try=W_m_minus(Cmp_ind);
        W_m_try=W_m_try/sum(W_m_try);

        for i=1:N  % Particle index
            r=W_m_try(1,1);
            rand_num=rand;
            for ii=1:M_try % Component index of the Gaussian mixture
                if rand_num<=r
                    tao=gamrnd(df/2,1/(df/2));% gamrnd(a,b)=gengamma(a,1/b);
                    X_try(i,:)=mvnrnd(Mu_try(ii,:),Sigma_try(:,:,ii)/tao);
                    break;
                else
                    r=r+W_m_try(1,ii+1);
                end
            end
        end

        proposal_try=t_mixture_pdf(X_try,M_try,W_m_try,Mu_try,Sigma_try,df);%Gaussian_Mixture_pdf(X(i,:),M,W_m,Mu,Sigma);
        q_try=Target_likelihood_function(X_try).^gama(j).*t_mixture_pdf(X_try,M0,W_m0,Mu0,Sigma0,df).^(1-gama(j))./proposal_try;% Target_pdf_v3(X,data).^gama(j)./proposal_try;
        q_try=q_try/sum(q_try);
        j2=j2+1
        q_max=max(q_try)
    end
    %Plot_GM(X,M,W_m,Mu',df/(df-2)*Sigma);grid on;
    %pause;

    ESS(j)=1/sum(q(:,j).^2);
    fprintf('v12_v5:j/Temp/ESS/M %i / %1.4f / %i /%i \r', j,gama(j),ESS(j), M);

end
%-------------------importance sampling for computing marginal likelihood-------------%
MarLik=0;
for i=1:N
    r=W_m(1,1);
    rand_num=rand;
    for ii=1:M
        if rand_num<=r
            tao=gamrnd(df/2,1/(df/2));% gamrnd(a,b)=gengamma(a,1/b);
            Y(i,:)=mvnrnd(Mu(ii,:),Sigma(:,:,ii)/tao);
            break;
        else
            r=r+W_m(1,ii+1);
        end
    end
end
proposal=t_mixture_pdf(Y,M,W_m,Mu,Sigma,df);
q=Target_likelihood_function(Y)./proposal;

MarLik=mean(q);  % The estimated marginal likelihood
MarLik           % Output the estimated marginal likelihood
s0_square=mean((q-MarLik).^2);
ci=3*sqrt(s0_square)/sqrt(N)  % credible interval defined by 3*standard error
q=q/sum(q);
ESS_final=1/sum(q.^2) % effective sample size

figure,plot(ESS);xlabel('Annealing time');ylabel('ESS');

% Scatter plot drawing %
for i = 1 : N
    u = rand; % uniform random number between 0 and 1
    qtempsum = 0;
    for ii = 1 : N
        qtempsum = qtempsum + q(ii);
        if qtempsum >= u
            Ytminus(i,:)=Y(ii,:);
            break;
        end
    end
end
Yt=Ytminus;
X=Yt;
figure,plotmatrix_v2(X);