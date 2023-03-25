function pdf =t_pdf(X,Mu,Sigma,df) %Ref: Robust mixture modelling using the t distribution;ML estimation of the t distribution using EM and its extensions,ECM and ECME
% X and Mu are row vectors
[n,p]=size(X);
pdf=zeros(n,1);
for k=1:n
pdf(k)=gamma((df+p)/2)*(det(Sigma))^(-.5)/((pi*df)^(p/2)*gamma(df/2)*(1+((X(k,:)-Mu)*inv(Sigma)*(X(k,:)-Mu)')/df)^((df+p)/2));
end