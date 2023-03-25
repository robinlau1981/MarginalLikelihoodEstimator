function pdf =Target_pdf_out_prod_density_v3(X);  % Target density being a outer product of 7 univariate densities % Ref: Posterior-guided importance sampling for calculating marginal likelihoods
[n,p]=size(X);
pdf=zeros(n,1);
for k=1:n
pdf(k)=1;
pdf(k)=pdf(k)*(3/5*gampdf(X(k,1)+10,2,3)+2/5*gampdf(10-X(k,1),2,5));
pdf(k)=pdf(k)*(3/4*dsn(X(k,2),3,1,5)+1/4*dsn(X(k,2),-3,3,-6));
pdf(k)=pdf(k)*mvtpdf(X(k,3),9,4);
pdf(k)=pdf(k)*(1/2*betapdf(X(k,4)+3,3,3)+1/2*mvnpdf(X(k,4),0,1));
pdf(k)=pdf(k)*(1/2*exppdf(X(k,5),1)+1/2*exppdf(-X(k,5),1));
pdf(k)=pdf(k)*dsn(X(k,6),0,8,-3);
pdf(k)=pdf(k)*(1/8*mvnpdf(X(k,7),-10,.1)+1/4*mvnpdf(X(k,7),0,.15)+5/8*mvnpdf(X(k,7),7,.2));
end
