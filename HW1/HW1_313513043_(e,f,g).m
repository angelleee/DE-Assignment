N1 = 1:1:256;
N2 = 2:1:256;
sigma_2 = 1;

e_E_alpha=1+(sigma_2./N1);
e_E_g=ones(1,256);
e_var_alpha=((4*sigma_2)./N2)+((2*sigma_2^2)./(N2.^2));
e_var_g=((4*sigma_2)./N2);

seed=101;
rng(seed);
x_4096=normrnd(1,sqrt(sigma_2),[1,4096]);
x_mean=zeros(1,4096);
x_var=zeros(1,4096);

for i=1:4096
    xm=sum(x_4096(1:i))/i;
    x_mean(i)=xm*xm;
    x_var(i)=(xm*xm-1)^2;
end

monte_e=(sum(x_mean)/4096)*ones(1,256);
monte_var=(sum(x_var)/4096)*ones(1,255);

figure
plot(N1,e_E_alpha,N1,e_E_g,N1,monte_e)
xlabel('N')
ylabel('expectation')
xlim([0 256])
ylim([0.9 2.1])
legend('E[a]','E[g(x)]','Monte Carlo')

figure
plot(N2,e_var_alpha,N2,e_var_g,N2,monte_var)
xlabel('N')
ylabel('variance')
legend('var[a]','var[g(x)]','Monte Carlo')

