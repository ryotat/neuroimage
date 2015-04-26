lambda=exp(linspace(log(100),log(0.01),20))
S=load('bcicomp2003iv/result_1stalpbet_ds.mat');
for ii=1:length(lambda)
time(ii)=S.memo(ii).C.status.time;
nm(ii)=normc(S.memo(ii).C, 'ds', 'v', S.v, 'div', 0);
end

S1=load('results/bcicomp2003iv/lrds_1stalpbet.mat')
for ii=1:length(lambda)
time1(ii)=S1.memo(ii).C.status.time(end);
nm1(ii)=normc(S1.memo(ii).C, 'ds', 'v', S1.memo(1).C.v, 'div', 0);
end
hold on;
plot(log(nm1), time1, 'r--', 'linewidth',2)
set(gca,'yscale','log')
grid on;
set(gca,'fontsize',16);
legend('ProjGrad','DAL')
xlabel('||W||_{DS}')
ylabel('computation time (s)')

