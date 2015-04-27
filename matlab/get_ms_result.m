function ixcv=get_ms_result(acc_cv, ixrep)

ixcv=zeros(1,length(ixrep));
for kk=1:length(ixrep)
  M=ixrep(kk);
  [mm,ix]=max(acc_cv(:,M));
  ixcv(kk) =medin(find(acc_cv(:,M)==mm));
end
