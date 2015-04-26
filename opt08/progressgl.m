function progressgl(kk, xx, fval, gg, dval, dpen, problem, bound, step, step_grad, opt)
switch(opt.display)
 case 2
  fprintf('[%d] xx=[%g %g] fval=%.20f dval=%.20f gnorm=%g step=%g\n',...
          kk, xx(1), xx(2), fval, dval, dpen, step);
 case 1
  fprintf('.');
end


countzero = sum(sqrt(sum(reshape(xx,[problem.ns,problem.nc]).^2))==0);

if isfield(problem,'fid') && ~isempty(problem.fid)
  fprintf(problem.fid, '%d\t%f\t%d\t%d\t%.20f\t%.20f\t%.20f\t%.20f\t%.20f\t%d\n',...
          kk, ...
          cputime - problem.t0, ...
          problem.nfeval, ...
          bound, ...
          fval, ...
          dval, ...
          dpen, ...
          max(abs(gg)), ...
          step,...
          countzero);
end


if any(kk==problem.klog)
  file = sprintf('xx%d.txt',kk);
  fid=fopen(file,'w');
  fprintf(fid, '%.32f ', xx);
  fprintf(fid, '\n');
  fprintf(fid, '%.32f ', gg);
  fprintf(fid, '\n');
  fclose(fid);
  
  
  fid=fopen(file);
  xx1=readxg(file);
  [fval1,gg1,dval1]=feval(problem.obj, xx1, problem, opt);
  if ~isequal(gg1, gg)
    keyboard;
  end

end

  