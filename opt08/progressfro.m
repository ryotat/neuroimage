function progressfro(kk, xx, fval, gg, dval, dpen, problem, bound, step, step_grad, opt)
 
gnorm = norm(gg);


switch(opt.display)
 case 2
  if problem.hasbias
    ix2 = length(xx);
  else
    ix2 = 2;
  end
  
  fprintf('[%d] xx=[%g %g] fval=%g dval=%g gnorm=%g step=%g\n',...
          kk, xx(1), xx(ix2), fval, dval, gnorm, step);
 case 1
  fprintf('.');
end


if isfield(problem,'fid') && ~isempty(problem.fid)
  fprintf(problem.fid, '%d\t%f\t%d\t%d\t%.8f\t%.8f\t%.8f\t%.8f\t%f\t%f\n',...
          kk, ...
          cputime - problem.t0, ...
          problem.nfeval, ...
          bound, ...
          fval, ...
          dval, ...
          dpen, ...
          gnorm, ...
          step, step_grad);
end

