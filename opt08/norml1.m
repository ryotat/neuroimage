function nm=norml1(xx, problem)

if problem.hasbias
  xx = xx(1:end-1);
end

nm=sum(abs(xx));
