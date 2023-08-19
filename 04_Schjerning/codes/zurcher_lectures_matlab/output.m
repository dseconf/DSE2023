classdef output
  methods (Static)
      function [results] = savemc(result_i, results, iMC);
          fields = fieldnames(result_i); 
          for j=1:numel(fields);
            if isnumeric(result_i.(char(fields(j)))) || islogical(result_i.(char(fields(j))))
              pvec=struct2vec(result_i,fields(j));
              results.(char(fields(j)))(1:numel(pvec), iMC)=pvec;
            end
          end
      end  

  function [] = table_mc(pnames, mp, results)
      fprintf('\nTable 1: Parameter estimates\n');
      for ir=1:numel(results); 
        r=results{ir};
        nMC=numel(r.n);
        fprintf('\nMethod: %s\n', r.title); 
        fprintf('------------------------------------------------------------------------------------\n');
        fprintf('%-14s %14s %14s %14s %14s\n', '', 'True value','Estimate', 'MCSD' ,'Bias');
        fprintf('------------------------------------------------------------------------------------\n');
        for ip=1:numel(pnames);
          trueval=struct2vec(mp,{pnames{ip}});
          theta_i=struct2vec(r,{pnames{ip}});
          theta_i=reshape(theta_i, numel(trueval), nMC);
          theta_i=theta_i(:,r.converged==1); 
          for j=1:numel(trueval)
            fprintf('%-14s %14.3f %14.3f %14.3f %14.4f \n',  pnames{ip}, trueval(j), mean(theta_i(j,:)), std(theta_i(j,:)) , mean(theta_i(j,:)- trueval(j))); 
          end
        end
      fprintf('%14s %14.3f \n\n', 'log-likelihood',  mean(r.llval(r.converged==1))); 
      end
      fprintf('------------------------------------------------------------------------------------\n');
    end

    function [] = table_np(results);
      fprintf('\nTable 2: Numerical performance\n');      
      for ir=1:numel(results); 
        r=results{ir};
        fprintf('\nMethod: %s\n', r.title); 
        nMC=numel(r.RC(1,:)); 
        fprintf('---------------------------------------------------------------------------------------\n');
        fprintf('%14s %14s  %14s %14s %14s\n', '', 'Runs Converged', 'CPU Time' ,'# of Major','# of Func.');
        fprintf('%14s %14s  %14s %14s %14s\n', '', sprintf(' (out of %g)',nMC),'(in sec.)','Iter','Eval.');
        fprintf('---------------------------------------------------------------------------------------\n');

        fprintf('%14s %14.3f %14.3f %14.1f %14.1f \n', '', ...
        sum(r.converged==1) , ...
        mean(r.cputime(r.converged==1)) , ...
        mean(r.MajorIter(r.converged==1)), ...
        mean(r.funcCount(r.converged==1))); 
      end
        fprintf('---------------------------------------------------------------------------------------\n');
    end
  function [p0] = estimates(p0, pnames, pvec, Avar, truevalue);
      p0=vec2struct(pvec, pnames, p0);
      se=sqrt(diag(Avar));

      k=numel(pnames);
      if nargin==5
          fprintf('    %-16s %4.4s %13s %13s %13s %13s\n','Param.',' ','True Value','Estimates','s.e.','t-stat');
      else
          fprintf('    %-16s %4.4s %13s %13s %13s\n','Param.',' ','Estimates','s.e.','t-stat');
      end
      fprintf('----------------------------------------------------------------------------------------\n');
       i=0;
       for iP=1:numel(pnames);
              for j=1:numel(p0.(char(pnames(iP))))
                  i=i+1;
                  myStr = ' ';
                  if numel(p0.(char(pnames(iP))))>1; 
                      myStr = sprintf('(%i)',j);
                  end;
                  if nargin==5
                      fprintf('    %-16s %4.4s %13.4f %13.4f %13.4f %13.4f\n', char(pnames(iP)), myStr, truevalue(i), p0.(char(pnames(iP)))(j), se(i), p0.(char(pnames(iP)))(j)/se(i));
                  else
                      fprintf('    %-16s %4.4s %13.4f %13.4f %13.4f\n', char(pnames(iP)), myStr, p0.(char(pnames(iP)))(j), se(i), p0.(char(pnames(iP)))(j)/se(i));
                  end
              end
      end
      fprintf('----------------------------------------------------------------------------------------\n');
    end % end of estimates
  end % end of methods
end % end of output class