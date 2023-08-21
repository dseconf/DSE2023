classdef estim
	properties (Constant)
    analytical=true;
    scaling=1;
    oldcarage=10;
    deriv_tol=1e-8;      % Tolerance on norm of error of derivatives; 
    momlist={'old_cars_tau_j', 'mean_age_tau_j', 'new_cars_tau_j'};

    theta={ 'u_a','u_a_sq', 'u_even', 'u_0','u_og','mum', 'phi', 'transcost','ptranscost','psych_transcost','psych_transcost_nocar', ...
              'pnew_notax','pscrap_notax','nocarsc','ptc_sale','tc_sale','tc_sale_age','tc_sale_even','sigma_s'};
    alpha={'acc_0','acc_a','acc_even'};
	end

	methods (Static)
    function [mp_hat, theta_hat, Avar_hat, sol_hat, pvec0] = mle(mp, s, dta, algorithm, MaxIter)
      % estim.mle: Estimates models using mle and produces estimation output. 
      %
      % syntax:  estim.mle(mp, s, dta, algorithm, MaxIter)
      % OUTPUT: 
      %    mp_hat structure with model parameters updated with estimated parameters
      %    theta_hat      parameter vector
      %    Avar_hat       variance-covariance matrix of theta
      %    sol_hat        model solution at mp_hat
      %    pvec0          sequence of parameter vectors for each algorithm used 
      %                   (first column contains starting values, last column contains theta_hat)
      %
      % INPUT: 
      %    mp     structure with model parameters (see trmodel.setparams)
      %    s      structure with indexes for model (see trmodel.index)
      %    dta    data set (matlab table)
      %
      % optional inputs: 
      %   algorithm (cell array) with element:  'bhhh', 'bfgs', 'newton_fd', 'nelder-mead'
      %
      %   MaxIter (vector) with elements:        max number of iterations for each algorithm
      %
      %   If algorithm and MaxIter have multiple elements, estim.mle use a sequence of 
      %   algorithms in the specified order
      %   If one algorithm exits with positive exitflag. It does not try the next one, 
      %   but terminates and reports results. 
      %
      % Example:
      %   estim.mle(mp, s, dta, {'nelder-mead', 'bhhh', 'bfgs'}, [10,20,30]
      %   starts with max 10 iterations of nelder-mead, then max 10 iterations using bhhh 
      %   and then 30 iterations of bfgs. 

      t0=tic;
      if nargin<4
        algorithm={'bhhh'}; % can be 'bhhh', 'bfgs', 'newton_fd', 'nelder-mead'
      end

      if nargin<5
        MaxIter=10*ones(numel(algorithm),1);
        MaxIter(end)=100;
      end

      % get parameter with starting values correspoding to values in mp.p0 and mp.pnames
      [pvec0,mp] = estim.mp2pvec(mp,mp.pnames); % vector of parameter starting values

      % evaluate likelihood
      f_mle= @(pvec) estim.nloglike(mp, s, dta, pvec);

      fprintf('\nMaximize likelihood using algorithms\n')
      for i=1:numel(algorithm) 
        fprintf(' Algorithm %d: %15s,    Max number of iterations: %d\n', i, algorithm{i}, MaxIter(i))
      end

      for i=1:numel(algorithm)
        opt_default=optimset('MaxIter',MaxIter(i), 'Display','iter', 'GradObj','on', 'TolFun',1E-6,'TolX',1E-6);
        opt=opt_default;
        switch algorithm{i}
          case 'nelder-mead'
          case 'bhhh'
            opt=optimset(opt_default, 'Hessian','on', 'Algorithm', 'trust-region');
          case 'newton_fd'
            opt=optimset(opt_default, 'Algorithm', 'trust-region');
          case 'bfgs'
            opt=optimset(opt_default, 'Algorithm', 'quasi-newton');
          case 'steepest descent'
            opt=optimset(opt_default, 'Algorithm', 'trust-region', 'HessUpdate', 'steepdesc', ...
                'TolFun', 1e-10, 'TolX', 1e-10); 
          otherwise
            fprintf('Requested algorithm, %s, not found: switching to default')
            opt=opt_default;
        end

        fprintf('\nStarting %s\n - algorithm %d out of %d \n', algorithm{i}, i, numel(algorithm) ); 

        t1=tic;
        if strcmp(algorithm{i}, 'nelder-mead')
          [pvec0(:,i+1), fval,exitflag] = fminsearch(f_mle,pvec0(:,i), opt); 
        else
          [pvec0(:,i+1), fval,exitflag] = fminunc(f_mle, pvec0(:,i), opt);
        end
        fprintf('Time spend using %s: %g (seconds) \n', algorithm{i}, toc(t1)); 

        if exitflag == 1
          fprintf('CONVERGENCE: Algorithm %d: %15s found solution. (flag=%d)\n', i, algorithm{i}, exitflag);
          break
        end
 
     end
     theta_hat=pvec0(:,end);

     [mp_hat] = estim.update_mp(theta_hat, mp);
     
      mp_hat.sp = trmodel.update_structural_par(mp_hat);
     
      % compute variance covariance matrix
      [f_hat,g_hat,h_hat]=f_mle(theta_hat);     
      Avar_hat=inv(h_hat)/sum(dta.count);
      g_hat=reshape(g_hat, numel(g_hat),1);
      estim.output(mp, theta_hat, Avar_hat);
      
      % final solution at mp_hat
      global p0
      [sol_hat]=equilibrium.solve(mp_hat, s, p0);

      fprintf('%30s: %g\n',       'Mean log likelihood ',   -f_hat); 
      fprintf('%30s: %g\n',       'Norm of gradient (g)',   max(abs(g_hat))); 
      fprintf('%30s: %g\n',       'g''inv(H)g ',            g_hat'*inv(h_hat)*g_hat); 
      fprintf('%30s: %d\n',       'Number of cells ',       numel(dta.count)); 
      fprintf('%30s: %d\n',       'Number of obs.',         sum(dta.count)); 
      fprintf('%30s: %g (min)\n', 'Total CPU-time ',        toc(t0)/60); 
      
    end

    function [] = output_err(mp, nm, an);
      pnames=mp.pnames;
      rel=nm./an-1;
      err=nm-an;
      k=numel(pnames);
      fprintf('\n');
      fprintf('\nDerivative Check: \n\n');
      fprintf('    %-16s %4.4s %13s %13s %13s %13s %25s %10s\n','Param.',' ','Numerical','Analytical','num-an', 'num./an-1', 'Consumer types', 'Car types');
      fprintf('---------------------------------------------------------------------------------------------------------------\n');
      i=0;
      for iP=1:numel(pnames);
        [nr,nc]=size(mp.p0.(char(pnames(iP)))); 
        if (nr > 1)
          hh_type='';
        else
          hh_type='all';
        end
        if (nc > 1)
          car_type='';
        else
          car_type='all';
        end
        for j=1:numel(mp.p0.(char(pnames(iP))))
          [r,c]=ind2sub([nr,nc],j);
          if (nr == mp.ntypes)
             hh_type=char(mp.lbl_types{r});
          end
          if (nc == mp.ncartypes)
             car_type=char(mp.lbl_cartypes{c});
          end
          i=i+1;
          myStr = ' ';
          if numel(mp.(char(pnames(iP))))>1; 
            myStr = sprintf('(%i)',j);
          end;
          fprintf('    %-16s %4.4s %13.4f %13.4f %13.4g %13.4g %25s %10s\n', char(pnames(iP)), myStr, nm(i), an(i), err(i), rel(i), hh_type, car_type);
        end
      end
      fprintf('----------------------------------------------------------------------------------------\n');
    end % end of estim.output

    function [mp] = output(mp, theta, Avar);
      pnames=mp.pnames;
      se=sqrt(diag(Avar));

      k=numel(pnames);
      fprintf('\n');
      fprintf('\nEstimation output: \n\n');
      fprintf('    %-16s %4.4s %13s %13s %13s %25s %10s\n','Param.',' ','Estimates','s.e.','t-stat','Consumer type','Car type');
      fprintf('-------------------------------------------------------------------------------------------------------------\n');
      i=0;
      for iP=1:numel(pnames);
        [nr,nc]=size(mp.p0.(char(pnames(iP)))); 
        if (nr > 1)
          hh_type='';
        else
          hh_type='all';
        end
        if (nc > 1)
          car_type='';
        else
          car_type='all';
        end
        for j=1:numel(mp.p0.(char(pnames(iP))))
          [r,c]=ind2sub([nr,nc],j);
          if (nr == mp.ntypes)
             hh_type=char(mp.lbl_types{r});
          end
          if (nc == mp.ncartypes)
             car_type=char(mp.lbl_cartypes{c});
          end
          i=i+1;
          myStr = ' ';
          if numel(mp.(char(pnames(iP))))>1; 
            myStr = sprintf('(%i)',j);
          end;
          fprintf('    %-16s %4.4s %13.4f %13.4f %13.4f %25s %10s\n', char(pnames(iP)), myStr, theta(i), se(i), theta(i)/se(i),hh_type,car_type);
        end
      end
      fprintf('----------------------------------------------------------------------------------------\n');
    end % end of estim.output

    function [mp] = outputt(mp, theta, theta_true, Avar);
      pnames=mp.pnames;
      se=sqrt(diag(Avar));

      k=numel(pnames);
      fprintf('\n');
      fprintf('\nEstimation output: \n\n');
      fprintf('    %-16s %4.4s %13s %13s %13s %13s\n','Param.',' ','True values','Estimates','s.e.','t-stat');
      fprintf('----------------------------------------------------------------------------------------\n');
      i=0;
      for iP=1:numel(pnames);
        for j=1:numel(mp.p0.(char(pnames(iP))))
          i=i+1;
          myStr = ' ';
          if numel(mp.(char(pnames(iP))))>1; 
            myStr = sprintf('(%i)',j);
          end;
          fprintf('    %-16s %4.4s %13.4f %13.4f %13.4f %13.4f\n', char(pnames(iP)), myStr, theta_true(i), theta(i), se(i), theta(i)/se(i));
        end
      end
      fprintf('----------------------------------------------------------------------------------------\n');
    end % end of estim.output

    function print_parameters(mp,pnames)
        % This produces the chunk of text that can be copied in as the
        % default parameter values .
        s2p = mp;
        s2pStr = 'mp';
        if numel(pnames)==0;
            fprintf('No parameters to print\n'); 
            return
        end

        % fprintf('--- Printing parameters --- \n');
        for i=1:numel(pnames);
            if numel(s2p.(char(pnames(i))))==1 && ~strcmp(pnames(i),'pnames');
                try
                    fprintf('%s.%-20s = %15.8g;\n',s2pStr,char(pnames(i)),s2p.(char(pnames(i))));
                catch err
                    keyboard;
                end;
            elseif ~iscell(s2p.(char(pnames(i)))) && ~strcmp(pnames(i),'pnames');
                for j=1:numel(s2p.(char(pnames(i))));
                    myStr = sprintf('%s(%i)',char(pnames(i)),j);
                    fprintf('%s.%-20s = %15.8f;\n',s2pStr,myStr,s2p.(char(pnames(i)))(j));
                end;
            elseif strcmp(pnames(i),'pnames');
                for j=1:numel(s2p.(char(pnames(i))));
                    myStr = sprintf('%s(%i)',char(pnames(i)),j);
                    fprintf('%s.%-20s = ''%s'';\n',s2pStr,myStr,s2p.(char(pnames(i))){j});
                end;
            else
                for j=1:numel(s2p.(char(pnames(i))));
                    myStr = sprintf('%s{%i}',char(pnames(i)),j);
                    fprintf('%s.%-16s = ''%s'';\n',s2pStr,myStr,char(s2p.(char(pnames(i)))(j)));
                end;
            end; % huge if-else
        end; % for i=1:numel(pnames)
    end % print_parameters     

    function [pvec0,mp] = mp2pvec(mp,pnames)
      %PARAM2VEC: Procedure that extract vector from a structure with matrix fields. 
      %
      % SYNTAX:
      %   pvec0 = mp2pvec(mp,pnames)
      % 
      %  or
      %
      %   [pvec0,mp]=mp2pvec(mp,pnames)
      % 
      % OUTPUTS
      %   pvec:      k dimentional vector of scalars, the values of parameters in the model structure
      %                mp with field names that match elements of pnames. k equals the sum of all parameters
      %                (including vectorized versions of parameters that are cell arrays in mp) over all the
      %                different structure field names in pnames
      %
      %   mp:         updated with field mp.pvec_info (defined below)
      %
      % mp.pvec_info: a struct with 9 fields pnames_type,theta_idx,theta_names,alpha_idx,alpha_names,beta_idx,beta_names,price_idx,price_names
      %               pnames_type is a cell array with the same dimension as pnames that specifies the type of each parameter name, i.e. 'theta', 'alpha', etc
      %               The *_idx fields are vectors of indices, (1 to k) corresponding to elements of pvec whose parameter values correspond
      %               to a subset of field names of mp that are designated to be "* parameters". For example if * is 'theta' these are
      %               parameters that enter the current utility function in the auto equilibrium application, and the corresponding vector
      %               theta_names is a cell array with same dimension as theta_idx with the names of the parameters for
      %               each element of theta_idx. Similarly alpha_idx and alpha_names are the indices and names of parameters that determine
      %               the accident/aging probability matrix, beta_idx and beta_names is the index and name for the discount factor,
      %               and price_idx and price_names are the indices and names for secondary market prices
      %
      % INPUTS
      %   mp:   Structure with fields taking scalar values (at least one field for each value element in pnames) 
      %
      %   pnames:    k dimensional cell array, that specify names of the fields in mp that should be written to pvec
      %
      %   pnames_map:    result of 'cellfun(@(x) numel(mp.(x)),pnames)'
      %
      % See also: pvec2mp

      if nargin < 2;
          pnames=mp.pnames;
      end

      fieldnotfound=0;

      for i=1:numel(pnames)
        if (~isfield(mp,char(pnames(i))))
           fprintf('Warning mp2pvec: field %s not found in mp structure, deleting it from pnames\n',char(pnames(i)));
           pnames(i)=[];
           fieldnotfound=1;
        end
      end

      if (fieldnotfound)
       mp.pnames=pnames;
      end

      if (numel(mp.pnames) == 0)
         fprintf('Error mp2pvec: no relevant elements in the pnames cell, nothing returned\n');
         pvec0=[];
         return;
      end;
      
      mp.p0=getfields(mp.p0, mp.pnames); % delete reduntant fields in mp.p0 
      

      unique_pnames=unique(pnames);
      if (numel(unique_pnames) ~= numel(pnames))
         fprintf('Warning: mp2pvec: one or more duplicated name in the pnames argument. Duplicates remove and revised pnames is sorted alphabetically.\n');
          pnames=unique_pnames;
          mp.pnames=pnames;
      end

      spnames=numel(pnames);
      pnames_map=zeros(spnames,1);
      for i=1:spnames
        if any(strcmp(pnames(i), estim.alpha))
            if (size(mp.p0.(char(pnames(i))),1) > 1)
               fprintf('warning: estim.mp2pvec:  pnames %s has more than 1 row. Only 1 row allowed, so using only first row for this field in mp.p0.\n',char(pnames(i)));
               fprintf('(Current version of code only allows accident probabities to vary by car type but not yet by consumer type)\n' );
               mp.p0.(char(pnames(i)))=mp.p0.(char(pnames(i)))(1,:);
            end
        end
        if (isnumeric(mp.p0.(char(pnames(i)))))
          pnames_map(i)=1;
        else
          pnames_map(i)=numel(vertcat(mp.p0.(char(pnames(i))){:}));
        end
      end

      function [idx, names]=update_idx(idx, names, j, np_i, pname)
        % private function that updates idx and names for j'th
          idx=[idx; (j+1:j+np_i)'];
          for tn=j+1:j+np_i
            ind=numel(names)+1;
            names{ind}=char(pname);
          end
      end

      pvec0=zeros(sum(pnames_map),1);  %preallocate memory, column vector
      alpha_idx=[];beta_idx=[];theta_idx=[]; price_idx=[];
      alpha_names={};theta_names={};beta_names={};price_names={};
      pnames_type={};
      j=0;
      for i=1:numel(pnames_map)
        if iscell(mp.p0.(char(pnames(i))))
          pvec0(j+1:j+pnames_map(i))=[mp.p0.(char(pnames(i))){:}]; %by column fit into the vector
        else
          pvec0(j+1:j+pnames_map(i))=mp.p0.(char(pnames(i)))(:); %by column fit into the vector
        end
        if any(strcmp(pnames(i), estim.theta))
          pnames_type{i}='theta';
          [theta_idx, theta_names]=update_idx(theta_idx, theta_names, j, pnames_map(i), pnames(i));
        elseif any(strcmp(pnames(i), estim.alpha))
          pnames_type{i}='alpha';
          [alpha_idx, alpha_names]=update_idx(alpha_idx, alpha_names, j, pnames_map(i), pnames(i));
        elseif (strcmp(pnames(i),'bet'))
          pnames_type{i}='beta';
          [beta_idx, beta_names]=update_idx(beta_idx, beta_names, j, pnames_map(i), pnames(i));
        end
        j=j+pnames_map(i);
      end

      mp.pvec_info=cell2struct({j, theta_idx,theta_names,alpha_idx,alpha_names,beta_idx,beta_names,price_idx,price_names,pnames_type}, ...
                {'np', 'theta_idx','theta_names','alpha_idx','alpha_names','beta_idx','beta_names','price_idx','price_names','pnames_type'}, 2);

      % final check to be sure that the union of the returned names is the same as in pnames
      return_names=unique([alpha_names'; beta_names'; theta_names'; price_names']); 
      if (sum(strcmp(unique_pnames',return_names)) ~= numel(pnames))
         fprintf('Error in return from mp2pvec: the union of returned parameter names does not equal the names in the pnames argument\n');
         fprintf('Sorted, uniquen names in pnames\n');
         unique_pnames
         fprintf('Sorted, unique names in alpha_names, beta_names, theta_names and price_names\n');
         return_names
      end        
    end % end of mp2pvec

    function mp_out = pvec2mp(pvec,mp,update_type);
      %VEC2STRUCT: Procedure to convert a vector to a structure with scalar fields. 
      %
      %SYNTAX:      
      %   mp = pvec2mp(pvec,pvec_info)
      %
      %OUTPUTS:
      %   mp:   Structure with fields taking scalar values (at least one for each value element in pnames) 
      %
      %INPUTS
      %   pvec:      k dimentional vector of scalars
      %
      %   pvec_info: a structure containing a map of the parameter values in pvec with the names of fields 
      %              (see struct2vec for details on its fields)
      %
      %   Note that when there are multiple names and thus multiple elements of pvec corresponding to a given parameter name,
      %   the convention is to use vertical concatenation (vertcat function) to order the individual elements of a cell array,
      %   as this is the same convention used to encode the pvec by the struct2vec function
      %
      % See also: struct2vec

      param_types={'theta_names','alpha_names','beta_names','price_names'};
      index_types={'theta_idx','alpha_idx','beta_idx','price_idx'};

      for pt=1:3
       unique_names=unique(mp.pvec_info.(param_types{pt}));
       if (strcmp(index_types{pt},sprintf('%s%s', update_type, '_idx')) | strcmp(update_type,'all'))
        if (numel(unique_names))
          % now loop through the elements and update the mp with the values in pvec
          stn=numel(unique_names);
          for i=1:stn
            if (iscell(mp.p0.(char(unique_names(i)))))  % variable corresponds to cell array
              [nr,nc]=size(mp.p0.(char(unique_names(i))));
              mp_out.(char(unique_names(i)))=num2cell(reshape(pvec(mp.pvec_info.(char(index_types{pt}))(find(strcmp(mp.pvec_info.(char(param_types{pt})),char(unique_names(i)))))),nr,nc));
            else                                                    % variable corresponds to  a scalar
              mp_out.(char(unique_names(i)))=pvec(mp.pvec_info.(char(index_types{pt}))(find(strcmp(mp.pvec_info.(char(param_types{pt})),char(unique_names(i))))));
            end
          end
        end
       end
      end

      % if prices are being updated, slightly different code is required because price_j is a cell array of vectors of possibly different lengths
      unique_names=unique(mp.pvec_info.(param_types{4}));
      if (numel(unique_names))
        [nr,nc]=size(mp.p0.(char(unique_names)));
        prices=pvec(mp.pvec_info.price_idx);
        sr=0;
        er=0;
        for i=1:nc
          nprices=numel(mp.p0.(char(unique_names)){i});
          sr=er+1;
          er=sr+nprices-1;
          mp_out.price_j{i}=prices(sr:er);
        end
      end      
    end % end of pvec2mp

    function [mp] = update_mp(pvec0, mp, update_type)
      % Updates mp parameter structure with values in pvec0 
      % Syntax: mp=estim.update_mp(pvec0, mp)
      % Inputs: 
      %   pvec0 and pvec_info:  generated by estim.mp2pvec
      %   mp:                   full mp struct to be updated
      %   update_type:          specifies which components of pvec0 to be updated:
      %                         ('theta_idx','alpha_idx','beta_idx','price_idx')
      %                         or if omitted, all components of pvec0 are used to update mp

      if (nargin < 3)
          update_type='all';
      end

      % get structure with updated parameters to be estimated (only those in pvec0)
      mp0=estim.pvec2mp(pvec0,mp,update_type);

      % add mp0 to complete parameter struct, mp
      mp=combinestructs(mp, mp0);

      % and update (this will fill out parameters for all consumer and car types if not populated)
      mp=trmodel.update_mp(mp);
    end

    function [g] = matgrad_mp(fmp, mp, pvec0, pnames_type)
      % takes derivatives of a matrix function that depends on model parameters mp
      % pvec0 is the vector where to evaluate derivatives. 
      % pnames_type specifies the subset of pvec0  components that derivatives are taken with respect to
      %             ('theta','alpha','beta','price') or if not provided, then all
      %             components of pvec0 are used

      % updat mp0
      mp=estim.update_mp(pvec0, mp, pnames_type);

      if ~strcmp(pnames_type, 'all')
        mp.pnames=mp.pnames(strcmp(mp.pvec_info.pnames_type,pnames_type));
      end 

      [pvec0,mp] = estim.mp2pvec(mp,mp.pnames);

      function [fvec] = mp_wrapper(fmp, pvec0, mp)
        mp=estim.update_mp(pvec0, mp);
        f=fmp(mp);
        fvec=reshape(f, numel(f),1); 
      end


      f=fmp(mp);

      szf=size(f);
      if sum(szf(2:end))==1;
        szf=numel(f);
      end 


      fvec=@(x) mp_wrapper(fmp, x, mp);

      g=gradp(fvec,pvec0);
      
      g=reshape(g, [szf, numel(pvec0)]);
    end

    function [f, g, h] = nloglike(mp, s, dta, pvec)
      
      logl = @(pvec) estim.loglike(mp, s, dta, pvec);
      
      N=dta.count;

      if nargout==1
          f=-sum(logl(pvec))/sum(N);
      elseif nargout>1
          if estim.analytical
              [logl_i, score]=logl(pvec);
          else
              [logl_i]=logl(pvec);
              score=gradp(logl, pvec);
          end
          f=-sum(logl_i)/sum(N);
          g=-sum(score)/sum(N);
      end
      if nargout>2
          score = bsxfun(@rdivide, score, sqrt(max(N,1)));
          % s=s./sqrt(max(N,1));
          h=score'*score/sum(N);
      end
    end

    function [df_err, df_nm, df_an] =check_deriv(fun, x0)
      [f_nm, df_an]=fun(x0); 
      df_nm=gradp(fun,x0); 
      df_err=df_an-df_nm;
      if any(norm(df_err)>estim.deriv_tol)
        warning('Analytical derivatives do not match at initial value, x0')
      else
        fprintf('Analytical and numerical derivatives match at initial value, x0\n')
      end 
    end % end of estim.check_deriv

    function [logl, score] = loglike(mp, s, dta, pvec)
      % loglike(): full loglikelihood, including discrete car choice and scrappage decision (where appropriate). 
      
      dta=dta((~isnan(dta.count)) & ~isnan(dta.count_scrap),:);
      assert(all(dta.is > 0), 'is must be base 1');
      assert(not(any(isnan(dta.is) | isnan(dta.id))), 'code cannot handle NaN values')
      
      % 1.update mp struct witb values in pvec
      [mp] = estim.update_mp(pvec, mp);

      save('results/estimation/mle_iter.mat','mp');
      
      % 2 solve model
      global p0
      [sol]=equilibrium.solve(mp, s, p0);
      p0=sol.p;
      
      np=0;
      if nargout==2
          [grad_dp_dmp, grad_dp, grad_dmp]=g_trmodel.total_deriv(mp, s, sol);
          np=size(grad_dp_dmp.ccp_tau{1}, 3);
      end
      
      % 3 evaluate likelihood
      ncell=size(dta,1);
      logl=nan(ncell, 1);
      score=nan(ncell, np);

      
      dta.nocar  = (dta.is == s.is.nocar)                           & (~isnan(dta.count)) & ~isnan(dta.count_scrap);
      dta.trade  = (dta.is ~= s.is.nocar) & (dta.id ~= s.id.keep)   & (~isnan(dta.count)) & ~isnan(dta.count_scrap);
      dta.keep   = (dta.is ~= s.is.nocar) & (dta.id == s.id.keep)   & (~isnan(dta.count)) & ~isnan(dta.count_scrap);
      dta.scrap  = (dta.d_scrap == 1)                               & (~isnan(dta.count)) & ~isnan(dta.count_scrap);


      for tau=1:mp.ntypes
          i_tau=dta.tau==tau;
          dt=dta(i_tau,:);
          idx = sub2ind([s.ns, s.nd], dt.is, dt.id); 
          is  = dt.is; 

          Ntau = size(dt, 1); 
          Pi = zeros(Ntau, 1); 

          % 4 Cases: 
          % case 1:  No car (including staying with no car)
          % case 2a: Trading car - scrapping is observed
          % case 2b: Trading car - scrapping NOT observed
          % case 3:  Keep existing car

          % ---- case 1 (no car)  ---  
          % People who are in the no car state 
          % Likelihood is just the choice probability 
          % since they have no car that could be in an accident or one they could scrap or sell. 
          I1=dt.nocar==1;
          Pi(I1)= sol.ccp_tau{tau}(idx(I1));

          % ---- case 2 (trading) ---
          % People who had a car and sold the car (purge or trade)
          % If car was scrapped, we do not know if this was a "forced scrappage" due to an accident 
          % (since accidents are unobserved) or a voluntary scrappage (since sell/scrap decisions are unobserved).   
                    
          % ---- case 2a: Trading car - scrapping is observed ---  
          %      The likelihood when we see the previous car was scrapped is given by equation (192) and involves 
          %      -  Probability of trading a car invent of accident (existing car becomes clunker)
          %      -  Probability of trading a car invent of no accident (exiting car ages)  
          %      -  the endogenous scrap probability for the existing car 
          %      -  the accident probability for the existing car (before aging)
          I2a=dt.trade==1 & dt.scrap==1; 

          is_clunker = [s.is.clunker{dt(I2a, :).s_car_type}]';
          idx_clunker = sub2ind([s.ns, s.nd], is_clunker, dt.id(I2a)); 

          Pcar_trade_scrap      = sol.ccp_tau{tau}(idx(I2a));
          Pclunker_trade_scrap  = sol.ccp_tau{tau}(idx_clunker);
          Ps_trade_scrap        = sol.ccp_scrap_tau{tau}(is(I2a)); 
          Pa_trade_scrap        = sol.F.accprob(is(I2a));  

          Pi(I2a) =  Pcar_trade_scrap.*Ps_trade_scrap.*(1-Pa_trade_scrap) ...
                   + Pclunker_trade_scrap.*Pa_trade_scrap ;

          % ---- case 2b: Trading car - scrapping NOT observed ---  
          %      The probability for that case is given by equation (193).  
          %      Product of the choice probability, the probability of no voluntary scrappage and the probability no accident occurred.
          I2b=dt.trade==1 & dt.scrap~=1;
          Pcar_trade_noscrap = sol.ccp_tau{tau}(idx(I2b));
          Ps_trade_noscrap=sol.ccp_scrap_tau{tau}(is(I2b)); 
          Pa_trade_noscrap=sol.F.accprob(is(I2b));  

          Pi(I2b) =  Pcar_trade_noscrap.*(1-Ps_trade_noscrap).*(1-Pa_trade_noscrap);

          % ---- case 3.  (keeping) ---
          % Finally we also have have the likelihood for keeping the previous car. 
          % This immediately implies there is no voluntary scrappage and there was no accident, 
          % and the likelihood for this is given in equation (194).  
          I3=dt.keep==1; 
          Pcar_keep= sol.ccp_tau{tau}(idx(I3));
          Pa_keep=sol.F.accprob(is(I3));  

          Pi(I3) =  Pcar_keep.*(1-Pa_keep);

          % ---- Collect results ---- 
          Pi=max(Pi, 1e-10);
          logl(i_tau,:)=log(Pi).*dt.count;

          if nargout==2
              dPi=nan(Ntau,np);
              dccpi=nan(numel(idx),np);
              dccpi_clunker=nan(numel(idx_clunker),np);
              for ip=1:np
                  dccp=grad_dp_dmp.ccp_tau{tau}(:,:,ip);
                  dccpi(:,ip)=dccp(idx);
                  dccpi_clunker(:,ip)=dccp(idx_clunker);
              end

              % The 3 Cases
              % ---- case 1 (no car)  ---  
              dPi(I1,:)= dccpi(I1,:);

              % ---- case 2a: Trading car - scrapping is observed ---  
              dPcar_trade_scrap       = dccpi(I2a,:);
              dPclunker_trade_scrap   = dccpi_clunker;
              dPs_trade_scrap         = grad_dp_dmp.ccp_scrap_tau{tau}(is(I2a),:); 
              dPa_trade_scrap         = grad_dmp.accprob(is(I2a),:);  

              dPi(I2a,:) =  ... 
              + dPcar_trade_scrap.*Ps_trade_scrap.*(1-Pa_trade_scrap) ... 
              + Pcar_trade_scrap.*dPs_trade_scrap.*(1-Pa_trade_scrap) ... 
              + Pcar_trade_scrap.*Ps_trade_scrap.*(-dPa_trade_scrap) ...
              + dPclunker_trade_scrap.*Pa_trade_scrap ...
              + Pclunker_trade_scrap.*dPa_trade_scrap ;

              % ---- case 2b: Trading car - scrapping NOT observed ---  
              dPcar_trade_noscrap = dccpi(I2b,:);
              dPs_trade_noscrap=grad_dp_dmp.ccp_scrap_tau{tau}(is(I2b),:);
              dPa_trade_noscrap=grad_dmp.accprob(is(I2b),:); 

              dPi(I2b,:) =  ... 
                dPcar_trade_noscrap.*(1-Ps_trade_noscrap).*(1-Pa_trade_noscrap) ...
                + Pcar_trade_noscrap.*(-dPs_trade_noscrap).*(1-Pa_trade_noscrap) ...
                + Pcar_trade_noscrap.*(1-Ps_trade_noscrap).*(-dPa_trade_noscrap);

              % ---- case 3.  (keeping) ---
              dPcar_keep = dccpi(I3,:);
              dPa_keep=grad_dmp.accprob(is(I3),:); 
              dPi(I3,:)    = dPcar_keep.*(1-Pa_keep) + Pcar_keep.*(-dPa_keep);  

              score(i_tau,:)=(dPi./Pi).*dt.count;

          end
      end
    end % end loglike()

    function [logl, score] = loglike_old(mp, s, dta, pvec)
      % loglike(): full loglikelihood, including discrete car choice 
      % and scrappage decision (where appropriate). 
      
      assert(all(dta.is > 0), 'is must be base 1');
      assert(not(any(isnan(dta.ipt) | isnan(dta.is) | isnan(dta.id))), 'code cannot handle NaN values')
      
      % 1.update mp struct witb values in pvec
      [mp] = estim.update_mp(pvec, mp);

      save('results/estimation/mle_iter.mat','mp');
      
      % 2 solve model
      global p0
      [sol]=equilibrium.solve(mp, s, p0);
      p0=sol.p;
      
      np=0;
      if nargout==2
          [grad_dp_dmp, grad_dp, grad_dmp]=g_trmodel.total_deriv(mp, s, sol);
          np=size(grad_dp_dmp.ccp_tau{1}, 3);
      end
      
      % 3 evaluate likelihood
      ncell=size(dta,1);
      logl=nan(ncell, 1);
      score=nan(ncell, np);
      
      dta.trade = (dta.is ~= s.is.nocar) & (dta.id ~= s.id.keep) & (~isnan(dta.count)) & ~isnan(dta.count_scrap);
      dta.keep  = (dta.id == s.id.keep)  & (~isnan(dta.count)) & ~isnan(dta.count_scrap);
      for tau=1:mp.ntypes
          i_tau=dta.tau==tau;
          dta_tau=dta(i_tau,:);
          Ntau = size(dta_tau, 1); 
          
          % 3.a --- no accident --- 
          % likelihood contribution from car choice
          ii_noacc = sub2ind([s.ns, s.nd], dta_tau.is, dta_tau.id); 
          Pi_car_noacc=sol.ccp_tau{tau}(ii_noacc);

          % likelihood contribution from scrapping 
          is_trade=dta_tau(dta_tau.trade,:).is; 
          dscrap_i=dta_tau.d_scrap(dta_tau.trade); 
          pscrap_i = sol.ccp_scrap_tau{tau}(is_trade); 
          Pi_scrap_noacc = ones(Ntau, 1); % for non-trading households, the scrappage contribution does not enter (so = 1)
          Pi_scrap_noacc(dta_tau.trade)= (1-pscrap_i) - dscrap_i + 2*dscrap_i.*pscrap_i; 

          % joint probability of the car and scrappage decisions conditional on no accident having occurred 
          Pi_noacc = Pi_car_noacc .* Pi_scrap_noacc; 

          % 3.b --- accident --- 
          % in the model, Pr(keep|is = clunker) = 0, so in this case, the household is forced to trade. 
          is_clunker = [s.is.clunker{dta_tau(dta_tau.trade, :).s_car_type}]';
          ii_acc = sub2ind([s.ns, s.nd], is_clunker, dta_tau.id(dta_tau.trade)); 
          
          Pi_car_acc = zeros(Ntau, 1); % non-traders cannot have had an accident: zero probability event
          Pi_car_acc(dta_tau.trade) = sol.ccp_tau{tau}(ii_acc); 
          
          % joint probability of the car and scrappage decisions conditional on accident 
          % scrappage choice: deterministic in this case, Pr(scrap|acc) = 100%. 
          Pi_acc = zeros(Ntau, 1); 
          Pi_acc(dta_tau.trade) = Pi_car_acc(dta_tau.trade) .* dscrap_i; 

          % 3.c --- mixed likelihood ---

          % pr_acc=sol.F.accprob(dta_tau.is);

          % keep:  pr_acc(dta_tau.keep)=sol.F.accprob(dta_tau.ipt(keep));
          % no car pr_acc(dta_tau.keep)=1;
                  
          % pr_acc=zeros(Ntau,1);            
          pr_acc = sol.F.accprob(dta_tau.ipt);
          % pr_acc=sol.F.accprob(dta_tau.is);

          Pi = Pi_noacc .* (1-pr_acc) + Pi_acc .* pr_acc; 
          Pi=max(Pi, 1e-10);

          logl(i_tau,:)=log(Pi).*dta_tau.count;

          if nargout==2
              dPi_car_noacc=zeros(Ntau,np);
              dPi_car_acc=zeros(Ntau,np);
              dPi_scrap_noacc=zeros(Ntau,np);
              dPi_acc = zeros(Ntau, np);
              
              for ip=1:np
                  dccp=grad_dp_dmp.ccp_tau{tau}(:,:,ip);
                  dPi_car_noacc(:,ip)=dccp(ii_noacc);
                  dPi_car_acc(dta_tau.trade,ip)=dccp(ii_acc);
              end

              dpscrap_i=grad_dp_dmp.ccp_scrap_tau{tau}(is_trade,:);
              % Derivatve of Pi_scrap_noacc(dta_tau.trade)  = (1-pscrap_i) - dscrap_i + 2*dscrap_i.*pscrap_i; 
              dPi_scrap_noacc(dta_tau.trade,:)= -dpscrap_i  + 2*dscrap_i.*dpscrap_i; 
              % Derivatve of Pi_noacc = Pi_car_noacc .* Pi_scrap_noacc; 
              dPi_noacc = dPi_car_noacc .* Pi_scrap_noacc + Pi_car_noacc .* dPi_scrap_noacc;  

              %Derivatve of Pi_acc(dta_tau.trade)    = Pi_car_acc(dta_tau.trade)   .* dscrap_i; 
              dPi_acc(dta_tau.trade,:) = dPi_car_acc(dta_tau.trade,:) .* dscrap_i; 

              % dpr_acc=zeros(Ntau,np); 
              dpr_acc=grad_dmp.accprob(dta_tau.ipt,:);
              % dpr_acc = grad_dmp.accprob(dta_tau.is,:);
              
              %Derivatve of Pi =  Pi_noacc .* (1-pr_acc) + Pi_acc .* pr_acc; 
              dPi = dPi_noacc .* (1-pr_acc) - Pi_noacc .* dpr_acc ...
                  + dPi_acc .* pr_acc + Pi_acc .* dpr_acc; 
              
              score(i_tau,:)=(dPi./Pi).*dta_tau.count;
          end
      end
    end % end loglike()
    
    function [out]=outcomes(mp, s, dta);
      % estim.outcomes: Estimates of model outcomes (based on tabulated choice data)
      %
      % SYNTAX: [out] = estim.outcomes(mp, s, dta);
      %
      % INPUT: 
      %   mp:   struct with model parameters (must have field "ntypes")
      %   s:    struct model indexes for states and deciusions (need to be consistent with dta and mp) 
      %   dta:  Table with choice data dta must have the following variables: "count", "tau", "is" and "id"'
      % 
      % OUTPUT: 
      %   out: structure with fields 
      %     'ccp_tau', 'ccp_scrap_tau', 'q_tau', 'q', 'h_tau', 'marketshares'
      [ccp_tau] = estim.ccp_tau(mp, s, dta);
      [ccp_scrap_tau] = estim.ccp_scrap_tau(mp, s, dta);
      [h_tau] = estim.h_tau(mp, s, dta);
      [q_tau, q] = estim.q_tau(mp, s, dta);
      [marketshares] = equilibrium.marketshares(mp, s, h_tau);
      out=cell2struct({ccp_tau, ccp_scrap_tau, q_tau, q, h_tau, marketshares},{'ccp_tau', 'ccp_scrap_tau', 'q_tau', 'q', 'h_tau', 'marketshares'}, 2);
    end

    function [ccp_scrap_tau_d] = ccp_scrap_tau(mp, s, dta) 
      % estim.ccp_scrap_tau: Frequency estimator of scrap probabilities (based on tabulated choice data)
      %                      ccp_scrap_tau is computed condtional on trading (i.e. not keeping) 
      %
      % SYNTAX: [ccp_scrap_tau_d] = estim.ccp_scrap_tau(mp, s, dta)
      %
      % INPUT: 
      %   mp:   struct with model parameters (must have field "ntypes")
      %   s:    struct model indexes for states and deciusions (need to be consistent with dta and mp) 
      %   dta:  Table with choice data dta must have the following variables: "count", "tau", "is" and "id"'
      % 
      % OUTPUT: 
      %   ccp_scrap_tau_d: (ntypes*T cell) Frequency estimator of scrap CCPs (conditional on trading (i.e. not keeping))
      ccp_scrap_tau = cell(mp.ntypes,1);  
        for tau=1:mp.ntypes  
          % tab=data.mean(dta, @(data) data.d_scrap, {'is'}, (dta.tau==tau & dta.id~=s.id.keep & dta.is~=s.is.nocar));
          tab=data.mean(dta, @(data) data.d_scrap, {'is'}, (dta.tau==tau & dta.id~=s.id.keep & dta.is~=s.is.nocar));
          ccp_scrap_tau_d{tau}=nan(s.ns, 1);
          ccp_scrap_tau_d{tau}(tab.is)=tab.mean;
        end
    end
        
    function [ccp_tau_d] = ccp_tau(mp, s, dta)
      % estim.ccp_tau: Frequyency estimator of CCPs based on tabulated choice data
      %
      % SYNTAX: [ccp_tau_d] = estim.ccp_tau(mp, s, dta)
      %
      % INPUT: 
      %   mp:   struct with model parameters (must have field "ntypes")
      %   s:    struct model indexes for states and deciusions (need to be consistent with dta and mp) 
      %   dta:  Table with choice data dta must have the following variables: "count", "tau", "is" and "id"'
      % 
      % OUTPUT: 
      %   ccp_tau_d: (ntypes*T cell) Frequyency estimator of CCPs
                 
      assert(all(ismember({'count', 'is','id', 'tau'}, dta.Properties.VariableNames)), 'Table must have the following variables: "count", "tau", "is" and "id"');

      ccp_tau_d=cell(mp.ntypes,1);
      for tau=1:mp.ntypes
        tab=data.tabulate(dta, {'is', 'id'}, {}, dta.tau == tau);
        i=sub2ind([s.ns, s.nd], tab.is, tab.id); 
        ccp_tau_d{tau}=nan(s.ns, s.nd);
        ccp_tau_d{tau}(i)=tab.count;
        ccp_tau_d{tau}=ccp_tau_d{tau}./nansum(ccp_tau_d{tau},2);
      end    
    end

    function [h_tau_d] = h_tau(mp, s, dta)
      % estim.h_tau: Frequency estimator of post trade holdings distribution
      %
      % SYNTAX: [h_tau_d] = estim.h_tau(mp, s, dta)
      %
      % INPUT: 
      %   mp:   struct with model parameters (must have field "ntypes")
      %   s:    struct model indexes for states and deciusions (need to be consistent with dta and mp) 
      %   dta:  Table with choice data dta must have the following variables: "count", "tau", "is" and "id"'
      % 
      % OUTPUT: 
      %   h_tau_d: (ntypes*T cell) Frequyency estimator of h_hau
                 
      assert(all(ismember({'count', 'ipt', 'tau'}, dta.Properties.VariableNames)), 'Table must have the following variables: "count", "tau", "ipt"');
      h_tau_d=cell(mp.ntypes,1);
      for tau=1:mp.ntypes
        tab=data.tabulate(dta(dta.tau == tau,:), {'ipt'});
        h_tau_d{tau}=nan(s.ns, 1);
        h_tau_d{tau}(tab.ipt)=tab.count;
        h_tau_d{tau}=h_tau_d{tau}./nansum(h_tau_d{tau},1)*mp.tw(tau);
      end    
      
    end

    function [q_tau_d, q_d] = q_tau(mp, s, dta)
      % estim.q_tau: Frequyency estimator of car age distribution (state)
      %
      % SYNTAX: [q_tau_d, q_d] = estim.q_tau(mp, s, dta)
      %
      % INPUT: 
      %   mp:   struct with model parameters (must have field "ntypes")
      %   s:    struct model indexes for states and deciusions (need to be consistent with dta and mp) 
      %   dta:  Table with choice data dta must have the following variables: "count", "tau", "is" and "id"'
      % 
      % OUTPUT: 
      %   q_tau_d: (ntypes*T cell) Frequyency estimator of q_tau
                 
      assert(all(ismember({'count', 'is','id', 'tau'}, dta.Properties.VariableNames)), 'Table must have the following variables: "count", "tau", "is" and "id"');
      q_tau_d=cell(mp.ntypes,1);
      for tau=1:mp.ntypes
        tab=data.tabulate(dta(dta.tau == tau,:), {'is'});
        i=sub2ind([s.ns, 1], tab.is); 
        q_tau_d{tau}=nan(s.ns, 1);
        q_tau_d{tau}(i)=tab.share;
      end    
      q_d=[q_tau_d{:}]*mp.tw(:);
    end

    function moments_w = moment_weights(moments_id, momlist, w_list, dta)
        % INPUTS: 
        %   moments_id: M-vector of moment values 
        %   momlist: num_mom cell array of names of moment blocks 
        %   w_list: w_list of weights (summing to 1.0) to be assigned to
        %   each moment block 
        
        % 0. dimensions
        num_mom = numel(momlist); 
        M = numel(moments_id); 
        if nargin < 3
            w_list = 1/num_mom * ones(size(momlist)); 
        end
        
        % 1. checks 
        assert(all(moments_id>=1) && all(moments_id <= num_mom), 'Unexpected values in moments_id'); 
        assert(numel(w_list) == numel(momlist), 'w_list should have same dim. as momlist.');
        assert(abs(sum(w_list)-1.0)< 1e-12, 'w_list must sum to 1');
        
        % 2. assign weights 
        moments_w = nan(M,1); 
        
        for i=1:num_mom
            % 2.a identify block belonging to momlist{i}
            I = moments_id == i; 
            
            % 2.b assign weights 
            moments_w(I) = w_list(i)/sum(I); 
            
            % 2.c exception for likelihood moment 
            switch momlist{i}
                case 'likelihood'
                    adj = dta.count ./ mean(dta.count); % mean(adj) = 1.0 
                    moments_w(I) = moments_w(I) .* adj; 
            end
        end
        
        s = sum(moments_w); 
        assert(abs(s-1.0) < 1e-8, 'Weights do not sum to 1: unexpected!'); 
    end
        
    function [moments, labels, moments_id] = moments_dta(mp, s, dta, momlist)
      % OUTPUT: 
      %   moments: M-vector of moment-values (M cannot be determined ex
      %   ante)
      %   labels: M-cell of labels 
      %   moments_id: M-vector of integers, identifying the "block" of
      %   moments each entry belongs to (for the purpose of creating
      %   weights 

      numInit = 10000; 
      moments=nan(numInit,1);
      moments_id=nan(numInit,1);
      labels=cell(numInit,1);
      im=1;
      
      % compute tau shares
      for i=1:numel(momlist)
          switch momlist{i}
              case 'scrap' 
                  pr_scrap = tabdat.get_scrap_data(mp, dta); 
                  [Ncartot, Nhh] = size(pr_scrap); 
                  lab_hh  = repmat(mp.lbl_types', Ncartot, 1); 
                  
                  lab_car = cell(Ncartot, 1); 
                  for j_ = 1:mp.ncartypes
                      jj_ = s.is.car{j_}; 
                      lab_car(jj_) = sprintfc([mp.lbl_cartypes{j_} ', age %d'], s.is.age(jj_)); 
                  end
                  lab_car = repmat(lab_car, 1, Nhh); 
                  
                  Idrop = [s.is.clunker{:}]'; 
                  pr_scrap(Idrop, :) = []; 
                  lab_car(Idrop, :) = []; 
                  
                  num_new = numel(pr_scrap); % not clunker  
                  ii = im : im+num_new-1; 
                  moments(ii) = pr_scrap(:); 
                  moments_id(ii) = i; 
                  for i_=1:numel(ii)
                      labels{ii(i_)} = sprintf('Scrap: %s, %s', lab_car{i_}, lab_hh{i_}); 
                  end
                  
                  % update 
                  im = im + num_new; 
              case 'q_tau'
                  % uses average over the years provided in the input data 
                  [q_tau, q_tau_avg] = tabdat.create_q_tau(dta); 
                  for tau=1:mp.ntypes
                      v = q_tau_avg{tau}; 
                      %I = not(isnan(s.ipt.age)); % isnan => nocar 
                      I = (s.is.age ~= 0) & not(isnan(s.ipt.age)); % isnan => nocar 
                      v = v(I); 
                      
                      % store 
                      ii = im : im+numel(v)-1; 
                      moments(ii) = v; 
                      labels(ii) = sprintfc(['q: ' mp.lbl_types{tau} ' car %d'], 1:numel(ii)); 
                      moments_id(ii) = i; 
                      
                      % iterate 
                      im = im+numel(ii); 
                  end
              case 'h_tau' 
                  % uses average over time 
                  [h_tau, h_tau_avg] = tabdat.create_h_tau(dta); 
                  for tau=1:mp.ntypes
                      v = h_tau_avg{tau}; 
                      I = (s.ipt.age ~= 0) & not(isnan(s.ipt.age)); % isnan => nocar 
                      v = v(I); 
                      
                      % store 
                      ii = im : im+numel(v)-1; 
                      moments(ii) = v; 
                      labels(ii) = sprintfc([mp.lbl_types{tau} ' car %d'], 1:numel(ii)); 
                      moments_id(ii) = i; 
                      
                      % iterate 
                      im = im+numel(ii); 
                  end
              case 'q_nocar'
                  [q_tau, q_tau_avg, q_tau_nocar] = tabdat.create_q_tau(dta); 
                  
                  % store 
                  ii = im:im+mp.ntypes-1; 
                  moments(ii) = mean(q_tau_nocar, 2); % avg. over time if any 
                  moments_id(ii) = i; 
                  for itau=1:mp.ntypes
                      i_ = ii(itau); 
                      labels(i_) = {sprintf('%s no car', mp.lbl_types{itau})}; 
                  end
                  
                  im = im+numel(ii); 
              
              case 'h_nocar'
                  [h_tau, h_tau_avg, h_tau_nocar] = tabdat.create_h_tau(dta); 
                  
                  % store 
                  ii = im:im+mp.ntypes-1; 
                  moments(ii) = mean(h_tau_nocar, 2); % avg. over time if any 
                  moments_id(ii) = i; 
                  for itau=1:mp.ntypes
                      i_ = ii(itau); 
                      labels(i_) = {sprintf('%s no car', mp.lbl_types{itau})}; 
                  end
                  
                  im = im+numel(ii); 
              
              case 'marketshares_tau_j'
                  % car market shares 
                  for j=1:mp.ncartypes
                      for tau=1:mp.ntypes
                          moments(im)=sum(dta.count(dta.s_car_type==j & dta.tau==tau))/sum(dta.count);
                          labels{im} = sprintf('market share (%s, %s)', mp.lbl_types{tau}, mp.lbl_cartypes{j});
                          moments_id(im) = i;
                          im=im+1;
                      end
                  end
                  
                  % no-car shares
                  for tau=1:mp.ntypes
                      moments(im)=sum(dta.count(dta.s_car_type==-1 & dta.tau==tau))/sum(dta.count);
                      labels{im} = sprintf('market share (%s, %s)', mp.lbl_types{tau}, 'no car');
                      moments_id(im) = i;
                      im=im+1;
                  end
              case 'likelihood'
                  % for each (tau, year, is), compute the fraction choosing
                  % each id 
                  
                  N = size(dta, 1); 
                  [prob, prob_tab] = estim.empirical_choice_probs(dta); 
                  
                  % write 
                  ii = im: im+N-1; 
                  moments(ii) = prob(:); 
                  labels(ii) = sprintfc('likelihood: (tau, year, is, id) = (%d, %d, %3d, %3d)', table2array(prob_tab(:, {'tau','year','is','id'}))); 
                  moments_id(ii) = i*ones(N,1); 
                  
                  % update 
                  im = im + numel(ii); 
              case 'pr_keep'
                  % for each (tau, year, is), compute the fraction choosing
                  % each id 
                  
                  %N = size(dta, 1); 
                  [prob, prob_tab] = estim.empirical_choice_probs(dta); 
                  I = prob_tab.id == s.id.keep; 
                  N = sum(I); 
                  prob = prob(I); 
                  prob_tab = prob_tab(I, :); 
                  
                  % write 
                  ii = im: im+N-1; 
                  moments(ii) = prob(:); 
                  labels(ii) = sprintfc('likelihood: (tau, year, is, id) = (%d, %d, %3d, %3d)', table2array(prob_tab(:, {'tau','year','is','id'}))); 
                  moments_id(ii) = i*ones(N,1); 
                  
                  % update 
                  im = im + numel(ii); 
              case 'mean_carage_tau_j'
                  tab=data.mean(dta, @(data) data.s_car_age ,{'tau', 's_car_type'}, dta.s_car_type>0);
                  for tau=1:mp.ntypes
                      for j=1:mp.ncartypes
                          moments(im)=tab.mean(tab.tau==tau & tab.s_car_type==j);
                          labels{im}= sprintf('mean car age (%s, %s)', mp.lbl_types{tau}, mp.lbl_cartypes{j});
                          moments_id(im) = i; 
                          im=im+1;
                      end
                  end
                  
              case 'old_cars'
                  holdings=data.tabulate(dta, {'d_car_age'}, {}, dta.s_car_type>0);
                  moments(im)=sum(holdings.share(holdings.d_car_age>=estim.oldcarage));
                  labels{im}= sprintf('holdings of cars >= %d y/o', estim.oldcarage);
                  moments_id(im) = i; 
                  im=im+1;
                  
              case 'old_cars_tau'
                  moments(im)=0;
                  for tau=1:mp.ntypes
                      moments(im)=0;
                      holdings=data.tabulate(dta, {'d_car_age'}, {}, dta.s_car_type>0 & dta.tau==tau);
                      moments(im)=sum(holdings.share(holdings.d_car_age>=estim.oldcarage));
                      labels{im}= sprintf('holdings of cars >= %d y/o (%s)', estim.oldcarage, mp.lbl_types{tau});
                      moments_id(im) = i; 
                      im=im+1;
                  end
                  
              case 'new_cars_tau'
                  moments(im)=0;
                  newcarage=2;
                  for tau=1:mp.ntypes
                      moments(im)=0;
                      holdings=data.tabulate(dta, {'d_car_age'}, {}, dta.s_car_type>0 & dta.tau==tau);
                      moments(im)=sum(holdings.share(holdings.d_car_age<=newcarage));
                      labels{im}= sprintf('holdings of cars 0-%d y/o (%s)', newcarage, mp.lbl_types{tau});
                      moments_id(im) = i; 
                      im=im+1;
                  end
              otherwise
                  disp('moment not implemented')
          end
      end
      
      % delete unused elements 
      ii = 1:im-1; 
      moments=moments(ii);
      labels = labels(ii); 
      moments_id = moments_id(ii); 
    end

    function [moments, dmoments] = moments_model(mp, s, momlist, dta)

      global p0;
      [sol]=equilibrium.solve(mp, s, p0);
      p0=sol.p;


      deriv=0; 
      dmoments=nan(100, 1);
      if nargout>1
        deriv=1;  
        [grad_dp_dmp, grad_dp, grad_dmp]=g_trmodel.total_deriv(mp, s, sol);
        dmarketshares=g_trmodel.dmarketshares_dp(mp, s, grad_dp_dmp.h_tau);
        np=size(grad_dp_dmp.h_tau{1},2); 
        dmoments=zeros(100, np);
      end
       
      
      moments=nan(100,1);
      im=1;
      for i=1:numel(momlist)
          switch momlist{i}
              case 'scrap' 
                  for tau=1:mp.ntypes
                      for j=1:mp.ncartypes
                          % pointers 
                          ii_read  = s.is.car_ex_clunker{j};
                          ii_write = im : im+numel(ii_read)-1;
                          
                          % write 
                          moments(ii_write) = ccp_scrap_tau{tau}(ii_read);
                          
                          if deriv
                              error 'derivs not implemented for scrap';
                          end
                          
                          % update 
                          im = im + numel(ii_write);
                      end
                  end
              case 'q_tau'
                  for tau=1:mp.ntypes
                      % where to read 
                      I = (s.is.age ~= 0) & not(isnan(s.is.age)); 
                      
                      % where to write
                      ii = im : im + sum(I)-1; 
                      
                      moments(ii, :) = q_tau{tau}(I); 
                      if deriv 
                          dmoments(ii, :) = grad_dp_dmp.q_tau{tau}(I, :); 
                      end
                      
                      % update 
                      im = im+numel(ii); 
                  end
              case 'q_nocar'
                  for tau=1:mp.ntypes
                      I = s.is.nocar;
                      moments(im,:) = q_tau{tau}(I); 
                      if deriv 
                          dmoments(im,:) = grad_dp_dmp.q_tau{tau}(I,:); 
                      end
                      im = im+1; 
                  end
              case 'marketshares_tau_j'
                  % car types... j=1,..,mp.ncartypes
                  for j=1:mp.ncartypes
                      for tau=1:mp.ntypes
                          moments(im)=marketshares(tau,j);
                          if deriv
                              dmoments(im,:)=dmarketshares(tau,j,:);
                          end
                          im=im+1;
                      end
                  end
                  % no car
                  for tau=1:mp.ntypes
                      moments(im)=marketshares(tau,end);
                      if deriv
                          dmoments(im,:)=dmarketshares(tau,end,:);
                      end
                      im=im+1;
                  end
                  
              case 'likelihood' 
                  % likelihood
                  N=size(dta,1);
                  like=nan(N, 1);
                  for tau=1:mp.ntypes
                      i_tau=dta.tau==tau;
                      dta_tau=dta(i_tau,:);
                      ii = sub2ind([s.ns, s.nd], dta_tau.is, dta_tau.id); 
                      %ii=(dta_tau.is-1)*s.nd+dta_tau.id; %through indexes of the chosen alternatives
                      Pi=ccp_tau{tau}(ii);
                      Pi=max(Pi, 1e-5);
                      like(i_tau) = Pi; 
                  end
                  
                  % write 
                  ii = im : im+N-1; 
                  moments(ii) = like(:); 
                  if deriv 
                      error 'not implemented'; 
                  end
                  
                  % update 
                  im = im + numel(ii); 
              case 'pr_keep'
                  % Pr(keep)
                  N = size(dta,1);
                  like = cell(mp.ntypes, 1); 
                  %like = nan(N, 1);
                  N_ = 0; 
                  for tau=1:mp.ntypes
                      i_tau = (dta.tau==tau) & (dta.id == s.id.keep);
                      N_ = N_ + sum(i_tau); 
                      dta_tau = dta(i_tau,:);
                      ii = sub2ind([s.ns, s.nd], dta_tau.is, dta_tau.id); 
                      %ii = (dta_tau.is-1)*s.nd+dta_tau.id; %through indexes of the chosen alternatives
                      Pi = ccp_tau{tau}(ii);
                      Pi = max(Pi, 1e-5);
                      like{tau} = Pi; 
                  end
                  
                  like_ = [like{:}]; 
                  assert(N_ == numel(like_), 'Unexpected size')
                  
                  N = numel(like_); 
                  % write 
                  ii = im : im+N-1; 
                  moments(ii) = like_(:); 
                  if deriv 
                      error 'not implemented'; 
                  end
                  
                  % update 
                  im = im + numel(ii); 

              case 'mean_carage_tau_j'
                  for tau=1:mp.ntypes
                      for j=1:mp.ncartypes
                          idx=s.ipt.car{j};
                          carage=s.ipt.age(idx)';
                          moments(im)=sum(carage.* h_tau{tau}(idx))/sum(h_tau{tau}(idx));
                          if deriv
                              dmoments(im,:)=sum(carage.* grad_dp_dmp.h_tau{tau}(idx,:),1)/sum(grad_dp_dmp.h_tau{tau}(idx,:),1);
                          end
                          im=im+1;
                      end
                  end
                  
              case 'old_cars'
                  moments(im)=0;
                  for tau=1:mp.ntypes
                      moments(im)=moments(im)+sum(h_tau{tau}(s.ipt.age>=estim.oldcarage));
                      if deriv
                          dmoments(im,:)=dmoments(im,:) + sum(grad_dp_dmp.h_tau{tau}(s.ipt.age>=estim.oldcarage,:),1);
                      end
                  end
                  im=im+1;
                  
              case 'old_cars_tau'
                  moments(im)=0;
                  for tau=1:mp.ntypes
                      moments(im)=sum(h_tau{tau}(s.ipt.age>=estim.oldcarage))/mp.tw(tau);
                      if deriv
                          dmoments(im,:)=sum(grad_dp_dmp.h_tau{tau}(s.ipt.age>=estim.oldcarage,:),1)/mp.tw(tau);
                      end
                      im=im+1;
                  end
                  
              case 'new_cars_tau'
                  moments(im)=0;
                  newcarage=2;
                  for tau=1:mp.ntypes
                      moments(im)=sum(h_tau{tau}(s.ipt.age<=newcarage))/mp.tw(tau);
                      if deriv
                          dmoments(im,:)=sum(grad_dp_dmp.h_tau{tau}(s.ipt.age<=newcarage,:),1)/mp.tw(tau);
                      end
                      im=im+1;
                  end
              otherwise
                  disp('moment not implemented')
          end
      end
      
      if numel(moments) > im-1
          % our preallocation was not large enough 
          moments=moments(1:im-1);
          if deriv 
              dmoments=dmoments(1:im-1,:);
          end
      end
    end

    function [f, g] = obj_moments(mp, s, moments_data, momlist, pvec, moments_w, dta)
      % obj_moments: returns the L2 moments for GMM estimation 
      % FIXME needs optimal weighting

      if nargin < 6 
          w = 1.0; 
      else 
          w = moments_w; 
      end
      
      save('pvec.mat', 'pvec'); 
      
      % 1.update mp struct witb values in pvec
      [mp] = estim.update_mp(pvec, mp);

      % 2. evaluate objective function and derivatives
      if nargout<2
          try 
              moments_model=estim.moments_model(mp, s, momlist, dta);
          catch me 
              warning('estim threw an error, returning inf! Printing pvec and message struct \n'); 
              disp(pvec)
              me
              moments_model = +inf * ones(size(moments_data)); 
          end
      else
        [moments_model, dmoments_model]=estim.moments_model(mp, s, momlist, dta);
        
        I = isnan(moments_data) == false; 
        g = 2*sum(w(I) .* (moments_model(I) - moments_data(I)).*dmoments_model(I, :), 1); 
      end

      I = isnan(moments_data) == false; 
      f = sum(w(I) .* (moments_model(I) - moments_data(I)).^2); 
    end

    function show_model_fit(mp, moments_model, moments_data, moments_w, moments_label, DOLATEX, fid, moments_id, momlist)
      % show_model_fit: FIXME: Add short description here
      % SYNTAX:
      % estim.show_model_fit(mp, moments_model, moments_data, moments_label, DOLATEX, fid)
      % estim.show_model_fit(mp, moments_model, moments_data, moments_label)

      if nargin < 6 
          DOLATEX = false; 
      end
      if nargin < 7
          fid = 1;
      end

      if DOLATEX 
          separator = ' & '; 
          endline = '\\\\ \n'; % double escape to get the "\\" that latex expects
      else
          separator = ' ';
          endline = '\n';
      end

      N = numel(moments_model); 
      printstr_header = sprintf('%%-50s %s %%18s %s %%18s %s', separator, separator, endline);
      printstr_result = sprintf('%%-50s %s %%18.4f %s %%18.4f %s', separator, separator, endline);
      fprintf(fid, printstr_header, 'Moment', 'Model', 'Data');

      if N < 1000 
          % 3-col setup    
          for i=1:N 
              fprintf(fid, printstr_result, moments_label{i}, moments_model(i), moments_data(i)); 
          end

      else
          assert(nargin >= 9, 'Long list of moments, please provide moments-id and momlist for a brief summary')
          mcats = unique(moments_id); 
          for k=1:numel(mcats)
              I = moments_id == mcats(k); 
              w = moments_w(I) ./ mean(moments_w(I)); 
              fprintf(fid, printstr_result, momlist{k}, mean(w .* moments_model(I)), mean(w .* moments_data(I))); 
          end
      end
          
      bar(moments_w .* [moments_model, moments_data]); 
      legend('Model', 'Data'); 
      ylabel('Raw moment value'); 
      xticklabels(moments_label); xtickangle(90); 
      set(gca, 'FontSize', 14); 
    end
    
    function [prob_vec, prob_table] = empirical_choice_probs(dta)
        % INPUTS 
        %   dta: table 
        % 
        % OUTPUT
        %   prob: N-vector of empirical frequencies, sorted the same way as
        %   dta is 
        %   prob_table: same, but in a table including the indexing
        %   variables (and nothing more); 
        
        % compute # of people in cells of state-decision
        % technically, this is a waste since the data should not
        % vary within units of (tau, year, is,id)
        g_ = groupsummary(dta(:, {'tau', 'year', 'is', 'id', 'Count'}), {'tau', 'year', 'is', 'id'}, 'sum');
        g_(:, 'GroupCount') = [];
        
        % compute total # of people at each state
        g__ = groupsummary(g_(:, {'tau', 'year', 'is', 'sum_Count'}), {'tau', 'year', 'is'}, 'sum');
        g__(:, 'GroupCount') = [];
        
        % merge back
        g___ = join(g_, g__);
        
        % Pr(id | is), frequency estimate
        g___.prob = g___.sum_Count ./ g___.sum_sum_Count;
        
        % to ensure they are sorted the same, we merge back on
        prob_table = join(dta(:, {'tau', 'year', 'is', 'id'}), g___(:, {'tau', 'year', 'is', 'id', 'prob'}));
        
        N = size(dta, 1);
        assert(size(prob_table,1) == N, 'Unexpected size; maybe empty cells or otherwise?');
        
        prob_vec = prob_table.prob; 
    end
    
    function print_copypasta_params_for_startvals(mphat, mp_name_string)
      % print copy-paste statement for starting values 
      % prints the estimated parameters in a format suitable to copy up for the
      % starting values directly into the code. 

      if nargin < 2
          mp_name_string = 'mp0'; 
      end

      fprintf('\nHere is a handy snippet that can be copied directly into this script to update starting values\n\n'); 

      for ip=1:numel(mphat.pnames)
          v = mphat.pnames{ip}; 
          if iscell(mphat.(v))
              fprintf('%s.p0.%s = { ', mp_name_string, v); 
              for i_r=1:size(mphat.p0.(v), 1)
                  if i_r>1
                      fprintf('; ');
                  end
                  for i_c=1:size(mphat.p0.(v), 2)
                      if i_c>1
                          fprintf(',');   
                      end
                      fprintf(' %12.8g ', mphat.(v){i_r, i_c}); 
                  end
              end
              fprintf('}; \n'); 
          elseif isscalar(mphat.(v))
              fprintf('%s.p0.%s = %12.8g; \n', mp_name_string, v, mphat.(v)); 
          else 
              error('Unexpected type of variable "%s".', v)
          end
      end
    end
                    
	end % end of methods
end % end of class

