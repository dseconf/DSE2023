function struct_diff = struct_diff(sol_d,sol_m)
  % struct_diff = struct_diff(sol_d, sol_m)
  struct_diff=sol_d; 
  fld=fields(sol_d); 
  for i=1:numel(fld)
    if iscell(sol_d.(fld{i}))
      struct_diff.(fld{i})=cell_diff(sol_d.(fld{i}), sol_m.(fld{i})); 
    else
      struct_diff.(fld{i})=sol_d.(fld{i}) - sol_m.(fld{i}); 
    end
  end
end