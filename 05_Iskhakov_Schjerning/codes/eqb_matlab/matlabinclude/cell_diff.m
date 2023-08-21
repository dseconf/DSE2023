  function d12 = cell_diff(c1, c2)
    assert(numel(c1)==numel(c2), 'c1 and c2 must have same number of elements')
    d12=cell(size(c1)); 
    for i=1:numel(c1);
      d12{i}=c1{i}-c2{i};
    end
  end