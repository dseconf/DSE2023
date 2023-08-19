function pstruct1 = getfields(pstruct0, fields);
    %GETFIELDS: Procedure to select a subset of fields of a structure and output them is a new structure.  
    %
    %SYNTAX:      
    %   pstruct1 = getfields(pstruct0, fields)
    %
    %OUTPUTS:
    %   pstruct1:   Structure with k fields corresponding to fieldnames given in the input "fields"
    %
    %INPUTS
    %   pstruct0:   
    %              Input structure with l>k fields
    %
    %   fields:    k dimensional cell array, that holds names of the fields to be selected from pstruct0
    %
    % See also: struct2vec, vec2struct
    k=numel(fields);
    for i=1:k;
     pstruct1.(char(fields(i)))=pstruct0.(char(fields(i)));
    end
end % end of vec2struct