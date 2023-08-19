function pstruct = vec2struct(pvec,pnames, pstruct);
    %VEC2STRUCT: Procedure to convert a vector to a structure with scalar fields. 
    %
    %SYNTAX:      
    %   pstruct = estim.vec2struct(pvec,pnames, pstruct)
    %
    %OUTPUTS:
    %   pstruct:   Structure with fields taking scalar values (at least one for each value element in pnames) 
    %
    %INPUTS
    %   pvec:      k dimentional vector of scalars
    %
    %   pnames:    k dimensional cell array, that holds names of the fields corresponding to each row in pvec
    %
    %   pstruct:   
    %              If this argument is passed, the fields pstruct corresponding to pnames will be be updated 
    %              with the values in pvec (if they exist already) or added if they do not already exists. 
    %
    % See also: estim.struct2vec
    k=numel(pnames);
    j = 1;
    for i=1:k;
        r_i = size(pstruct.(char(pnames(i))),1);
        c_i = size(pstruct.(char(pnames(i))),2);
        for c=1:c_i;
            for r=1:r_i;
                pstruct.(char(pnames(i)))(r, c)=pvec(j,1);
                j = j + 1; 
            end
        end    
    end
end % end of estim.vec2struct
