function [structure_out] = append_structure(struture,string,value)
test = struct(string,value);
structure_out = catstruct(struture,test);
end

