function struct_out = concatStruct(struct1,struct2,concatdim);
%function struct_out = concatStruct(struct1,struct2,concatdim);
if nargin<3
    concatdim = 1;
end

fields1 = fieldnames(struct1);
fields2 = fieldnames(struct2);

if ~all(ismember(fields1,fields2))
    error('Structure fields don''t match!'); end

struct_out = struct;
for ifield = 1:length(fields1)
    struct_out.(fields1{ifield}) = cat(concatdim, struct1.(fields1{ifield}), struct2.(fields1{ifield}));
end
end