function model = addCytosolMet(model,chloroMetId, newComp, AnnoFieldList)
%ADDCYTOSOLMET adds a cytosol version of the given chloroplast metabolite
%to the model. Copies all annotations from the given chloroplas met.

if nargin < 3
    newComp = 'c';
end
if nargin < 4
    AnnoFieldList = {'metPubChemID', 'metKEGGID', 'metChEBIID',...
        'metMetaNetXID', 'metBioCycID'};
end
% find index
chlMetIx = strcmp(model.mets, chloroMetId);

% get met ID without [h]
if length(chloroMetId) >= 3 && strcmp(chloroMetId(end-2:end), '[h]')
    metID = chloroMetId(1:end-3);
else
        % Return the original string unchanged if suffix is not present
    metID = chloroMetId;
end
newID = strcat(metID, ['[' newComp ']']);

model = addMetabolite(model, newID, 'metName', ...
    model.metNames(chlMetIx), 'metFormula', model.metFormulas(chlMetIx), ...
    'Charge', model.metCharges(chlMetIx));
for fieldNameCell = AnnoFieldList
    fieldName = char(fieldNameCell);    
    model.(fieldName)(strcmp(model.mets, newID)) = ...
        model.(fieldName)(chlMetIx);
end

end

