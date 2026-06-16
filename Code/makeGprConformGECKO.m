function gpr = makeGprConformGECKO(gpr)
%makeGprConformGECKO Re-phrase the gpr to be of the form 
% (U1 and U2 and ...) or (...) or (...) so GECKO can read the gpr correctly

% convert to symbolic friendly string
gpr = regexprep(gpr, ' or ', ' | ', 'ignorecase');
gpr = regexprep(gpr, ' and ', ' & ', 'ignorecase');

% expand and export to char. No brackets but good enough
symGpr = expand(str2sym(gpr));
gpr = char(symGpr);

% catch empty gpr
if strcmp(gpr, '[]')
    gpr = char.empty;
else
    gpr = regexprep(gpr, ' \| ', ' or ');
    gpr = regexprep(gpr, ' \& ', ' and ');
end
end

