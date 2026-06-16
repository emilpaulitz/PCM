function gpr = simplifyGpr(gpr)
%SIMPLIFYGPR Logically simplify the given gpr using Symbolic Math Toolbox
%   if any gene is replaced by 'false', it will be removed respecting rules

% convert to symbolic friendly string
gpr = regexprep(gpr, ' or ', ' | ', 'ignorecase');
gpr = regexprep(gpr, ' and ', ' & ', 'ignorecase');

% simplify
symGpr = str2sym(gpr);

% convert back to cobra representation
gpr = char(symGpr);

% catch empty gpr
if strcmp(gpr, '[]')
    gpr = char.empty;
else
    gpr = regexprep(gpr, ' \| ', ' or ');
    gpr = regexprep(gpr, ' \& ', ' and ');
end
end

