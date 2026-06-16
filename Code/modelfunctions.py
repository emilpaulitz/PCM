import cobra
import sys
import io
from cobra import Reaction

class silence():
    def __enter__(self):
        self._stdout = sys.stdout
        self._stderr = sys.stderr

        self.text_trap = io.StringIO()
        sys.stdout = self.text_trap
        sys.stderr = self.text_trap
        return self

    def __exit__(self, *args):
        sys.stdout = self._stdout
        sys.stderr = self._stderr

def get_gid(model, gid):
    return model.genes.get_by_id(gid)

def get_mid(model, mid):
    return model.metabolites.get_by_id(mid)

def get_rid(model, rid):
    return model.reactions.get_by_id(rid)

# extract the reaction ID that the given ecGEM rid would be in a conventional model
def get_ori_rid(rid):
    new_rid = rid
    if '_EXP_' in rid:
        new_rid = rid[:rid.find('_EXP_')]
    return new_rid.removesuffix('_REV')
    
def find_group(model, r):
    if type(r) == str:
        r = get_rid(model, r)
    
    gs = list()
    for g in model.groups:
        if r in g.members:
            gs.append(g)

    return gs

# return_val can be one of:
#   'flux': only return number of maximum flux
#   'sol': return cobraPy solution object
#   'model': return solved model object
def check_production(model, mid, add_export = list(), add_import = list(), exclude_rxns = list(), 
                     consumption = False, return_val = 'flux', perform_pfba = True,
                     import_ub = 1000, export_lb = 0):
    if type(mid) == str:
        target_met = get_mid(model, mid)
    with model:
        for rid in exclude_rxns:
            get_rid(model, rid).knock_out()

        exp_rxns = list()
        for exp in add_export:
            if type(export_lb) == dict:
                curr_lb = export_lb[exp] if exp in export_lb else 0
            else:
                curr_lb = export_lb
            rxn = Reaction('export' + exp, lower_bound = curr_lb)
            rxn.add_metabolites({get_mid(model, exp): -1})
            exp_rxns.append(rxn)
        model.add_reactions(exp_rxns)

        imp_rxns = list()
        for imp in add_import:
            if type(import_ub) == dict:
                curr_ub = import_ub[imp] if imp in import_ub else 1000
            else:
                curr_ub = import_ub
            rxn = Reaction('import' + imp, upper_bound = curr_ub)
            rxn.add_metabolites({get_mid(model, imp): 1})
            imp_rxns.append(rxn)
        model.add_reactions(imp_rxns)

        biomass_reaction = Reaction('BIOMASS_tmp')
        biomass_reaction.name = 'BIOMASS_tmp'
        biomass_reaction.add_metabolites({target_met: 1 if consumption else -1})
        model.add_reactions([biomass_reaction])
        model.objective = 'BIOMASS_tmp'
        if return_val == 'sol':
            result = cobra.flux_analysis.pfba(model) if perform_pfba else model.optimize()
        elif return_val == 'model':
            if perform_pfba:
                cobra.flux_analysis.pfba(model)
            else:
                model.optimize()
            result = model.copy()
        else:
            result = model.slim_optimize()
    
    return result