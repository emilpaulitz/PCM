import cobra
import sys
import io

def get_rid(model, rid):
    return model.reactions.get_by_id(rid)
    
def find_group(model, r):
    if type(r) == str:
        r = get_rid(model, r)
    
    gs = list()
    for g in model.groups:
        if r in g.members:
            gs.append(g)

    return gs

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