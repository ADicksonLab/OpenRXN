
def join_tuple(tup):
    return '_'.join([str(a) for a in tup])

def make_ID(array_ID,comp_ID):
    if array_ID is not None:
        tag = array_ID + '-'
    else:
        tag = ""
        
    if isinstance(comp_ID, tuple) or isinstance(comp_ID, list):
        return tag + join_tuple(comp_ID)
    elif isinstance(comp_ID, str):
        return tag + comp_ID
    elif isinstance(comp_ID, int):
        return tag + str(comp_ID)
