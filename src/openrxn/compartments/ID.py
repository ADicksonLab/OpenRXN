
def join_tup(tup):
    return '_'.join([str(a) for a in tup])

def makeID(array_ID,comp_ID):
    if array_ID is not None:
        tag = array_ID + '-'
    else:
        tag = ""
        
    if type(comp_ID) is tuple or type(comp_ID) is list:
        return tag + join_tup(comp_ID)
    elif type(comp_ID) is str:
        return tag + comp_ID
    elif type(comp_ID) is int:
        return tag + str(comp_ID)
