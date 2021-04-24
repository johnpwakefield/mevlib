

from mevlib.parsing.sbl import parse_sensible


# file parsing

# check if is zero; ideally this will be rewritten to be type-agnostic
def iszero(x):
    return x == 0

def parse_chemkin(fn):
    raise NotImplementedError("Chemkin parser not written.")

def parse_ini(fn):
    raise NotImplementedError("Ini parser not written.")

file_types = {
#   '.inp'  : ("Chemkin-II", parse_chemkin),
#   '.ini'  : ("ini", parse_ini),
    '.sbl'  : ("sensible", parse_sensible)
}

def parse_attempt(f, ext, verb):
    name, parser = file_types[ext]
    try:
        if verb:
            print("Attempting to parse input as a {} file.".format(name))
        params = parser(f, verb=verb)
        if verb:
            print("Input successfully parsed as a {} file.".format(name))
        return params
    except (FileTypeException,) as err:
        print("File could not be parsed as a {} file.".format(name))
        print(err)
        return None

def parse_dynamic(f, verb):
    ext = f[:-4] if f[:-4] in file_types.keys() else None
    if ext is not None:
        if verb:
            print("Inferring file type from file extension.")
        params = parse_attempt(f, ext, verb)
        if params is not None:
            return params
    for k in file_types.keys():
        if k != ext:
            params = parse_attempt(f, k, verb)
            if params is not None:
                return params
    print("Unable to parse input file.")
    raise Exception("Unable to parse input file.")
