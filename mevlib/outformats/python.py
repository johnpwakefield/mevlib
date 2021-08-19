

from pickle import dump


#TODO should really be an `all' since the pickle format is pretty flexible

def write_mat_pkl(outfile, spcsyms, temperatures, matrices, verb=False):
    print("Generating a pickle file for use with 'mevlookup.py'...")
    with open(outfile, 'wb') as f:
        dump({'mats': matrices, 'temps': temperatures}, f)


