

from pickle import dump


def pickle(outfile, matrices, temperatures, verb=False):
    print("Generating a pickle file for use with 'mevlookup.py'...")
    with open(outfile, 'wb') as f:
        dump({'mats': matrices, 'temps': temperatures}, f)


