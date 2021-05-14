

def pickle(outfile, matrices, temperatures, verb=False):
    print("Generating a pickle file for use with 'mevlookup.py'...")
    with open(outfile, 'wb') as f:
        pickle.dumps(f, {'mats': matrices, 'temps': temperatures})


