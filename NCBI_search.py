from Bio import Entrez

def NCBI_esearch(db="nucleotide", term="(RAD1[Gene Name]) AND Human[Organism]"):
    '''
    db:(str) database default is nucleotide
    term: (str) keywords for searching default is (RAD1[Gene Name]) AND Human[Organism]
    Searching by keywords and return a record
    '''

    Entrez.email = "your_email_address"  # Always tell NCBI who you are

    # Downloading...
    handle = Entrez.esearch(db=db, term=term, retmax=10, idtype='acc')
    record = Entrez.read(handle)
    handle.close()

    return record 
