import os
from Bio import Entrez

def NCBI_downloader(query_ls, folder):
    ''' 
    query_ls: a list of NM_ transcript IDs 
    folder: a string to specify the target folder to download
    Download the genbank files by NM_ transcript IDs to target folder   
    '''
    Entrez.email = "d0789100@gmail.com"  # Always tell NCBI who you are
    files = [file for file in os.listdir(folder)]
    
    for NM_ in query_ls: 
        file = NM_ + '.gb'
        filename = folder + file
        try:
            if file not in files:
                # Downloading...
                input_handle = Entrez.efetch(
                    db="nucleotide", id=NM_, rettype="gb", retmode="text"
                )
                out_handle = open(filename, "w")
                out_handle.write(input_handle.read())
                out_handle.close()
                input_handle.close()
        except:
            print("Interrupted!!There are some problems when handling "+ NM_)