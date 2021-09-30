import os
from Bio import SeqIO
from Bio.Seq import Seq 


def retrieve_seq(NM_, mutation_pos=None, folder="snp_downloaded/", save_length_nt=150):
    '''
    NM_:(str) NCBI transcript IDs
    mutation_pos:(str)mutation position
    folder:(str)from which folder 
    save_length_nt:(int) length for extracting the sequence
    Extract sequence from genbank file
    Error code: 'error', 'Mismatch stop codon', 'No stop codon in readthrough'
    '''
    
    n = save_length_nt//3
    
    filename = folder + NM_ + '.gb'
    record = SeqIO.read(filename, "gb")
        
    for feature in record.features: #Coding region
        if feature.type == "CDS":
            location = feature.location
    start = int(location.start)
    end = int(location.end)
    error = 'error'
    try: 
        if mutation_pos!=None:

            mut_pos = int(mutation_pos.split('>')[0].split('.')[1][:-1])
            new_nt = mutation_pos[-1]

            normal = record.seq[start:end] #mRNA coding seq include stop codon 

            if abs(mut_pos - len(record.seq[start:end])) < 3 and mut_pos <= len(record.seq[start:end]):
                readthrough_end = record.seq[start:start+mut_pos-1] + Seq(new_nt) + record.seq[start+mut_pos:] #updated mRNA seq from start to end
            else:
                error = "Mismatch stop codon"
                raise KeyError

            normal_p = normal.translate()
            readthrough_p = readthrough_end.translate()

            if readthrough_p[len(normal_p)-1:].find('*') != -1: #include mutated first stop codon
                readthrough_len = readthrough_p[len(normal_p)-1:].find('*') 
            else:
                error = "No stop codon in readthrough"
                raise KeyError
                #readthrough_len = len(readthrough_p[len(normal_p)-1:]) #if there is no stop codon in readthrough region

            readthrough = readthrough_end[:end-start+3*readthrough_len] #include the second stop codon

            save_normal = normal[-n*3:] 
            save_readthrough = readthrough[-3*(readthrough_len+1):] #include the second and mutated first stop codon
            save_normal_p = save_normal.translate()
            save_readthrough_p = save_readthrough.translate()

        else:
            normal = record.seq[start:end] #mRNA coding seq include stop codon 
            readthrough_end = record.seq[end-3:] #mRNA readthrough seq included first stop codon

            normal_p = normal.translate()
            readthrough_p = readthrough_end.translate()

            if readthrough_p[1:].find('*') != -1: 
                readthrough_len = readthrough_p[1:].find('*')+1 #include first stop codon
            else:
                readthrough_len = len(readthrough_p)

            readthrough = readthrough_end[:3*readthrough_len+3] #include second stop codon

            save_normal = normal[-n*3:]
            save_readthrough = readthrough
            save_normal_p = save_normal.translate()
            save_readthrough_p = save_readthrough.translate()

        return str(save_normal), str(save_normal_p), str(save_readthrough), str(save_readthrough_p)
    except:
        return error