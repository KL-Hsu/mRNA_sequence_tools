from Bio.Seq import Seq 
import pyensembl as pen

genome = pen.EnsemblRelease(102)

def position_of_codon(transcript, stop=True): 
    
    if stop == False:
        pos = transcript.stop_codon_positions[-1]
        distance = len(transcript.five_prime_utr_sequence)
    else:
        intervals = transcript.coding_sequence_position_ranges
        n = 0 
        for i in intervals:
            n+= (i[1] - i[0]+1)
        distance = len(transcript.five_prime_utr_sequence) + n
    
    return distance 


def retrieve_seq_from_ensemble(ID, mutation_pos, save_length_nt=150):
    '''
    ID: (string)Ensembl transcript ID, mutation_pos: (string) c.1452T>C
    save_length_nt: default 150 nt
    Extract sequneces from Ensembl by given mutations
    Error code: No transcript, No cds!, Mismatch stop codon, No stop codon in readthrough 
    '''
    n = save_length_nt//3
    mut_pos = int(mutation_pos.split('>')[0].split('.')[1][:-1]) 
    old_nt = mutation_pos[-3]
    new_nt = mutation_pos[-1] 

    error = "No transcript"
    try:
        transcript = genome.transcript_by_id(ID)

        try:
            mRNA_cds = transcript.coding_sequence
        except:
            error = "No cds!"
            raise KeyError     

        save_mRNA = mRNA_cds[-3*n:]
        save_protein = Seq(save_mRNA).translate()

        start = position_of_codon(transcript, stop=False)
        stop = position_of_codon(transcript, stop=True)

        if abs(mut_pos - len(transcript.coding_sequence)) <= 3 and mut_pos <= len(transcript.coding_sequence):

            updated_cds = transcript.coding_sequence[:mut_pos-1] + new_nt + transcript.coding_sequence[mut_pos:]
            mutated_stop = updated_cds[-3:]
        else: 
            error = "Mismatch stop codon"
            raise KeyError

        UTR_3 = transcript.three_prime_utr_sequence
        readthrough = mutated_stop + UTR_3 
        readthrough_p = Seq(readthrough).translate()

        index = readthrough_p.find('*')
        if index != -1:
            save_readthrough = readthrough[:(index+1)*3]
            save_readthrough_p = readthrough_p[:index+1]
        else:
            error = "No stop codon in readthrough"
            raise KeyError

        return str(save_mRNA), str(save_protein), str(save_readthrough), str(save_readthrough_p)

    except:
        return error

def retrieve_normal_seq_from_ensemble(ID, save_length_nt=150):

    n = save_length_nt//3
    error = "No transcript"
    try:
        transcript = genome.transcript_by_id(ID)

        try:
            mRNA_cds = transcript.coding_sequence
        except:
            error = "No cds!"
            raise KeyError     

        save_mRNA = mRNA_cds[-3*n:]
        save_protein = Seq(save_mRNA).translate()

        return str(save_mRNA), str(save_protein)

    except:
        return error