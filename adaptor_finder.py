#import required modules
from Bio import SeqIO
import argparse

#Define a function to count frequencies of substrings in a FASTQ file
def find_adaptors(args): 
    """
    Get frequencies of substrings of length k across all sequences in a 
    fastq file(f).
    
    Args:
        f: FASTQ file
        k: length of substring
        
    Returns:
        Substrings and counts, sorted by counts
    """
    #Initiate a dictionary, to hold substrings as key and count as values
    kmers={}
    #Loop over each sequence in fastq file
    for record in SeqIO.parse(args.f, "fastq"):   
        #Retrieve and save raw sequence line as variable seq
        seq=record.seq
        #Pull out substring of required length from sequence
        for i in range(len(seq) - args.k + 1):
            kmer = seq[i:i+args.k]
            #If substring already exists, increase count
            if kmer in kmers:
                kmers[kmer] += 1
            #If substring is new, store new key
            else:
                kmers[kmer] = 1 
    #Print substrings and counts, sorted by values(counts)
    for s in sorted(kmers, key=kmers.get, reverse=True):
        print (s, kmers[s])

#Parse arguments from command line
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-k",action="store",dest="k",type=int, help="length of substring",required=True)
    parser.add_argument("-f",action="store",dest="f", help="fasta file",required=True)
    args = parser.parse_args()  
    find_adaptors(args)
