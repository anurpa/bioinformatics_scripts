#Import required modules
from cigar import Cigar
import pandas as pd
import argparse

#Define a function to make a table of trancriptome to genomic position
#mapping.Return genomic positon of query transcript position, from table.

def map_pos(dna_pos,cigar_val,rna_query):
    """
    Return genomic positon of a transcript position.
    
    Args:
        dna_pos: read mapping start position on a chromosome
        cigar_val: cigar string
        rna_query: transcript position
    
    Returns:
        Genomic positon
    """
    #Split cigar using cigar module
    c = Cigar(cigar_val)
    c_split=list(c.items())
    
    #Initiate variables
    rna_pos=0
    dna_pos=dna_pos
    
    #Initiate list
    rna_map=[]
    dna_map=[]
    
    #Using cigar string, build transcript to genomic position mapping table
    for i,(c_len, c_type) in enumerate(c_split):
        
        #Define action for each type of cigar string 
        
        #Cigar type: match, mismatch
        if c_type == "M" or c_type == "=" or c_type == "X":
            rna_map=rna_map+list(range(rna_pos, rna_pos+c_len))
            dna_map=dna_map+list(range(dna_pos, dna_pos+c_len))
            rna_pos=rna_map[-1]+1
            dna_pos=dna_map[-1]+1
        
        #Cigar type: Soft clip
        elif c_type == "S":
            dna_pos=dna_pos-c_len
            rna_map=rna_map+list(range(rna_pos, rna_pos+c_len))
            dna_map=dna_map+list(range(dna_pos, dna_pos+c_len))
            rna_pos=rna_map[-1]+1
            dna_pos=dna_map[-1]+1
        
        #Cigar type: Hard clip    
        elif c_type == "H":
            rna_pos=rna_pos
            dna_pos=dna_pos
        
        #Cigar type: deletion              
        elif c_type == "D":
            rna_map=rna_map+[str(rna_pos)+'D']*c_len
            dna_map=dna_map+list(range(dna_pos, dna_pos+c_len))
            dna_pos=dna_map[-1]+1
        
        #Cigar type: Skipped region in the read    
        elif c_type == "N":
            rna_map=rna_map+[str(rna_pos)+'N']*c_len
            dna_map=dna_map+list(range(dna_pos, dna_pos+c_len))
            dna_pos=dna_map[-1]+1
        
        #Cigar type: insertion in the read                            
        elif c_type == "I" :
            rna_map=rna_map+list(range(rna_pos, rna_pos+c_len))
            dna_map=dna_map+[str(dna_pos)+'I']*c_len
            rna_pos=rna_map[-1]+1
        
        #Cigar type: padding 
        elif c_type == "P" :
            rna_map=rna_map+list(range(rna_pos, rna_pos+c_len))
            dna_map=dna_map+[str(dna_pos)+'P']*c_len
            rna_pos=rna_map[-1]+1
                    
    #Convert list to data frame
    pos_map_df=pd.DataFrame(list(zip(rna_map,dna_map)),columns =['rna','dna'])
    
    #Get genomic position for transcript position query
    dna_val=pos_map_df[pos_map_df['rna']==rna_query]
    #Return genomic position only
    return(dna_val['dna'].values[0])

def translate(args):
    """
    Return table of transcript to genomic mapping
    
    Args:
        genomic_mapping: file with four columns, with genomic mapping
        transcript_coord : two column transcript query file, with transcript name
        and positions
    
    Returns:
        Transcript to genomic mapping
    
    """
    #Exception handling if file not found
    try:
    #Read in tab delimited genomic mapping file, add column names
        file1=pd.read_csv(args.g,delimiter='\t',names=['tr_name','chr_name','chr_pos','cigar_val'])
        file2=pd.read_csv(args.t,delimiter='\t',names=['tr_name','tr_pos'])
    except FileNotFoundError:
        print ("File not found. Check file name and path.")
        # If file not found, exit script
        return
        
    #Merge the two input files
    pos_df=pd.merge(file2,file1,how='outer',on='tr_name')
    
    #Initiate a column for returned genomic position
    pos_df['dna_val']=""
    
    #Iterate each row in table
    for i, row in pos_df.iterrows():
        rna_query=row['tr_pos']
        dna_pos=row['chr_pos']
        cigar_val=row['cigar_val']
        
        #get mapping position from map function
        map_df=map_pos(dna_pos,cigar_val,rna_query)
        #print(map_df)
        pos_df.loc[i,'dna_val']=map_df
        
    #Subset data frame with required columns    
    sub_pos_df=pos_df[['tr_name','tr_pos','chr_name','dna_val']]
    #Return data frame
    print(sub_pos_df)

#To parse arguments from command line
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-g",action="store",dest="g", help="four column file with genomic mapping",required=True)
    parser.add_argument("-t",action="store",dest="t", help="two column file with transcript coordinates",required=True)
    args = parser.parse_args()  
    translate(args)
