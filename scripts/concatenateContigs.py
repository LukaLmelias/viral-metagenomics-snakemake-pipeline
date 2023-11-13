import os
import argparse

"""
contig_concatenation Module

This module contains functions to concatenate contig files from (meta)spades directories and rename contig headers.

Functions:
- concatenate_fasta(spades_path, all_contigs_filename, sample_ids): Concatenate contig files and rename headers.
- main(): Entry point for the script to run contig concatenation.
"""

#pick samples from metaspades directory

def concatenate_fasta(spades_path, all_contigs_filename,sample_ids):

    """
    Concatenate contig files from (meta)spades directories and rename contig headers.

    Parameters:
    spades_path (str): Path to the (meta)spades directory containing sample subdirectories.
    all_contigs_filename (str): Base name for the file with all concatenated contigs.
    sample_ids (str): Base name for the file containing sample IDs assigned to shorten their names.

    Returns:
    None: This function saves concatenated contigs and sample ID information as files.
    """
    samples = os.listdir(spades_path)


    # assign samples id to shorten their names
    samplesIDs = {}
    for index, sample in enumerate(samples):
        samplesIDs[f'sample_{index + 1}']= sample
    with open( sample_ids, "w") as txt:
        txt.write(str(samplesIDs))

    #rename all contig headers to sample_id_original_contig_name
    with open(f"{ all_contigs_filename}", "w") as fn:
        for id,  sample in zip(samplesIDs.keys(), samples):
            contigPath = os.path.join(spades_path, sample, "contigs.fasta")

            with open(contigPath, "r") as contigs:
                for line in contigs:                
                    if line.startswith(">"):
                        #print(line)

                        newHeader = f'>{id}_{line[1:]}'# replace the header with sample info
                        # print(newHeader)
                        fn.write(newHeader)
                    else:
                        fn.write(line)
                    
    return None
        
def main():
    """
    Entry point for the script to run contig concatenation.

    Usage: python script.py spades_path all_contigs_filename sample_id_filename
    """

    #parse paths
    parser = argparse.ArgumentParser()
    parser.add_argument("spades_path", help="path to the (meta)spades directory")
    parser.add_argument("all_contigs_filename", help="Base name for the file with all contigs")
    parser.add_argument("sample_id_filename", help="Base name for the file ids assigned to samples")
    args = parser.parse_args()


    metaspadesPath = args.spades_path
    all_contigs_filename = args.all_contigs_filename
    sample_ids = args.sample_id_filename

    #run concate function
    concatenate_fasta(metaspadesPath, all_contigs_filename, sample_ids)

if __name__ == "__main__":
    main()

    
   
    