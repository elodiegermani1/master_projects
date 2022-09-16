#! /usr/bin/python3
# -*- coding: utf-8 -*-
"""
This program maps reads stored in a fasta file on a genome indexed \
    with FM-Index method (see index.py).

Author : Elodie GERMANI, Kévin CHATEAU 
Version : 21/12/2020
"""
import sys
import os
import getopt
import pickle
import time

from index import FMI

def seed_and_extend(read, k, h, index, genome):
    """ Try to seed the kmer into the genome file indexed by our FMI \
        object.

    Parameters
    ----------
    read : str
        Read we try to map on the genome.

    k : int
        Size of the kmers to seed.

    h : int
        Maximum value of mismatch in the mapping.

    index : FMI
        FMI which stores the index of the genome in which we want to map \
        the read

    genome : str
        Sequence indexed in the FMI Object

    Returns
    -------
    position_mapping : int
        Index on the genome where the best mapping is obtained with this read.\
        Length of the genome if the read cannot be mapped on the genome with \
        these parameters.

    nb_mismatch : int
        Number of mismatches in the alignment. h if no mapping found.

    list_mismatch : listl
        List of the position on the genome of the mismatches detected during \
        this alignment. Empty list if no mapping found or no mismatch.

    """

    list_mapping_read = [] # List containing the positions tested to map the read on the genome
    #(will be used to not try to align a read twice at the same position)

    # Variables which will be returned
    position_mapping = len(genome) # Optimal position of mapping for the read
    nb_mismatch = int(h) + 1 # Number of mismatch in this mapping
    list_mismatch = [] # List of mismatch positions on the genome

    for kmer_index in range(len(read)-int(k)+1):
        kmer = read[kmer_index:kmer_index + int(k)]
        # For each kmer, tries to find the optimal position of mapping
        # for the read with this kmer as seed.
        position_mapping_kmer = len(genome)
        nb_mismatch_kmer = int(h) + 1
        list_mismatch_kmer = []

        list_occurences = sorted(index.get_occurences(kmer))

        if not list_occurences:
            continue

        for occurences in list_occurences:

            nb_mismatch_occu = 0 # For each occurence of the kmer,
            # count the number of mismatch during alignment

            list_mismatch_occu = [] # List of mismatch seen during alignment
            # of read with this occurence of the kmer

            index_char_genome = occurences - kmer_index # Index where to map in the genome
            index_char_read = 0 # Index of the character to compare

            if index_char_genome in list_mapping_read: # If position already tested,
                #do not test it a second time.
                continue
            else:
                list_mapping_read.append(index_char_genome) # Add this position to the list
                # so it won't be tested a second time for this read

            while nb_mismatch_occu <= int(h) \
                 and index_char_read < len(read) \
                      and index_char_genome < len(genome):
                if genome[index_char_genome] != read[index_char_read]:
                    nb_mismatch_occu += 1
                    list_mismatch_occu.append(index_char_genome)

                index_char_genome += 1
                index_char_read += 1


            # If the mapping of the read with this occurence of the read
            # is better than the previous one (less mismatch) : optimal values for kmer stored
            if nb_mismatch_occu < nb_mismatch_kmer:
                nb_mismatch_kmer = nb_mismatch_occu
                list_mismatch_kmer = list_mismatch_occu
                position_mapping_kmer = occurences - kmer_index

        # If the best mapping found for this kmer is better than the mapping
        # found with the previous kmer : optimal values for read stored
        if nb_mismatch_kmer < nb_mismatch \
             or nb_mismatch_kmer == nb_mismatch \
                  and position_mapping_kmer < position_mapping:
            nb_mismatch = nb_mismatch_kmer
            list_mismatch = list_mismatch_kmer
            position_mapping = position_mapping_kmer

    return position_mapping, nb_mismatch, list_mismatch



def reverse_read(read):
    """
    Creates the equivalent of this read on the reverse strand.
    Modify "A" to "T" and inverse and "C" to "G" and inverse.

    Parameters
    ----------
    read : str
        Read to reverse.

    Raises
    ------
    ValueError
        If read contains characters not in the alphabet

    Returns
    -------
    reversed_read : str
        Read reversed.
    """
    reversed_read = ""
    for i in range(len(read)-1, -1, -1):
        if read[i] == "A":
            reversed_read += "T"
        elif read[i] == "T":
            reversed_read += "A"
        elif read[i] == "G":
            reversed_read += "C"
        elif read[i] == "C":
            reversed_read += "G"
        else:
            raise ValueError("One of the read contains wrong characters.")

    return reversed_read


def mapping(reads_list, k, h, index, genome):
    """ Maps the reads stored in reads_list on the genome using seed_and_extend function
        and creates a dictionary which will contains the snps found during mapping :
     - keys will be the index in the genome where the snps was found
     - values will be the list of characters found in reads mapped at this position.

    Parameters
    ----------
    reads_list : list of str
        list of reads to be mapped

    K_VALUE : int
        Size of the kmers to seed.

    H_VALUE : int
        Maximum value of mismatch in the mapping.

    index : FMI
        FMI which stores the index of the genome in which we want to map the read

    genome : str
        Sequence indexed in the FMI Object

    Returns
    -------
    snps_dict : Dictionary
        dictionary which will contains the snps found during mapping :
            - keys will be the index in the genome where the snps was found
            - values will be the list of characters found in reads mapped at this position.
    """
    snps_dict = {}
    # Map the read on the genome and store the snps found
    for read in reads_list:
        reversed_read = reverse_read(read)
        reverse = False
        list_mapping = seed_and_extend(read, k, h, index, genome)
        if list_mapping[0] < len(genome):
            reverse = False
            if VERBOSE:
                print("Read number : ", reads_list.index(read) + 1, \
                    "\n Mapping at position :", list_mapping[0], \
                        " on straight strand. \n With ", list_mapping[1], \
                            "substitutions at positions :", list_mapping[2])
        else:
            list_mapping = seed_and_extend(reversed_read, k, h, index, genome)
            if list_mapping[0] < len(genome):
                reverse = True
                if VERBOSE:
                    print("Read number : ", reads_list.index(read) + 1, \
                        "\n Mapping at position :", list_mapping[0], \
                            " on reverse strand. \n With ", list_mapping[1], \
                                "substitutions at positions :", list_mapping[2])
            else:
                reverse = False
                if VERBOSE:
                    print("No mapping found for read number :", reads_list.index(read) + 1)
        if list_mapping[0] < len(genome):
            for mismatch in list_mapping[2]:
                if reverse == False:
                    if mismatch in snps_dict.keys():
                        snps_dict[mismatch].append(read[mismatch - list_mapping[0]])
                    else:
                        snps_dict[mismatch] = [read[mismatch - list_mapping[0]]]
                else:
                    if mismatch in snps_dict.keys():
                        snps_dict[mismatch].append(reversed_read[mismatch - list_mapping[0]])
                    else:
                        snps_dict[mismatch] = [reversed_read[mismatch - list_mapping[0]]]

    return snps_dict


def write_vcf(snps_dict):
    """ Write vcf file containing the list of Single Nucleotid Polymorphisms found during mapping

    Parameters
    ----------
    snps_dict : Dictionnary
        dictionnary which contains the snps found during mapping :
            - keys will be the index in the genome where the snps was found
            - values will be the list of characters found in reads mapped at this position.

    Returns
    -------
    No return but write a VCF file containing the list of Single Nucleotid Polymorphisms found during mapping.

    VCF file example :
    #READS: my_reads.fa
    #K: 14
    #MAX_SUBST: 5
    #MIN_ABUNDANCE: 10
    129836 A T 26
    145831 C G 25

    """
    # Header of the vcf file
    header = f"""#REF: {REFERENCE_FILE}
#READS: {READS_FILE}
#K: {K_VALUE}
#MAX_SUBST: {H_VALUE}
#MIN_ABUNDANCE: {M_VALUE}
"""

    with open(OUTPUT_FILE, 'w') as vcf:
        vcf.write(header)
        for position in sorted(snps_dict.keys()): # For each snp position found,
            # count for each nucleotid the number of time it was found in reads mapped
            # at this position
            nA = 0
            nT = 0
            nC = 0
            nG = 0
            for nucleotid in snps_dict[position]:
                if nucleotid == "A":
                    nA += 1
                elif nucleotid == "T":
                    nT += 1
                elif nucleotid == "G":
                    nG += 1
                else:
                    nC += 1
            if nA >= int(M_VALUE): # If the same nucleotid was found more than M_VALUE time
                # in reads mapped at this position, write it in the vcf file.
                vcf.write(f"{position}\t{GENOME[position]}\tA\t{nA}\n")
            if nT >= int(M_VALUE):
                vcf.write(f"{position}\t{GENOME[position]}\tT\t{nT}\n")
            if nG >= int(M_VALUE):
                vcf.write(f"{position}\t{GENOME[position]}\tG\t{nG}\n")
            if nC >= int(M_VALUE):
                vcf.write(f"{position}\t{GENOME[position]}\tC\t{nC}\n")


if __name__ == "__main__":
    # Create variables which will contain options input by the user
    REFERENCE_FILE = None
    OUTPUT_FILE = None
    INDEX_FILE = None
    READS_FILE = None
    K_VALUE = None
    H_VALUE = None
    M_VALUE = None
    VERBOSE = False

    try:
        OPTIONS, REMAINDER = getopt.getopt(sys.argv[1:], \
            'o:r:i:k:m:vh', ['out=', 'ref=', 'index=', 'reads=', 'K_VALUE=', 'max_hamming=', 'min_abundance=', 'VERBOSE', 'help'])

    except getopt.GetoptError as err:
        print(err)
        sys.exit(2)

    # Replace variables depending on options
    for opt, arg in OPTIONS:
        if opt in ('-o', '--out'):
            OUTPUT_FILE = arg
        elif opt in ('-r', '--ref'):
            REFERENCE_FILE = arg
        elif opt in ('-v', '--verbose'):
            VERBOSE = True
        elif opt in ('--reads'):
            READS_FILE = arg
        elif opt in ('-i', '--index'):
            INDEX_FILE = arg
        elif opt in ('-k', '--K_VALUE'):
            K_VALUE = arg
        elif opt in ('--max_hamming'):
            H_VALUE = arg
        elif opt in ('-m', '--min_abundance'):
            M_VALUE = arg
        elif opt in ('-h', '--help'):
            os.system("vi README.md")
            sys.exit()

    print('OPTIONS   :', OPTIONS)

    if VERBOSE:
        print('VERBOSE   :', VERBOSE)
        print('OUTPUT    :', OUTPUT_FILE)
        print('REFERENCE :', REFERENCE_FILE)
        print('READS :', READS_FILE)
        print('INDEX :', INDEX_FILE)
        print('K VALUE :', K_VALUE)
        print('MAX HAMMING :', H_VALUE)
        print('MIN ABUNDANCE :', M_VALUE)
        print('REMAINING :', REMAINDER)

    # If all options are filed -> execute program.
    if REFERENCE_FILE is not None \
        and OUTPUT_FILE is not None \
            and READS_FILE is not None \
                and INDEX_FILE is not None \
                    and K_VALUE is not None \
                        and M_VALUE is not None \
                            and H_VALUE is not None:
        try:
            # Read genome stored in REFERENCE_FILE
            GENOME_FILE = open(REFERENCE_FILE, 'r')
            GENOME = GENOME_FILE.readlines()[-1]
            # Read reads stored in READS_FILE
            READS_TEXT = open(READS_FILE, "r")
            READS_LIST = []
            for line in READS_TEXT.readlines():
                if line[0] != ">":
                    READS_LIST.append(line[:-1])
            # Read index stored in INDEX_FILE
            INDEX = pickle.load(open(str(INDEX_FILE), "rb"))

            T1 = time.time()

            SNPS_DICT = mapping(READS_LIST, K_VALUE, H_VALUE, INDEX, GENOME)

            T2 = time.time()

            write_vcf(SNPS_DICT)

            if VERBOSE:
                print(f"The mapping of {READS_FILE} on {REFERENCE_FILE} using {INDEX_FILE} as index took {T2-T1} second(s)")

        except OSError as error:
            print(error)
            sys.exit(2) # If file not found, raise an error and exit program.
