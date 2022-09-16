#! /usr/bin/python3
# -*- coding: utf-8 -*-

"""
Index.py takes a genome stored in a fasta file and indexed it with FM-Index method.
- Stores the BWT (Burrows Wheeler Transform) of the genome sequence, the F column,
the SA (Suffix Array), the n Dictionnary and the r list.
- Will be used to map reads on the genome.

Author : Elodie GERMANI, Kévin CHATEAU
Version : 21/12/2020
"""
import sys
import os
import getopt
import pickle
import time

import tools_karkkainen_sanders as tks

# Create class which will contain the FMI object
class FMI:
    """
        Stores the FMI index built from the genome.
        Object with 5 attributes : BWT, F, SA, N and R which are calculated
        automatically from the genome.

    """

    def __init__(self, genome):
        """ Defines and stores initial values.

        Parameters
        ----------
        genome : str
            String which will be indexed in the FMI Object.

        Attributes
        ----------
        SA : list
            Suffix array of the genome indexed in the object

        BWT : list
            List containing the second letter of each suffix in alphabetical order
            (Burrows Wheeler Transform)

        f : list
            List containing the first letter of each suffix in alphabetical order

        n : Dictionnary
            Dictionnary containing in keys the alphabetical letter
            and in values the number of occurences of this letter in the genome

        r : list
            List of the same size as the BWT which tell for each caracter of the BWT
            its index (1st A, 1st T...)
        """
        # Suffix array of the genome indexed in the object
        self.sa = tks.simple_kark_sort(genome)

        # BWT of the FMI : list containing the second letter of each suffix in alphabetical order
        self.bwt = [''] * len(self.sa)

        for i in range(len(self.sa)):
            if self.sa[i] == 0:
                self.bwt[i] = "$"
            else:
                self.bwt[i] = genome[self.sa[i]-1]

        # F column of the FMI : aka the first letter of each element of the SA
        self.f = [''] * len(self.sa)

        for i in range(len(self.sa)):
            self.f[i] = genome[self.sa[i]:][0]

        # Dictionnary n containing the number of occurences of each caracter
        # of the alphabet in the genome
        self.n = {}

        # List r of the same size as the BWT which tell for each caracter its index
        # (1st A, 1st T...)
        self.r = []

        self.n["A"] = 0
        self.n["C"] = 0
        self.n["G"] = 0
        self.n["T"] = 0
        self.n["$"] = 0
        for i in range(len(self.bwt)):
            self.n[self.bwt[i]] += 1
            self.r.append(self.n[self.bwt[i]])

    def left_first(self, alpha, k):
        """Return i such as the k-th suffix beginning with an alpha letter is \
            the i-th suffix in the SA

        Parameters
        ----------
        alpha : char,
            letter to find in the FMI
        k : int
            index of the alpha letter to find in the FMI

        Raises
        ------
        ValueError
            if alpha is not in the alphabet. Must be one of ['A', 'T', 'G', 'C', '$']

        Returns
        -------
        i : int
            i such as the k-th suffix beginning with an alpha letter is the i-th suffix in the SA

        """
        if alpha == "$":
            return 0
        elif alpha == "A":
            return int(self.n["$"]) + int(k) - 1
        elif alpha == "C":
            return int(self.n["$"]) + int(self.n["A"]) + int(k) - 1
        elif alpha == "G":
            return int(self.n["$"]) + int(self.n["A"]) + int(self.n["C"]) + int(k) - 1
        elif alpha == "T":
            return int(self.n["$"]) + int(self.n["A"]) + int(self.n["C"]) + int(self.n["G"]) \
                + int(k) - 1
        else:
            raise ValueError("Alpha is not in the alphabet. \
                Must be one of ['A', 'T', 'G', 'C', '$'].")

    def get_down(self, alpha, start, stop):
        """
        Detects the first occurrence of alpha in bwt for i in [start, stop].

        From start go down in the bwt as long as bwt[line] != alpha and line <= stop
          - if bwt[line] == alpha, returns the corresponding line
          - if line > stop: returns -1

        Parameters
        ----------
        alpha : int
            letter to search in the bwt

        start : int
            index where to begin the search in the bwt

        stop : int
            index where to stop the search in the bwt

        Returns
        -------
        line : int
            index where you can find the last occurence of alpha in bwt

        -1
            if alpha not in the bwt[start:stop] list
        """
        line = start
        while line <= stop:
            if self.bwt[line] == alpha:
                return line
            line += 1
        return -1

    def get_up(self, alpha, start, stop):
        """
        Detects the first occurrence of alpha in bwt for i in [start, stop].

        From stop go up in the bwt as long as bwt[line] != alpha and line >= stop
          - if bwt[line] == alpha, returns the corresponding line
          - if line > stop: returns -1

        Parameters
        ----------
        alpha : int
            letter to search in the bwt

        start : int
            index where to begin the search in the bwt

        stop : int
            index where to stop the search in the bwt

        Returns
        -------
        line : int        # F column of the FMI : aka the first letter of each element of \
        the SA
        self.f = [''] * len(self.sa)

        for i in range(len(self.sa)):
            self.f[i] = genome[self.sa[i]:][0]
        -1
            if alpha not in the bwt[start:stop] list
        """
        line = stop
        while line >= start:
            if self.bwt[line] == alpha:
                return line
            line -= 1
        return -1

    def get_occurences(self, pattern):
        """
        Returns the list of positions where the pattern P can be mapped to in the object.

        Parameters
        ----------
        pattern : str
            pattern to search in the FMI

        Raises
        ------
        ValueError
            if P contains characters that are not in the alphabet.

        Returns
        -------
        list_occu : list
            list of position index of occurences of the pattern in the genome represented by \
                the FMI
        """

        #start = 0
        #stop = len(self.bwt)-1
        last_char = pattern[-1]

        if last_char not in ['A', 'T', 'G', 'C', '$']:
            raise ValueError("P contains characters that are not in the alphabet used here. \
                Must be one of ['A', 'T', 'G', 'C', '$'].")

        o = self.n[last_char]
        start = self.left_first(last_char, 1)
        stop = self.left_first(last_char, o)
        # read the pattern from right to left
        for pos_pattern in range(len(pattern)-2, -1, -1):
            current_char = pattern[pos_pattern]

            if current_char not in ['A', 'T', 'G', 'C', '$']:
                raise ValueError("P contains characters that are not \
                    in the alphabet used here. Must be one of ['A', 'T', 'G', 'C', '$'].")

            new_start = self.get_down(current_char, start, stop)
            if new_start == -1:
                return []
            new_stop = self.get_up(current_char, start, stop)
            start = self.left_first(self.bwt[new_start], self.r[new_start])
            stop = self.left_first(self.bwt[new_stop], self.r[new_stop])
        list_occu = []
        for i in range(start, stop+1):
            list_occu.append(self.sa[i])
        return list_occu

if __name__ == "__main__":
    # Create variables which will contain options input by the user
    REFERENCE_FILE = None
    OUTPUT_FILE = None
    READ_FILES = None
    VERBOSE = False

    try:
        OPTIONS, REMAINDER = getopt.getopt(sys.argv[1:], 'o:r:vh', ['out=',
                                                             'ref=',
                                                             'verbose',
                                                             'help'
                                                             ])
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
        elif opt in ('-h', '--help'):
            os.system('vi README.md')
            sys.exit(0)

    if VERBOSE:
        print('OPTIONS   :', OPTIONS)
        print('VERBOSE   :', VERBOSE)
        print('OUTPUT    :', OUTPUT_FILE)
        print('REFERENCE :', REFERENCE_FILE)

    # Reading of the reference file only if REFERENCE_FILE and OUTPUT_FILE are filled out
    if REFERENCE_FILE is not None and OUTPUT_FILE is not None:
        try:
            T1 = time.time()
            GENOME_FILE = open(REFERENCE_FILE, 'r')
            GENOME = GENOME_FILE.readlines()[-1]
            GENOME_FMI = FMI(GENOME)
            pickle.dump(GENOME_FMI, open(OUTPUT_FILE, "wb"))
            T2 = time.time()
            if VERBOSE:
                print(f"{GENOME_FILE} indexed with success in {OUTPUT_FILE}. \
                    \nTime of execution : {T2 - T1} seconds.")

        except OSError as error:
            print(error)
            sys.exit(2) # If file not found, raise an error and exit program.
