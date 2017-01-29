# -*- coding: utf-8 -*-
"""
YOUR HEADER COMMENT HERE

@author: Seungin Lyu

"""

import random
from amino_acids import aa, codons, aa_table   # you may find these useful
from load import load_seq


def shuffle_string(s):
    """Shuffles the characters in the input string
        NOTE: this is a helper function, you do not
        have to modify this in any way """
    return ''.join(random.sample(s, len(s)))

# YOU WILL START YOUR IMPLEMENTATION FROM HERE DOWN ###


def get_complement(nucleotide):
    """ Returns the complementary nucleotide

        nucleotide: a nucleotide (A, C, G, or T) represented as a string
        returns: the complementary nucleotide
    >>> get_complement('A')
    'T'
    >>> get_complement('C')
    'G'

    #checks the complements backward
    >>> get_complement(get_complement('A'))
    'A'
    >>> get_complement(get_complement('C'))
    'C'
    """
    complementary_nucelotide = ''

    if (nucleotide == 'A'):
        complementary_nucelotide = 'T'
    elif (nucleotide == 'T'):
        complementary_nucelotide = 'A'
    elif (nucleotide == 'G'):
        complementary_nucelotide = 'C'
    elif (nucleotide == 'C'):
        complementary_nucelotide = 'G'

    return complementary_nucelotide


def get_reverse_complement(dna):
    """ Computes the reverse complementary sequence of DNA for the specfied DNA
        sequence

        dna: a DNA sequence represented as a string
        returns: the reverse complementary DNA sequence represented as a string
    >>> get_reverse_complement("ATGCCCGCTTT")
    'AAAGCGGGCAT'
    >>> get_reverse_complement("CCGCGTTCA")
    'TGAACGCGG'
    """
    reverse_dna = []
    for i in range(len(dna)):
        reverse_dna.append(get_complement(dna[i]))
    reverse_dna.reverse()
    return ''.join(reverse_dna)


def rest_of_ORF(dna):
    """ Takes a DNA sequence that is assumed to begin with a start
        codon and returns the sequence up to but not including the
        first in frame stop codon.  If there is no in frame stop codon,
        returns the whole string.

        dna: a DNA sequence
        returns: the open reading frame represented as a string

    >>> rest_of_ORF("ATGTGAA")
    'ATG'
    >>> rest_of_ORF("ATGAGATAGG")
    'ATGAGA'

    # if there is no dna input
    >>> rest_of_ORF("")
    ''

    # if there is no frame stop codons
    >>> rest_of_ORF("ATGAAAAAAAAAAAAAA")
    'ATGAAAAAAAAAAAAAA'

    """
    end_codon = ["TAG", "TAA", "TGA"]
    framed_dna = []
    result = []
    begin = 0
    end = 3
    # reorganizes dna into a list with interval of 3
    while (True):
        if(begin >= len(dna)):
            break
        framed_dna.append(dna[begin:end])
        begin += 3
        end += 3
    # creates rest_of_ORF
    for codon in framed_dna:
        if(codon not in end_codon):
            result.append(codon)
        else:
            break
    return ''.join(result)


def find_all_ORFs_oneframe(dna):
    """ Finds all non-nested open reading frames in the given DNA
        sequence and returns them as a list.  This function should
        only find ORFs that are in the default frame of the sequence
        (i.e. they start on indices that are multiples of 3).
        By non-nested we mean that if an ORF occurs entirely within
        another ORF, it should not be included in the returned list of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_oneframe("ATGCATGAATGTAGATAGATGTGCCC")
    ['ATGCATGAATGTAGA', 'ATGTGCCC']

    # dna that doesn't start with the start codon
    >>> find_all_ORFs_oneframe("TATGCATGAATG")
    ['ATG']

    # dna that doesn't have any ORFs
    >>> find_all_ORFs_oneframe("")
    []

    # dna that doesn't have any ORFs (version 2)
    >>> find_all_ORFs_oneframe("TATTATGGGGAAAA")
    []

    """
    start_codon = "ATG"
    framed_dna = []
    result = []
    begin = 0  # default frame index
    end = 3  # default frame index

    # divides the dna into frames of 3
    while (True):
        if(begin >= len(dna)):
            break
        framed_dna.append(dna[begin:end])
        begin += 3
        end += 3

    k = 0  # codon index
    result = []  # initialization of result list
    while (True):
        if(k > len(framed_dna) - 1):
            break
        if(framed_dna[k] != start_codon):
            k += 1
        elif(framed_dna[k] == start_codon):
            result.append(rest_of_ORF(''.join(framed_dna[k:])))
            k = k + len(rest_of_ORF(''.join(framed_dna[k:])))/3
    return result


def find_all_ORFs(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence in
        all 3 possible frames and returns them as a list.  By non-nested we
        mean that if an ORF occurs entirely within another ORF and they are
        both in the same frame, it should not be included in the returned list
        of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs

    >>> find_all_ORFs("ATGCATGAATGTAG")
    ['ATGCATGAATGTAG', 'ATGAATGTAG', 'ATG']

    # no dna input should return empty list
    >>> find_all_ORFs("")
    []

    # custom test
    >>> find_all_ORFs("ATGAATGAAATG")
    ['ATGAATGAAATG', 'ATGAAATG']

    """
    result = []
    for i in range(3):  # frame changes from 0~2
        result.extend(find_all_ORFs_oneframe(dna[i:]))
    return result


def find_all_ORFs_both_strands(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence on both
        strands.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_both_strands("ATGCGAATGTAGCATCAAA")
    ['ATGCGAATG', 'ATGCTACATTCGCAT']
    """
    result = []
    result.extend(find_all_ORFs(dna))
    result.extend(find_all_ORFs(get_reverse_complement(dna)))
    return result


def longest_ORF(dna):
    """ Finds the longest ORF on both strands of the specified DNA and returns it
        as a string
    >>> longest_ORF("ATGCGAATGTAGCATCAAA")
    'ATGCTACATTCGCAT'
    """
    # TODO: implement this
    pass


def longest_ORF_noncoding(dna, num_trials):
    """ Computes the maximum length of the longest ORF over num_trials shuffles
        of the specfied DNA sequence

        dna: a DNA sequence
        num_trials: the number of random shuffles
        returns: the maximum length longest ORF """
    # TODO: implement this
    pass


def coding_strand_to_AA(dna):
    """ Computes the Protein encoded by a sequence of DNA.  This function
        does not check for start and stop codons (it assumes that the input
        DNA sequence represents an protein coding region).

        dna: a DNA sequence represented as a string
        returns: a string containing the sequence of amino acids encoded by the
                 the input DNA fragment

        >>> coding_strand_to_AA("ATGCGA")
        'MR'
        >>> coding_strand_to_AA("ATGCCCGCTTT")
        'MPA'
    """
    # TODO: implement this
    pass


def gene_finder(dna):
    """ Returns the amino acid sequences that are likely coded by the specified dna

        dna: a DNA sequence
        returns: a list of all amino acid sequences coded by the sequence dna.
    """
    # TODO: implement this
    pass

if __name__ == "__main__":
    import doctest
    #doctest.testmod()
    doctest.run_docstring_examples(get_complement, globals())
    doctest.run_docstring_examples(get_reverse_complement, globals())
    doctest.run_docstring_examples(rest_of_ORF, globals())
    doctest.run_docstring_examples(find_all_ORFs_oneframe, globals())
    doctest.run_docstring_examples(find_all_ORFs, globals())
    doctest.run_docstring_examples(find_all_ORFs_both_strands, globals())
