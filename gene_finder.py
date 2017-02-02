# -*- coding: utf-8 -*-
"""
The gene_finder function recieves a DNA sequence as input,
then outputs snippets of DNA that are likely to be protein-coding genes

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
        character must me uppercase
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

    # code from Oliver (using dictionary)
    # COMPLEMENTS = {
    #     'C': 'G',
    #     'G': 'C',
    #     'A': 'T',
    #     'T': 'A'
    # }
    # return COMPLEMENTS[nucleotide]

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
    # reorganizes dna string into a list of codons
    while (begin < len(dna)):
        framed_dna.append(dna[begin:end])
        begin += 3
        end += 3
    # creates rest_of_ORF
    for codon in framed_dna:
        if(codon in end_codon):
            break
        else:
            result.append(codon)
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

    # reorganizes dna string into a list of codons
    while (begin < len(dna)):
        framed_dna.append(dna[begin:end])
        begin += 3
        end += 3

    k = 0  # codon index
    result = []  # initialization of result list
    while (k < len(framed_dna)):
        if(framed_dna[k] != start_codon):
            k += 1
        elif(framed_dna[k] == start_codon):
            new_ORF = rest_of_ORF(''.join(framed_dna[k:]))
            result.append(new_ORF)
            k = k + len(new_ORF)/3
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

    >>> find_all_ORFs_both_strands("ATGCATGAATGTAG")
    ['ATGCATGAATGTAG', 'ATGAATGTAG', 'ATG', 'ATGCAT']
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

    >>> longest_ORF("ATGCATGAATGTAG")
    'ATGCATGAATGTAG'
    """
    ORFs = find_all_ORFs_both_strands(dna)
    longest = ''
    for ORF in ORFs:
        if(len(ORF) > len(longest)):
            longest = ORF
    return longest


def longest_ORF_noncoding(dna, num_trials):
    """ Computes the maximum length of the longest ORF over num_trials shuffles
        of the specfied DNA sequence

        dna: a DNA sequence
        num_trials: the number of random shuffles
        returns: the maximum length longest ORF """

    max_length = 0
    for i in range(num_trials):
        shuffled_dna = shuffle_string(dna)  # shuffles the dna randomly
        longest = longest_ORF(shuffled_dna)
        if(len(longest) > max_length):
                max_length = len(longest)
    return max_length


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
    framed_dna = []
    begin = 0  # default frame index
    end = 3  # default frame index

    # reorganizes dna string into a list of codons
    while (begin < len(dna)):
        framed_dna.append(dna[begin:end])
        begin += 3
        end += 3

    result = []
    for codon in framed_dna:
        if(len(codon) == 3):  # there might exist incomplete codon
            result.append(aa_table[codon])
    return ''.join(result)


def gene_finder(dna):
    """ Returns the amino acid sequences that are likely coded by the specified dna

        dna: a DNA sequence (assumed to be all lowercase)
        returns: a list of all amino acid sequences coded by the sequence dna.
    """
    threshold = longest_ORF_noncoding(dna, 1500)  # computes threshold
    ORFs = find_all_ORFs_both_strands(dna)
    result = []
    for ORF in ORFs:
        if(len(ORF) >= threshold):  # only encodes ORFs longer than threshold
            result.append(coding_strand_to_AA(ORF))
    return result


if __name__ == "__main__":
    import doctest
    doctest.testmod()
    dna = load_seq("./data/X73525.fa")
    result = gene_finder(dna)
    for gene in result:
        print(gene)
