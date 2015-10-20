#!/usr/bin/python3
# helicase - Simulate transcription, translation, mutation, and repair on lists of DNA bases

# Imports
import sys
import string
import logging
import argparse

# Lookup Tables
# Ref   L    3let   Full Name
Phe = ('F', 'Phe', 'Phenylalanine')
Leu = ('L', 'Leu', 'Leucine')
Ile = ('I', 'Ile', 'Isoleucine')
Met = ('M', 'Met', 'Methionine')
Val = ('V', 'Val', 'Valine')
Ser = ('S', 'Ser', 'Serine')
Pro = ('P', 'Pro', 'Proline')
Thr = ('T', 'Thr', 'Threonine')
Ala = ('A', 'Ala', 'Alanine')
Tyr = ('Y', 'Tyr', 'Tyrosine')
His = ('H', 'His', 'Histidine')
Gln = ('Q', 'Gln', 'Glutamine')
Asn = ('N', 'Asn', 'Asparagine')
Lys = ('K', 'Lys', 'Lysine')
Asp = ('D', 'Asp', 'Aspartic acid')
Glu = ('E', 'Glu', 'Glutamic acid')
Cys = ('C', 'Cys', 'Cysteine')
Trp = ('W', 'Trp', 'Tryptophan')
Arg = ('R', 'Arg', 'Arginine')
Gly = ('G', 'Gly', 'Glycine')

base_names = {'a': 'Adenine', 't': 'Thymine', 'g': 'Guanine', 'c': 'Cytosine'}

codons_to_amino_acids = {'ttt': Phe, 'ttc': Phe,
                      'tta': Leu, 'ttg': Leu, 'ctt': Leu, 'ctc': Leu, 'cta': Leu, 'ctg': Leu,
                      'att': Ile, 'atc': Ile, 'ata': Ile,
                      'atg': Met,
                      'gtt': Val, 'gtc': Val, 'gta': Val, 'gtg': Val,
                      'tct': Ser, 'tcc': Ser, 'tca': Ser, 'tcg': Ser,
                      'cct': Pro, 'ccc': Pro, 'cca': Pro, 'ccg': Pro,
                      'act': Thr, 'acc': Thr, 'aca': Thr, 'acg': Thr,
                      'gct': Ala, 'gcc': Ala, 'gca': Ala, 'gcg': Ala,
                      'tat': Tyr, 'tac': Tyr,
                      'cat': His, 'cac': His,
                      'caa': Gln, 'cag': Gln,
                      'aat': Asn, 'aac': Asn,
                      'aaa': Lys, 'aag': Lys,
                      'gat': Asp, 'gac': Asp,
                      'gaa': Glu, 'gag': Glu,
                      'tgt': Cys, 'tgc': Cys,
                      'tgg': Trp,
                      'cgt': Arg, 'cgc': Arg, 'cga': Arg, 'cgg': Arg, 'aga': Arg, 'agg': Arg,
                      'agt': Ser, 'agc': Ser,
                      'ggt': Gly, 'ggc': Gly, 'gga': Gly, 'ggg': Gly}
stop_codons = ['taa', 'tag', 'tga']
start_codon = 'atg'


def load_bases_from_file(filename):
    """
    Load bases from the file specified by filename.
    Produce a list of strings, where each string is a valid
    """
    allowed = {'a', 'c', 't', 'g'}
    converted_strands = []

    try:
        bases_file = open(filename, 'r')
    except OSError:
        sys.exit(1) # Exit with an "other" error.

    for line in bases_file:
        lower_line = line.lower().strip("\n\t ")
        for letter in lower_line:
            if not letter in allowed:
                raise ValueError("File " + filename + " contains invalid base name" + letter + ". Remember: DNA not RNA, so no Uracil.")
        # If we reach this point, the strand is valid
        converted_strands.append(lower_line)
    bases_file.close()
    logging.info("Validated and loaded " + str(len(converted_strands)) + " strands from " + filename + ".")
    return converted_strands


def frame_strand(strand):
    """
    Take an unframed strand of DNA and frame it, discarding bases before the START marker
    :param strand: The unframed strand
    :return:
    """
    logging.info("Framing strand: " + strand)
    framed_strand = []
    rolling_frame = ['','','']
    frame_begin = -1
    for i in range(0, len(strand)):
        rolling_frame[0] = rolling_frame[1]
        rolling_frame[1] = rolling_frame[2]
        rolling_frame[2] = strand[i]
        if rolling_frame[0]+rolling_frame[1]+rolling_frame[2] == start_codon:
            # We have the frame at this point. Prune the beginning of the string:
            # We are at a pos+2 (c ccc atg ccc ccc c) so cut back 2 and then break
            #                            ^
            frame_begin = i - 2
            logging.debug("Found start codon, strand is now framed at " + str(frame_begin) + ".")
            break
    if frame_begin == -1:
        # In this case, there is no valid frame in the strand. Return an empty list.
        logging.info("No valid frame in strand.")
        return []
    pruned_strand = strand[frame_begin:]
    framed_strand = [pruned_strand[i:i+3] for i in range(0, len(pruned_strand), 3)] # Make triples
    logging.debug("Framed strand is: " + str(framed_strand))
    if len(framed_strand[-1]) < 3:
        logging.info("Newly framed sequence terminates with non-codon: " + framed_strand[-1])
    return framed_strand


def transcribe(strand):
    """
    Transcribe a strand into its complementary bases
    :param strand:
    :return:
    """
    transcribed_strand = ''
    for i in range(0, len(strand)):
        if strand[i] == 'g':
            transcribed_strand += 'c'
        elif strand[i] == 'c':
            transcribed_strand += 'g'
        elif strand[i] =='a':
            transcribed_strand += 't'
        elif strand[i] == 't':
            transcribed_strand += 'a'
        else:
            raise ValueError("Tried to transcribe a strand with invalid bases.")
    return transcribed_strand


def transcribe_to_rna(strand):
    """
    Transcribe a strand into its complementary bases of RNA
    :param strand:
    :return:
    """
    transcribed_strand = ''
    for i in range(0, len(strand)):
        if strand[i] == 'g':
            transcribed_strand += 'c'
        elif strand[i] == 'c':
            transcribed_strand += 'g'
        elif strand[i] =='a':
            transcribed_strand += 'u'
        elif strand[i] == 't':
            transcribed_strand += 'a'
        else:
            raise ValueError("Tried to transcribe a strand with invalid bases.")
    return transcribed_strand


def translate_unframed_strand(strand):
    """
    Frame and translate a strand into a polypeptide
    :param strand: The unframed strand to translate
    :return: A list of amino acids
    """
    framed_strand = frame_strand(strand)
    return translate_framed_strand(framed_strand)


def translate_framed_strand(framed_strand):
    """
    Translate a strand into a polypeptide
    :param framed_strand: The framed strand to translate
    :return: A list of amino acids
    """

    polypeptide = []
    # Iterate through all the codons in our strand and look up their resultant amino acid
    logging.info("Translating strand: " + str(framed_strand))
    for codon in framed_strand:
        if codon in stop_codons:
            # We are at the end of the sequence, so cut.
            logging.debug("Stopping at codon " + codon)
            break
        logging.debug("Adding polypeptide " + codons_to_amino_acids[codon][2])
        polypeptide.append(codons_to_amino_acids[codon])
    return polypeptide


def represent_polypeptide(polypeptide, level=0):
    """
    Represent a polypeptide as a string.
    :param polypeptide: A list of peptide 3-tuples
    :param level: Representation verbosity level. 0 = PP, 1 = Pro/Pro, 2 = Proline, Proline
    :return: A string
    """
    output_string = ""
    separator = ""
    separator_backspace = 0
    if level == 0:
        separator = ""
        separator_backspace = 0
    elif level == 1:
        separator = "/"
        separator_backspace = 1
    elif level == 2:
        separator = ", "
        separator_backspace = 2
    else:
        raise ValueError("Representation verbosity level must be one of: 0, 1, 2.")

    for amino_acid in polypeptide:
        output_string += amino_acid[level]
        output_string += separator

    if level > 0:
        output_string = output_string[:-separator_backspace]
    return output_string


if __name__ == "__main__":
    # This code will be skipped if we are a library.
    # Define command line args
    parser = argparse.ArgumentParser(prog="helicase.py", description="Convert sequences of DNA nucleotide bases into "+
                                                                     "polypeptides or mRNA bases.")
    # Determine procedure
    # Get input
    # Return output
    pass

