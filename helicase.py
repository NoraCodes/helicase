#!/usr/bin/python3
# helicase - Simulate transcription, translation, mutation, and repair on lists of DNA bases

# Imports
import sys
import string

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
Ser = ('S', 'Ser', 'Serine')
Gly = ('G', 'Gly', 'Glycine')

base_names = {'a': 'Adenine', 't': 'Thymine', 'g': 'Guanine', 'c': 'Cytosine'}

bases_to_proteins = {'ttt': Phe, 'ttc': Phe,
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
        lower_line = line.lower()
        if not set(lower_line) <= allowed:
            raise ValueError("File " + filename + " contains invalid base names. Remember: DNA not RNA, so no Uracil.")
        # If we reach this point, the strand is valid
        converted_strands.append(lower_line)
    return converted_strands

def frame_strand(strand):
    """
    Take an unframed strand of DNA and frame it, discarding bases before the START marker
    :param strand: The unframed strand
    :return:
    """

    framed_strand = []
    rolling_frame = ["","",""]
    frame_begin = 0
    for i in range(0, len(strand)):
        rolling_frame[0] = rolling_frame[1]
        rolling_frame[1] = rolling_frame[2]
        rolling_frame[2] = strand[i]
        if rolling_frame[0]+rolling_frame[1]+rolling_frame[2] == start_codon:
            # We have the frame at this point. Prune the beginning of the string:
            # We are at a pos+2 (c ccc atg ccc ccc c) so cut back 2 and then break
            #                            ^
            frame_begin = i - 2
            break

    if frame_begin == 0:
        # In this case, there is no valid frame in the strand. Return an empty list.
        return []

    pruned_strand = strand[frame_begin:]
    framed_strand = [pruned_strand[i:i+3] for i in range(0, len(pruned_strand), 3)] # Make triples
    # NOTE: Perhaps add a notification on a non-
    return framed_strand

def transcribe(strand):
    """
    Transcribe a strand into its complementary bases
    :param strand:
    :return:
    """
    transcribed_strand = strand
    for i in range(0, len(strand)):
        if strand[i] == 'g':
            transcribed_strand[i] = 'c'
        elif strand[i] == 'c':
            transcribed_strand[i] = 'g'
        elif strand[i] =='a':
            transcribed_strand[i] = 't'
        elif strand[i] == 't':
            transcribed_strand[i] = 'a'
        else:
            raise ValueError("Tried to transcribe a strand with invalid bases.")
    return transcribed_strand