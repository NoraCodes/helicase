import helicase
import logging

logging.basicConfig(level=logging.WARNING)

print("Loading from File:")
strands = helicase.load_bases_from_file("example.dna")
print(strands)

print("\nFraming:")
framed_strands = []
for strand in strands:
    framed_strands.append(helicase.frame_strand(strand))
print(framed_strands)

print("\nTranscription:")
transcribed_strands = []
for strand in strands:
    transcribed_strands.append(helicase.transcribe(strand))
print(transcribed_strands)

print("\nRNA Transcription:")
transcribed_rna_strands = []
for strand in strands:
    transcribed_rna_strands.append(helicase.transcribe_to_rna(strand))
print(transcribed_rna_strands)


print("\nTranslation:")
polypeptides = []
for strand in strands:
    polypeptides.append(helicase.represent_polypeptide(helicase.translate_unframed_strand(strand),helicase.NORMAL))
print(polypeptides)

print("\nVerbose Translation:")
verbose_polypeptides = []
for strand in strands:
    verbose_polypeptides.append(helicase.represent_polypeptide(helicase.translate_unframed_strand(strand),helicase.VERBOSE))
print(verbose_polypeptides)

print("\nTerse Translation:")
terse_polypeptides = []
for strand in strands:
    terse_polypeptides.append(helicase.represent_polypeptide(helicase.translate_unframed_strand(strand),helicase.SINGLE_CHAR))
print(terse_polypeptides)

print("\n\nStrand -> Transcribed")
for i in range(0, len(strands)):
    print(str(strands[i]) + " -> " + str(transcribed_strands[i]))

print("\n\nStrand -> Transcribed RNA")
for i in range(0, len(strands)):
    print(str(strands[i]) + " -> " + str(transcribed_rna_strands[i]))


print("\nStrand -> Polypeptide")
for i in range(0, len(strands)):
    print(str(strands[i]) + " -> " + str(polypeptides[i]))

print("\nStrand -> Verbose Polypeptide")
for i in range(0, len(strands)):
    print(str(strands[i]) + " -> " + str(verbose_polypeptides[i]))

print("\nStrand -> Terse Polypeptide")
for i in range(0, len(strands)):
    print(str(strands[i]) + " -> " + str(terse_polypeptides[i]))
