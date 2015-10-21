# helicase
A Python3 script to simulate genetic transcription, translation, mutation, and repair on lists of DNA bases.

## Example:

Assuming a file example.dna with the contents:

```
CATGTAACC
catgccccccccctaatct
```

helicase could be used thus:

```
print("Loading from File:")
strands = helicase.load_strands_from_file("example.dna")
print(strands)

framed_strands = []
for strand in strands:
    framed_strands.append(helicase.frame_strand(strand))
print(framed_strands)

polypeptides = []
for strand in strands:
    polypeptides.append(helicase.represent_polypeptide(helicase.translate_unframed_strand(strand),1))
print(polypeptides)
```

to produce:

```
['catgtaacc', 'catgccccccccctaatct']
[['atg', 'taa', 'cc'], ['atg', 'ccc', 'ccc', 'ccc', 'taa', 'tct']]
['Met', 'Met/Pro/Pro/Pro']
```

A better example can be found at example.py.
