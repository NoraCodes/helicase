# Tests for helicase
import unittest
import helicase
import tempfile
import os
import logging

logging.basicConfig(level=logging.DEBUG)


class TestHelicase(unittest.TestCase):
    def test_transcription(self):
        print("\nBeginning Transcription Test.")
        coding =   "atcg"
        template = "tagc"
        self.assertEqual(template, helicase.transcribe(coding))

    def test_transcription_to_rna(self):
        print("\nBeginning RNA Transcription Test.")
        coding =   "atcg"
        template = "uagc"
        self.assertEqual(template, helicase.transcribe_to_rna(coding))

    def test_cannot_transcribe_invalid_bases(self):
        print("\nBeginning Invalid Transcription Test.")
        coding = "aaaccctttnggg"
        with self.assertRaises(ValueError):
            helicase.transcribe(coding)

    def test_cannot_transcribe_invalid_bases_to_rna(self):
        print("\nBeginning Invalid RNA Transcription Test.")
        coding = "aaaccctttnggg"
        with self.assertRaises(ValueError):
            helicase.transcribe_to_rna(coding)

    def test_frame_strand(self):
        print("\nBeginning Framing Test.")
        self.assertEqual(helicase.frame_strand('catgccccccccctaatct'), (['atg', 'ccc', 'ccc', 'ccc', 'taa', 'tct'], 1))

    def test_translate_framed_strand(self):
        print("\nBeginning Framed Translation Test.")
        framed_strand = helicase.frame_strand('catgccccccccctaatct')[0]
        self.assertEqual(helicase.translate_framed_strand(framed_strand), [helicase.Met, helicase.Pro, helicase.Pro,
                                                                           helicase.Pro])

    def test_translate_unframed_strand(self):
        print("\nBeginning Unframed Translation Test.")
        self.assertEqual(helicase.translate_unframed_strand('catgccccccccctaatct'), [helicase.Met, helicase.Pro,
                                                                                     helicase.Pro, helicase.Pro])

    def test_load_from_string(self):
        self.assertEqual(helicase.load_strands_from_string("catgtaacc\ncatgccccccccctaatct"), ['catgtaacc', 'catgccccccccctaatct'])
        self.assertEqual(helicase.load_strands_from_string("catgtaacc\tcatgccccccccctaatct", seperator = '\t'), ['catgtaacc', 'catgccccccccctaatct'])

    def test_load_from_file(self):
        print("\nBeginning File Loading Test.")
        (handle, filename) = tempfile.mkstemp()
        handle = open(filename, "w")
        handle.write("catgtaacc\ncatgccccccccctaatct")
        handle.close()
        self.assertEqual(helicase.load_strands_from_file(filename), ['catgtaacc', 'catgccccccccctaatct'])
        handle.close()
        os.unlink(filename)

    def test_represent_polypeptide_raises_valueerror(self):
        print("\nBeginning Invalid Verbosity Level Test.")
        with self.assertRaises(ValueError):
            helicase.represent_polypeptide([helicase.Met, helicase.Pro], 99)

    def test_represent_polypeptide(self):
        print("\nBeginning Representation Test.")
        polypeptide = [helicase.Met, helicase.Pro, helicase.Cys]
        self.assertEqual(helicase.represent_polypeptide(polypeptide, helicase.IUPAC_1), "MPC")
        self.assertEqual(helicase.represent_polypeptide(polypeptide, helicase.IUPAC_3), "Met/Pro/Cys")
        self.assertEqual(helicase.represent_polypeptide(polypeptide, helicase.FULL_NAME), "Methionine, Proline, Cysteine")

    def test_comprehensive(self):
        logging.debug("\nBeginning Comprehensive Test.-----------")
        (handle, filename) = tempfile.mkstemp()
        handle = open(filename, "w")
        handle.write("CATGTAACC\ncatgccccccccctaatct")
        handle.close()
        strands = helicase.load_strands_from_file(filename)
        self.assertEqual(strands, ['catgtaacc', 'catgccccccccctaatct'])
        framed_strands = []
        for strand in strands:
            framed_strands.append(helicase.frame_strand(strand))
        self.assertEqual(framed_strands, [(['atg', 'taa', 'cc'], 1), (['atg', 'ccc', 'ccc', 'ccc', 'taa', 'tct'], 1)])
        polypeptides = []
        for strand in framed_strands:
            polypeptides.append(helicase.translate_framed_strand(strand[0]))
        self.assertEqual(polypeptides, [[('M', 'Met', 'Methionine')], [('M', 'Met', 'Methionine'),
                                                                       ('P', 'Pro', 'Proline'), ('P', 'Pro', 'Proline'),
                                                                       ('P', 'Pro', 'Proline')]])
        polypeptide_representations = []
        for polypeptide in polypeptides:
            polypeptide_representations.append(helicase.represent_polypeptide(polypeptide, 1))
        self.assertEqual(polypeptide_representations, ['Met', 'Met/Pro/Pro/Pro'])
        logging.debug("Done Comprehensive Test.-----------\n")



if __name__ == "__main__":
    unittest.main()
