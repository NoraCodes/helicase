# Tests for helicase
import unittest
import helicase
import tempfile
import os
import logging

logging.basicConfig(level=logging.DEBUG)

class TestHelicase(unittest.TestCase):
    def test_transcription(self):
        coding =   "atcg"
        template = "tagc"
        self.assertEqual(template, helicase.transcribe(coding))

    def test_cannot_transcribe_invalid_bases(self):
        coding = "aaaccctttnggg"
        with self.assertRaises(ValueError):
            helicase.transcribe(coding)

    def test_frame_strand(self):
        self.assertEqual(helicase.frame_strand('catgccccccccctaatct'), ['atg', 'ccc', 'ccc', 'ccc', 'taa', 'tct'])

    def test_translate_framed_strand(self):
        framed_strand = helicase.frame_strand('catgccccccccctaatct')
        self.assertEqual(helicase.translate_framed_strand(framed_strand), [helicase.Met, helicase.Pro, helicase.Pro,
                                                                           helicase.Pro])
    def test_translate_unframed_strand(self):
        self.assertEqual(helicase.translate_unframed_strand('catgccccccccctaatct'), [helicase.Met, helicase.Pro,
                                                                                     helicase.Pro, helicase.Pro])

    def test_load_from_file(self):
        (handle, filename) = tempfile.mkstemp()
        handle = open(filename, "w")
        handle.write("catgtaacc\ncatgccccccccctaatct")
        handle.close()
        self.assertEqual(helicase.load_bases_from_file(filename), ['catgtaacc', 'catgccccccccctaatct'])
        handle.close()
        os.unlink(filename)

    def test_represent_polypeptide_raises_valueerror(self):
        with self.assertRaises(ValueError):
            helicase.represent_polypeptide([helicase.Met, helicase.Pro], 99)

    def test_represent_polypeptide(self):
        polypeptide = [helicase.Met, helicase.Pro, helicase.Cys]
        self.assertEqual(helicase.represent_polypeptide(polypeptide, 0), "MPC")
        self.assertEqual(helicase.represent_polypeptide(polypeptide, 1), "Met/Pro/Cys")
        self.assertEqual(helicase.represent_polypeptide(polypeptide, 2), "Methionine, Proline, Cysteine")

    def test_comprehensive(self):
        logging.debug("\nBeginning Comprehensive Test.-----------")
        (handle, filename) = tempfile.mkstemp()
        handle = open(filename, "w")
        handle.write("CATGTAACC\ncatgccccccccctaatct")
        handle.close()
        strands = helicase.load_bases_from_file(filename)
        self.assertEqual(strands, ['catgtaacc', 'catgccccccccctaatct'])
        framed_strands = []
        for strand in strands:
            framed_strands.append(helicase.frame_strand(strand))
        self.assertEqual(framed_strands, [['atg', 'taa', 'cc'], ['atg', 'ccc', 'ccc', 'ccc', 'taa', 'tct']])
        polypeptides = []
        for strand in framed_strands:
            polypeptides.append(helicase.translate_framed_strand(strand))
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
