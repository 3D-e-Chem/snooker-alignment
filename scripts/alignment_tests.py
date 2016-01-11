"""
Tests for alignment script

Usage:

nosetests scripts/

Or with coverage:

nosetests --with-coverage --cover-branches --cover-html \
--cover-package=alignment scripts/

"""
import cPickle
from unittest import TestCase

import alignment


class fix_wrong_lengthsTestCase(TestCase):
    def setUp(self):
        self.score_matrix = cPickle.load(open('scripts/fixtures/score_matrix.pickle'))

    def seq_oklength_nochanges_test(self):
        seq = 'MENETVSELNQTQLQPRAVVALEYQVVTILLVLIICGLGIVGNIMVVLVVMRTKHMRTPT'
        seq += 'NCYLVSLAVADLMVLVAAGLPNITDSIYGSWVYGYVGCLCITYLQYLGINASSCSITAFT'
        seq += 'IERYIAICHPIKAQFLCTFSRAKKIIIFVWAFTSLYCMLWFFLLDLNISTYKDAIVISCG'
        seq += 'YKISRNYYSPIYLMDFGVFYVVPMILATVLYGFIARILFLNPIPSDPKENSKTWKNDSTH'
        seq += 'QNTNLNVNTSNRCFNSTVSSRKQVTKMLAVVVILFALLWMPYRTLVVVNSFLSSPFQENW'
        seq += 'FLLFCRICIYLNSAINPVIYNLMSQKFRAAFRKLCNCKQKPTEKPANYSVALNYSVIKES'
        seq += 'DHFSTELDDITVTDTYLSATKVSFDDTCLASEVSFSQS'
        tmseq = 'GCLCITYLQYLGINASSCSITAFTIERYIA'
        aln = {
            'TRFR_HUMAN': {
                3: {
                    'helixseq': tmseq,
                    'sequence': seq,
                    'firstaa': 96,
                    'lastaa': 126
                }
            }
        }
        tms = {3: 30}
        preferred_gaps = {3: []}
        alignment.fix_wrong_lengths(aln, tms, self.score_matrix, preferred_gaps)

        expected_seq = tmseq
        self.assertEquals(aln['TRFR_HUMAN'][3]['helixseq'], expected_seq)
        self.assertEquals(aln['TRFR_HUMAN'][3]['firstaa'], 96)
        self.assertEquals(aln['TRFR_HUMAN'][3]['lastaa'], 126)

    def seq_toosmall_alnextended_test(self):
        seq = 'MENETVSELNQTQLQPRAVVALEYQVVTILLVLIICGLGIVGNIMVVLVVMRTKHMRTPT'
        seq += 'NCYLVSLAVADLMVLVAAGLPNITDSIYGSWVYGYVGCLCITYLQYLGINASSCSITAFT'
        seq += 'IERYIAICHPIKAQFLCTFSRAKKIIIFVWAFTSLYCMLWFFLLDLNISTYKDAIVISCG'
        seq += 'YKISRNYYSPIYLMDFGVFYVVPMILATVLYGFIARILFLNPIPSDPKENSKTWKNDSTH'
        seq += 'QNTNLNVNTSNRCFNSTVSSRKQVTKMLAVVVILFALLWMPYRTLVVVNSFLSSPFQENW'
        seq += 'FLLFCRICIYLNSAINPVIYNLMSQKFRAAFRKLCNCKQKPTEKPANYSVALNYSVIKES'
        seq += 'DHFSTELDDITVTDTYLSATKVSFDDTCLASEVSFSQS'
        tmseq = 'LCITYLQYLGINASSCSITAFTIE'
        aln = {
            'TRFR_HUMAN': {
                3: {
                    'helixseq': tmseq,
                    'sequence': seq,
                    'firstaa': 98,
                    'lastaa': 122
                }
            }
        }
        tms = {3: 30}
        preferred_gaps = {3: []}
        alignment.fix_wrong_lengths(aln, tms, self.score_matrix, preferred_gaps)


        expected_seq = 'GCLCITYLQYLGINASSCSITAFTIERYIA'
        self.assertEquals(aln['TRFR_HUMAN'][3]['helixseq'], expected_seq)
        self.assertEquals(aln['TRFR_HUMAN'][3]['firstaa'], 96)
        self.assertEquals(aln['TRFR_HUMAN'][3]['lastaa'], 126)

    def seq_waytoosmall_alnextended_test(self):
        seq = 'MENETVSELNQTQLQPRAVVALEYQVVTILLVLIICGLGIVGNIMVVLVVMRTKHMRTPT'
        seq += 'NCYLVSLAVADLMVLVAAGLPNITDSIYGSWVYGYVGCLCITYLQYLGINASSCSITAFT'
        seq += 'IERYIAICHPIKAQFLCTFSRAKKIIIFVWAFTSLYCMLWFFLLDLNISTYKDAIVISCG'
        seq += 'YKISRNYYSPIYLMDFGVFYVVPMILATVLYGFIARILFLNPIPSDPKENSKTWKNDSTH'
        seq += 'QNTNLNVNTSNRCFNSTVSSRKQVTKMLAVVVILFALLWMPYRTLVVVNSFLSSPFQENW'
        seq += 'FLLFCRICIYLNSAINPVIYNLMSQKFRAAFRKLCNCKQKPTEKPANYSVALNYSVIKES'
        seq += 'DHFSTELDDITVTDTYLSATKVSFDDTCLASEVSFSQS'
        tmseq = 'QYLGINASSCSI'
        aln = {
            'TRFR_HUMAN': {
                3: {
                    'helixseq': tmseq,
                    'sequence': seq,
                    'firstaa': 104,
                    'lastaa': 116
                }
            }
        }
        tms = {3: 30}
        preferred_gaps = {3: []}
        alignment.fix_wrong_lengths(aln, tms, self.score_matrix, preferred_gaps)

        expected_seq = 'GCLCITYLQYLGINASSCSITAFTIERYIA'
        self.assertEquals(aln['TRFR_HUMAN'][3]['helixseq'], expected_seq)
        self.assertEquals(aln['TRFR_HUMAN'][3]['firstaa'], 96)
        self.assertEquals(aln['TRFR_HUMAN'][3]['lastaa'], 126)

    def seq_oklengthwithgaps_nochanges_test(self):
        seq = 'MENETVSELNQTQLQPRAVVALEYQVVTILLVLIICGLGIVGNIMVVLVVMRTKHMRTPT'
        seq += 'NCYLVSLAVADLMVLVAAGLPNITDSIYGSWVYGYVGCLCITYLQYLGINASSCS-TAFT'
        seq += 'IERYIAICHPIKAQFLCTFSRAKKIIIFVWAFTSLYCMLWFFLLDLNISTYKDAIVISCG'
        seq += 'YKISRNYYSPIYLMDFGVFYVVPMILATVLYGFIARILFLNPIPSDPKENSKTWKNDSTH'
        seq += 'QNTNLNVNTSNRCFNSTVSSRKQVTKMLAVVVILFALLWMPYRTLVVVNSFLSSPFQENW'
        seq += 'FLLFCRICIYLNSAINPVIYNLMSQKFRAAFRKLCNCKQKPTEKPANYSVALNYSVIKES'
        seq += 'DHFSTELDDITVTDTYLSATKVSFDDTCLASEVSFSQS'
        tmseq = 'GCLCITYLQYLGINASSCS-TAFTIERYIA'
        aln = {
            'TRFR_HUMAN': {
                3: {
                    'helixseq': tmseq,
                    'sequence': seq,
                    'firstaa': 96,
                    'lastaa': 126
                }
            }
        }
        tms = {3: 30}
        preferred_gaps = {3: []}
        alignment.fix_wrong_lengths(aln, tms, self.score_matrix, preferred_gaps)

        expected_seq = tmseq
        self.assertEquals(aln['TRFR_HUMAN'][3]['helixseq'], expected_seq)
        self.assertEquals(aln['TRFR_HUMAN'][3]['firstaa'], 96)
        self.assertEquals(aln['TRFR_HUMAN'][3]['lastaa'], 126)

    def seq_toosmallwithgapsinside_alnextended_test(self):
        seq = 'MENETVSELNQTQLQPRAVVALEYQVVTILLVLIICGLGIVGNIMVVLVVMRTKHMRTPT'
        seq += 'NCYLVSLAVADLMVLVAAGLPNITDSIYGSWVYGYVGCLCITYLQYLGINASSCS-TAFT'
        seq += 'IERYIAICHPIKAQFLCTFSRAKKIIIFVWAFTSLYCMLWFFLLDLNISTYKDAIVISCG'
        seq += 'YKISRNYYSPIYLMDFGVFYVVPMILATVLYGFIARILFLNPIPSDPKENSKTWKNDSTH'
        seq += 'QNTNLNVNTSNRCFNSTVSSRKQVTKMLAVVVILFALLWMPYRTLVVVNSFLSSPFQENW'
        seq += 'FLLFCRICIYLNSAINPVIYNLMSQKFRAAFRKLCNCKQKPTEKPANYSVALNYSVIKES'
        seq += 'DHFSTELDDITVTDTYLSATKVSFDDTCLASEVSFSQS'
        tmseq = 'LCITYLQYLGINASSCS-TAFTIE'
        aln = {
            'TRFR_HUMAN': {
                3: {
                    'helixseq': tmseq,
                    'sequence': seq,
                    'firstaa': 98,
                    'lastaa': 122
                }
            }
        }
        tms = {3: 30}
        preferred_gaps = {3: []}
        alignment.fix_wrong_lengths(aln, tms, self.score_matrix, preferred_gaps)

        expected_seq = 'GCLCITYLQYLGINASSCS-TAFTIERYIA'
        self.assertEquals(aln['TRFR_HUMAN'][3]['helixseq'], expected_seq)
        self.assertEquals(aln['TRFR_HUMAN'][3]['firstaa'], 96)
        self.assertEquals(aln['TRFR_HUMAN'][3]['lastaa'], 126)

    def add_preferred_gaps_test(self):
        seq = 'MAPWPHENSSLAPWPDLPTLAPNTANTSGLPGVPWEAALAGALLALAVLATVGGNLLVIV'
        seq += 'AIAWTPRLQTMTNVFVTSLAAADLVMGLLVVPPAATLALTGHWPLGATGCELWTSVDVLC'
        seq += 'VTASIETLCALAVDRYLAVTNPLRYGALVTKRCARTAVVLVWVVSAAVSFAPIMSQWWRV'
        seq += 'GADAEAQRCHSNPRCCAFASNMPYVLLSSSVSFYLPLLVMLFVYARVFVVATRQLRLLRG'
        seq += 'ELGRFPPEESPPAPSRSLAPAPVGTCAPPEGVPACGRRPARLLPLREHRALCTLGLIMGT'
        seq += 'FTLCWLPFFLANVLRALGGPSLVPGPAFLALNWLGYANSAFNPLIYCRSPDFRSAFRRLL'
        seq += 'CRCGRRLPPEPCAAARPALFPSGVPAARSSPAQPRLCQRLDGASWGVS'
        tmseq = 'YVLLSSSVSFYLPLLVMLFVYARVFV'
        aln = {
            'ADRB3_HUMAN': {
                5: {
                    'helixseq': tmseq,
                    'sequence': seq,
                    'firstaa': 203,
                    'lastaa': 229
                }
            }
        }
        tms = {5: 31}
        preferred_gaps = {5: [5, 11]}
        alignment.fix_wrong_lengths(aln, tms, self.score_matrix, preferred_gaps)
        expected_seq = 'MPYV-LLSSS-VSFYLPLLVMLFVYARVFVV'
        self.assertEquals(aln['ADRB3_HUMAN'][5]['helixseq'], expected_seq)
        self.assertEquals(aln['ADRB3_HUMAN'][5]['firstaa'], 201)
        self.assertEquals(aln['ADRB3_HUMAN'][5]['lastaa'], 230)


class cut_out_aligned_partTestCase(TestCase):

    def nogaps_test(self):
        hmmaln = '................................................................'
        hmmaln += '....................menetvselnqtqlqpravvaleyqvvtillvliicglgivgn'
        hmmaln += 'imvvlvvmrtkhmrtptncylvslavadlmvlvaaglpnitdsiygswvygyvGCLCITYL.Q'
        hmmaln += '.YLGIN.A..SSCSI.T..AFTIERYIAichpikaqflctfsrakkiiifvwaftslycmlwf'
        hmmaln += 'flldlnistykdaiviscgykisrnyyspiylmdfgvfyvvpmilatvlygfiarilflnpip'
        hmmaln += 'sdpkensktwkndsthqntnlnvntsnrcfnstvssrkqvtkmlavvvilfallwmpyrtlvv'
        hmmaln += 'vnsflsspfqenwfllfcriciylnsainpviynlmsqkfraafrklcnckqkptekpanysv'
        hmmaln += 'alnysvikesdhfsteldditvtdtylsatkvsfddtclasevsfsqs.........'

        aln = {
            'TRFR_HUMAN': {
                3: {
                    'sequence': hmmaln,
                }
            }
        }
        tm_lengths = {3: 30}

        alignment.cut_out_aligned_part(aln, tm_lengths)

        seq = 'menetvselnqtqlqpravvaleyqvvtillvliicglgivgn'
        seq += 'imvvlvvmrtkhmrtptncylvslavadlmvlvaaglpnitdsiygswvygyvGCLCITYLQ'
        seq += 'YLGINASSCSITAFTIERYIAichpikaqflctfsrakkiiifvwaftslycmlwf'
        seq += 'flldlnistykdaiviscgykisrnyyspiylmdfgvfyvvpmilatvlygfiarilflnpip'
        seq += 'sdpkensktwkndsthqntnlnvntsnrcfnstvssrkqvtkmlavvvilfallwmpyrtlvv'
        seq += 'vnsflsspfqenwfllfcriciylnsainpviynlmsqkfraafrklcnckqkptekpanysv'
        seq += 'alnysvikesdhfsteldditvtdtylsatkvsfddtclasevsfsqs'
        tmseq = 'GCLCITYLQYLGINASSCSITAFTIERYIA'
        expected = {
            'TRFR_HUMAN': {
                3: {
                    'helixseq': tmseq,
                    'helixlength': 30,
                    'sequence': seq,
                    'seqlength': 398,
                    'firstaa': 96,
                    'lastaa': 126,
                    'cut_off': True
                }
            }
        }
        self.assertEquals(aln, expected)

    def singlegap_test(self):
        hmmaln = '...mdgsnvtsfvveeptnistgrnasvgnahrqipivhwvimsispvgfvengillwflcfrmrrnpft-VYITHLSIA..DISLLF..CIFILSI--dya'
        hmmaln += 'ldyelssghyytivtlsvtflfgyntglylltaisverclsvlypiwyrchrpkyqsalvcallwalsclvttmeyvmcidreeeshsrndcraviifiailsflvftplmlvsstilvvkirkntwashssklyivimvtiiiflifampmrllyllyyeywstfgnlhhisllfstinssanpfiyffvgsskkkrfk'
        hmmaln += 'eslkvvltrafkdemqprrqkdncntvtvetvv.......................................................................................................................................................................'

        aln = {
            'MAS_HUMAN': {
                2: {
                    'sequence': hmmaln,
                }
            }
        }
        tm_lengths = {2: 26}

        alignment.cut_out_aligned_part(aln, tm_lengths)

        seq = 'mdgsnvtsfvveeptnistgrnasvgnahrqipivhwvimsispvgfvengillwflcfrmrrnpftVYITHLSIADISLLFCIFILSIdya'
        seq += 'ldyelssghyytivtlsvtflfgyntglylltaisverclsvlypiwyrchrpkyqsalvcallwalsclvttmeyvmcidreeeshsrndcraviifiailsflvftplmlvsstilvvkirkntwashssklyivimvtiiiflifampmrllyllyyeywstfgnlhhisllfstinssanpfiyffvgsskkkrfk'
        seq += 'eslkvvltrafkdemqprrqkdncntvtvetvv'
        tmseq = 'VYITHLSIADISLLFCIFILSI'
        # TODO expected seq to contain gaps in TM section + tmseq=TVYITHLSIADISLLF-CIFILSIDY
        expected = {
            'MAS_HUMAN': {
                2: {
                    'helixseq': tmseq,
                    'helixlength': 22,
                    'sequence': seq,
                    'seqlength': 325,
                    'firstaa': 67,
                    'lastaa': 89,
                    'cut_off': True
                }
            }
        }
        self.maxDiff = None
        self.assertEquals(aln, expected)

    def doublegap_test(self):
        hmmaln = '...mravfiqgaeehpaafcyqvngscprtvhtlgiqlviylacaagmliivlgnvfvafavsyfkalhtptnflllslaladmflgllvlplstirsvescwffgdflcrlhtyldtlfcltsifhlcfisidrhcaicdpllypskftvrvalryilagwgvpaaytslflytdvvetrlsqwleempc'
        hmmaln += 'vgscqlllnk-FWG.WLNFP.L-FFVPCLIMISLYVKIFV-vatrqaqqittlskslagaakherkaaktlgiavgiyllcwlpftidtmvdsllhfitpplvfdifiwfayfnsacnpiiyvfsyqwfrkalkltlsqkvfspqtrtvdlyqe...........'

        aln = {
            'TAAR5_HUMAN': {
                5: {
                    'sequence': hmmaln,
                }
            }
        }
        tm_lengths = {5: 26}

        alignment.cut_out_aligned_part(aln, tm_lengths)

        seq = 'mravfiqgaeehpaafcyqvngscprtvhtlgiqlviylacaagmliivlgnvfvafavsyfkalhtptnflllslaladmflgllvlplstirsvescwffgdflcrlhtyldtlfcltsifhlcfisidrhcaicdpllypskftvrvalryilagwgvpaaytslflytdvvetrlsqwleempc'
        seq += 'vgscqlllnkFWGWLNFPLFFVPCLIMISLYVKIFVvatrqaqqittlskslagaakherkaaktlgiavgiyllcwlpftidtmvdsllhfitpplvfdifiwfayfnsacnpiiyvfsyqwfrkalkltlsqkvfspqtrtvdlyqe'
        tmseq = 'FWGWLNFPLFFVPCLIMISLYVKIFV'
        # TODO expected seq to contain gaps in TM section + tmseq=NKFW-GWLNF-PLFFVPCLIMISLYVKIFVV
        expected = {
            'TAAR5_HUMAN': {
                5: {
                    'helixseq': tmseq,
                    'helixlength': 26,
                    'sequence': seq,
                    'seqlength': 337,
                    'firstaa': 198,
                    'lastaa': 224,
                    'cut_off': True
                }
            }
        }
        self.maxDiff = None
        self.assertEquals(aln, expected)

