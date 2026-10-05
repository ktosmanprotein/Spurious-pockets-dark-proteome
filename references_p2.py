# -*- coding: utf-8 -*-
"""Reference list for Paper 2, in Royal Society (Vancouver) style.

Every entry here was retrieved from PubMed and its DOI, journal, volume and
page numbers copied from the record rather than recalled. Nothing in this file
is written from memory: a fabricated or mis-transcribed citation is the single
worst failure available in a manuscript, and it is not recoverable after
publication.

PMIDs are kept alongside each entry so any of them can be re-checked in one
step. They are not printed in the manuscript.

Royal Society style: authors as surname then initials, up to 10 before et al.,
year, title, abbreviated journal in italics, volume in bold, pages, DOI in
parentheses.
"""

# key -> (authors, year, title, journal, volume, pages, doi, pmid)
REFS = [
    ('gardner2002',
     'Gardner MJ, Hall N, Fung E, White O, Berriman M, Hyman RW, Carlton JM, '
     'Pain A, Nelson KE, Bowman S et al.',
     2002, 'Genome sequence of the human malaria parasite Plasmodium falciparum',
     'Nature', '419', '498-511', '10.1038/nature01097', '12368864'),

    # The origin of the term this paper uses throughout. Cited at first use of
    # "dark proteome" rather than left as received vocabulary.
    ('perdigao2015',
     'Perdigão N, Heinrich J, Stolte C, Sabir KS, Buckley MJ, Tabor B, '
     'Signal B, Gloss BS, Hammang CJ, Rost B et al.',
     2015, 'Unexpected features of the dark proteome',
     'Proc. Natl Acad. Sci. USA', '112', '15898-15903',
     '10.1073/pnas.1508380112', '26578815'),

    # Misannotation by homology transfer, quantified. The upstream version of
    # this paper's "putative" result: function prediction overpredicts.
    ('schnoes2009',
     'Schnoes AM, Brown SD, Dodevski I, Babbitt PC',
     2009, 'Annotation error in public databases: misannotation of molecular '
     'function in enzyme superfamilies', 'PLoS Comput. Biol.', '5', 'e1000605',
     '10.1371/journal.pcbi.1000605', '20011109'),

    # The empirical compact-globule relation Rg = 2.2 N^0.38 used in section 4.4.
    # Verified via Crossref; this journal is not indexed in PubMed, so there is
    # no PMID.
    ('hong2009',
     'Hong L, Lei J',
     2009, 'Scaling law for the radius of gyration of proteins and its '
     'dependence on hydrophobicity',
     'J. Polym. Sci. B Polym. Phys.', '47', '207-214', '10.1002/polb.21634',
     None),

    ('li2003',
     'Li L, Stoeckert CJ, Roos DS',
     2003, 'OrthoMCL: identification of ortholog groups for eukaryotic genomes',
     'Genome Res.', '13', '2178-2189', '10.1101/gr.1224503', '12952885'),

    ('zhang2005',
     'Zhang Y, Skolnick J',
     2005, 'TM-align: a protein structure alignment algorithm based on the '
     'TM-score', 'Nucleic Acids Res.', '33', '2302-2309', '10.1093/nar/gki524',
     '15849316'),

    ('depledge2007',
     'Depledge DP, Lower RPJ, Smith DF',
     2007, 'RepSeq - a database of amino acid repeats present in lower '
     'eukaryotic pathogens', 'BMC Bioinformatics', '8', '122',
     '10.1186/1471-2105-8-122', '17428323'),

    ('leguilloux2009',
     'Le Guilloux V, Schmidtke P, Tuffery P',
     2009, 'Fpocket: an open source platform for ligand pocket detection',
     'BMC Bioinformatics', '10', '168', '10.1186/1471-2105-10-168', '19486540'),

    ('zhang2018',
     'Zhang M, Wang C, Otto TD, Oberstaller J, Liao X, Adapa SR, Udenze K, '
     'Bronner IF, Casandra D, Mayho M et al.',
     2018, 'Uncovering the essential genes of the human malaria parasite '
     'Plasmodium falciparum by saturation mutagenesis', 'Science', '360',
     'eaap7847', '10.1126/science.aap7847', '29724925'),

    ('wang2020',
     'Wang Y, Yang HJ, Harrison PM',
     2020, 'The relationship between protein domains and homopeptides in the '
     'Plasmodium falciparum proteome', 'PeerJ', '8', 'e9940',
     '10.7717/peerj.9940', '33062426'),

    ('muhleip2021',
     'Mühleip A, Kock Flygaard R, Ovciarikova J, Lacombe A, Fernandes P, '
     'Sheiner L, Amunts A',
     2021, 'ATP synthase hexamer assemblies shape cristae of Toxoplasma '
     'mitochondria', 'Nat. Commun.', '12', '120',
     '10.1038/s41467-020-20381-z', '33402698'),

    ('jumper2021',
     'Jumper J, Evans R, Pritzel A, Green T, Figurnov M, Ronneberger O, '
     'Tunyasuvunakool K, Bates R, Žídek A, Potapenko A et al.',
     2021, 'Highly accurate protein structure prediction with AlphaFold',
     'Nature', '596', '583-589', '10.1038/s41586-021-03819-2', '34265844'),

    ('amos2022',
     'Amos B, Aurrecoechea C, Barba M, Barreto A, Basenko EY, Bażant W, '
     'Belnap R, Blevins AS, Böhme U, Brestelli J et al.',
     2022, 'VEuPathDB: the eukaryotic pathogen, vector and host bioinformatics '
     'resource center', 'Nucleic Acids Res.', '50', 'D898-D911',
     '10.1093/nar/gkab929', '34718728'),

    ('zhang2022',
     'Zhang C, Shine M, Pyle AM, Zhang Y',
     2022, 'US-align: universal structure alignments of proteins, nucleic '
     'acids, and macromolecular complexes', 'Nat. Methods', '19', '1109-1115',
     '10.1038/s41592-022-01585-1', '36038728'),

    ('vankempen2023',
     'van Kempen M, Kim SS, Tumescheit C, Mirdita M, Lee J, Gilchrist CLM, '
     'Söding J, Steinegger M',
     2023, 'Fast and accurate protein structure search with Foldseek',
     'Nat. Biotechnol.', '42', '243-246', '10.1038/s41587-023-01773-0',
     '37156916'),

    ('varadi2024',
     'Varadi M, Bertoni D, Magana P, Paramval U, Pidruchna I, Radhakrishnan M, '
     'Tsenkov M, Nair S, Mirdita M, Yeo J et al.',
     2024, 'AlphaFold Protein Structure Database in 2024: providing structure '
     'coverage for over 214 million protein sequences',
     'Nucleic Acids Res.', '52', 'D368-D375', '10.1093/nar/gkad1011',
     '37933859'),

    # The genome-scale novel-fold survey this paper is a correction to: the
    # established way of asking the AlphaFold database for novelty.
    ('durairaj2023',
     'Durairaj J, Waterhouse AM, Mets T, Brodiazhenko T, Abdullah M, '
     'Studer G, Tauriello G, Akdel M, Andreeva A, Bateman A et al.',
     2023, 'Uncovering new families and folds in the natural protein universe',
     'Nature', '622', '646-653', '10.1038/s41586-023-06622-3', '37704037'),

    # Direct predecessor for the reassignment half: the same task, same
    # organism, structure search against the PDB.
    ('behrens2024',
     'Behrens HM, Spielmann T',
     2024, 'Identification of domains in Plasmodium falciparum proteins of '
     'unknown function using DALI search on AlphaFold predictions',
     'Sci. Rep.', '14', '10527', '10.1038/s41598-024-60058-x', '38719885'),

    # Direct predecessor for the druggability half. Ligand evidence is
    # transferred from homologues, so the method is blind to the dark set by
    # construction - which is the complementarity argument in section 3.
    ('godinezmacias2025',
     'Godinez-Macias KP, Chen D, Wallis JL, Siegel MG, Adam A, Bopp S, '
     'Carolino K, Coulson LB, Durst G, Thathy V et al.',
     2025, 'Revisiting the Plasmodium falciparum druggable genome using '
     'predicted structures and data mining',
     'npj Drug Discov.', '2', '3', '10.1038/s44386-025-00006-5', '40066064'),

    # Pocket geometry on predicted models is accurate; docking into it is not.
    # Supports the limitation stated in section 3.
    ('karelina2023',
     'Karelina M, Noh JJ, Dror RO',
     2023, 'How accurately can one predict drug binding modes using AlphaFold '
     'models?', 'eLife', '12', 'e89386', '10.7554/eLife.89386', '38131311'),
]

# The companion methods paper. Not on PubMed: it is a preprint deposited on
# Zenodo and under review. Filled in by hand and flagged so that it is checked
# before submission rather than assumed.
COMPANION = ('osman2026',
             'Osman KO', 2026,
             'Three correctable biases inflate estimates of structural '
             'novelty from predicted-structure databases',
             'Zenodo', '', '', '10.5281/zenodo.23046522', None)

# Order of first citation in the text, which is what Vancouver numbering
# requires. Verified against the built manuscript by check_citation_order()
# below, which fails the build if the text and this list disagree.
def check_companion_title():
    """The companion entry is the one reference with no external record to
    verify it against, so verify it against the manuscript itself. An earlier
    build carried an invented title here."""
    import os, glob
    here = os.path.dirname(os.path.abspath(__file__))
    for d in (os.path.join(here, '..', 'paper1_methods'),):
        for f in glob.glob(os.path.join(d, '*.docx')):
            if 'Paper1' not in os.path.basename(f):
                continue
            try:
                import subprocess
                txt = subprocess.run(['pandoc', '-t', 'plain', f],
                                     capture_output=True, text=True).stdout
            except Exception:
                return None
            head = ' '.join(txt.strip().split('\n')[:2]).strip()
            want = COMPANION[3]
            return head.lower().startswith(want.split(' inflate')[0].lower()[:40])
    return None


ORDER = ['gardner2002', 'perdigao2015',                      # 1-2   intro, para 1
         'jumper2021', 'varadi2024', 'durairaj2023',         # 3-5   intro, para 2
         'behrens2024', 'godinezmacias2025',                 # 6-7   prior art, Pf
         'schnoes2009', 'osman2026',                         # 8-9   intro, para 3
         'zhang2005', 'vankempen2023',                       # 10-11 2.1
         'zhang2022', 'muhleip2021', 'zhang2018',            # 12-14 2.2-2.4
         'leguilloux2009',                                   # 15    2.7
         'depledge2007', 'wang2020',                         # 16-17 2.8
         'karelina2023',                                     # 18    3, limitations
         'hong2009',                                         # 19    4.4
         'li2003', 'amos2022']                               # 20-21 4.5
assert sorted(ORDER) == sorted([k for k, *_ in REFS] + [COMPANION[0]]), (
    'ORDER must list every reference exactly once')
BY_KEY = {r[0]: r for r in REFS}
BY_KEY[COMPANION[0]] = COMPANION


def check_citation_order(text):
    """Every reference must first appear in ascending numerical order."""
    import re
    seen, first = set(), []
    for m in re.finditer(r'\[(\d+(?:,\d+)*)\]', text):
        for n in sorted(int(x) for x in m.group(1).split(',')):
            if n not in seen:
                seen.add(n)
                first.append(n)
    bad = [(a, b) for a, b in zip(first, first[1:]) if b < a]
    missing = sorted(set(range(1, len(ORDER) + 1)) - seen)
    return first, bad, missing


def number(key):
    """1-based reference number for in-text citation."""
    return ORDER.index(key) + 1


def formatted(key):
    k, authors, year, title, journal, vol, pages, doi, _pmid = BY_KEY[key]
    # a title ending in '?' or '!' already carries its stop
    s = '%s%s %d %s%s' % (authors, '' if authors.endswith('.') else '.', year,
                          title, '' if title[-1] in '?!.' else '.')
    if journal == 'Zenodo':
        return s + ' Zenodo preprint, not peer reviewed. (doi:%s)' % doi
    bits = ' %s %s' % (journal, vol)
    if pages:
        bits += ', %s' % pages
    return s + bits + '. (doi:%s)' % doi


def all_formatted():
    return [(i + 1, formatted(k)) for i, k in enumerate(ORDER)]


if __name__ == '__main__':
    for n, t in all_formatted():
        print('%2d. %s' % (n, t))
    npm = sum(1 for r in REFS if r[8])
    print('\n%d references; %d verified against PubMed, %d verified via '
          'Crossref, %d filled by hand'
          % (len(ORDER), npm, len(REFS) - npm, 1))
