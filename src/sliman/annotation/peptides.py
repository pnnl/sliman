"""
sliman/annotation/peptides.py

Dylan Ross (dylan.ross@pnnl.gov)

    Utilities for representing and performing calculations with peptides
"""


import re

from mzapy.isotopes import MolecularFormula, monoiso_mass


# map 1-letter AA codes to R-group formulas
_AA_R_FORMULAS = {
    'A': MolecularFormula(C=1, H=3), 
    'C': MolecularFormula(C=1, H=3, S=1),
    'D': MolecularFormula(C=2, H=3, O=2),
    'E': MolecularFormula(C=3, H=5, O=2),
    'F': MolecularFormula(C=7, H=7),
    'G': MolecularFormula(H=1),
    'H': MolecularFormula(C=4, H=5, N=2),
    'I': MolecularFormula(C=4, H=9),
    'K': MolecularFormula(C=4, H=10, N=1),
    'L': MolecularFormula(C=4, H=9),
    'M': MolecularFormula(C=3, H=7, S=1),
    'N': MolecularFormula(C=2, H=4, O=1, N=1),
    'P': MolecularFormula(C=3, H=5),
    'Q': MolecularFormula(C=3, H=6, O=1, N=1),
    'R': MolecularFormula(C=4, H=10, N=3),
    'S': MolecularFormula(C=1, H=3, O=1),
    'T': MolecularFormula(C=2, H=5, O=1),
    'U': MolecularFormula(C=1, H=3, Se=1),
    'V': MolecularFormula(C=3, H=7),
    'W': MolecularFormula(C=9, H=8, N=1),
    'X': MolecularFormula(C=3, H=3, O=1),  # pyroglutamic acid
    'Y': MolecularFormula(C=7, H=7, O=1),
}


# map between 1-letter and 3-letter AA codes
_AA_1LETTER_TO_3LETTER = {
    'A': 'Ala', 'C': 'Cys', 'D': 'Asp', 'E': 'Glu', 'F': 'Phe',
    'G': 'Gly', 'H': 'His', 'I': 'Ile', 'K': 'Lys', 'L': 'Leu', 
    'M': 'Met', 'N': 'Asn', 'P': 'Pro', 'Q': 'Gln', 'R': 'Arg',
    'S': 'Ser', 'T': 'Thr', 'U': 'Sec', 'V': 'Val', 'W': 'Trp',
    'X': 'Pyr', 'Y': 'Tyr',
}
_AA_3LETTER_TO_1LETTER = {v: k for k, v in _AA_1LETTER_TO_3LETTER.items()}


class _RGroup():
    """

    Attributes
    ----------
    aa : ``str``
        1-letter AA code
    formula : ``dict(str:int)``
        molecular formula of R group
    """

    def __init__(self, aa):
        """

        Parameters
        ----------
        aa : ``str``
            1-letter or 3-letter code specifying amino acid
        """
        if aa in _AA_R_FORMULAS:
            self.aa = aa
            self.formula = _AA_R_FORMULAS[self.aa]
        elif aa in _AA_3LETTER_TO_1LETTER:
            self.aa = _AA_3LETTER_TO_1LETTER[aa]
            self.formula = _AA_R_FORMULAS[self.aa]
        else:
            msg = '_RGroup: __init__: amino acid code "{}" not recognized'
            raise ValueError(msg.format(aa))


class _AminoAcid():
    """
    object encapsulating an amino acid residue

                 R
                 |
    <nterm>--NH--CH--CO--<cterm>

    Attributes
    ----------
    n_aa : ``_AminoAcid`` or ``mzapy.isotopes.MolecularFormula``
        reference to previous amino acid (or MolecularFormula(H=1) if this AA is the N-terminus)
    c_aa : ``_AminoAcid`` or ``mzapy.isotopes.MolecularFormula``
        reference to next amino acid (or MolecularFormula(O=1, H=1) if this AA is the C-terminus)
    idx_nc : ``int``
        index of amino acid in peptide (numbered from N-terminus to C-terminus)
    idx_cn : ``int``
        index of amino acid in peptide (numbered from C-terminus to N-terminus)
    mod : ``None``
        chemical modification
        TODO (Dylan Ross): not implemented yet
    name1 : ``str``
        1-letter code name
    name3 : ``str``
        3-letter code name
    """

    def __init__(self, aa, idx_nc, idx_cn, n_aa=None, c_aa=None):
        """
        inits a new instance of _AminoAcid with the AA type and indexing information specified

        Parameters
        ----------
        aa : ``str``
            1-letter or 3-letter code specifying amino acid
        idx_nc : ``int``
            index of amino acid in peptide (numbered from N-terminus to C-terminus)
        idx_cn : ``int``
            index of amino acid in peptide (numbered from C-terminus to N-terminus)
        n_aa : ``_AminoAcid``, optional
            reference to previous amino acid (set to MolecularFormula(H=1) if not provided)
        c_aa : ``_AminoAcid`` or ``mzapy.isotopes.MolecularFormula``
            reference to next amino acid (set to MolecularFormula(O=1, H=1) if not provided)
        """
        self.mod=None
        self._r = _RGroup(aa)
        self.name1, self.name3 = self._r.aa, _AA_1LETTER_TO_3LETTER[self._r.aa]
        self.idx_nc, self.idx_cn = idx_nc, idx_cn  # validate?
        self.formula = MolecularFormula(C=2, H=2, O=1, N=1) + self._r.formula
        if n_aa is not None:
            self.n_aa = n_aa
            # automatically set the c_aa reference on the other _AminoAcid instance to self
            self.n_aa.c_aa = self
        else:
            self.n_aa = MolecularFormula(H=1)
        if c_aa is not None:
            self.c_aa = c_aa
            # automatically set the n_aa reference on the other _AminoAcid instance to self
            self.c_aa.n_aa = self
        else:
            self.c_aa = MolecularFormula(O=1, H=1)

    def __repr__(self):
        return '_AminoAcid("", idx_nc={}, idx_cn={})'.format(self.name1, self.idx_nc, self.idx_cn)

    def __str__(self):
        s = '-NH-({})-CO-'.format(self.name3)
        if isinstance(self.n_aa, MolecularFormula):
            s = 'H' + s
        else:
            s = '<...>' + s
        if isinstance(self.c_aa, MolecularFormula):
            s += 'OH'
        else:
            s += '<...>'
        return s

    def _formula_with_nterm(self):
        """ returns the molecular formula including this AA and anything on the N-terminal side """
        n_formula = self.formula.copy()
        if isinstance(self.n_aa, MolecularFormula):
            return n_formula + self.n_aa
        else:
            return n_formula + self.n_aa._formula_with_nterm()

    def _formula_with_cterm(self):
        """ returns the molecular formula including this AA and anything on the C-terminal side """
        c_formula = self.formula.copy()
        if isinstance(self.c_aa, MolecularFormula):
            return c_formula + self.c_aa
        else:
            return c_formula + self.c_aa._formula_with_cterm()

    def get_fragments(self, fragment_types):
        """
        produces molecular formulas for typical peptide fragments (a, b, c, x, y, z)

        the produced fragments are produced from the neutral peptide formula and have +1 charge
        
        Parameters
        ----------
        fragment_types : ``list(str)``
            list of fragment types to produce formulas for, if not provided defaults to:
                ['a', 'b', 'c', 'x', 'y', 'z']

        Returns
        -------
        fragments : ``dict(str:mzapy.isotopes.MolecularFormula)``
            dict mapping fragment type and index (e.g., "b2" or "y3") to fragment formula
        """
        # map fragment type to the modification for each type of fragment
        # modifications are relative to NEUTRAL peptide formula, but produce +1 charged fragments
        ft = {
            'a': {'C': -1, 'O': -1}, 
            'b': {},  # no change 
            'c': {'N': 1, 'H': 3},
            'x': {'C': 1, 'O': 1}, 
            'y': {'H': 2}, 
            'z': {'N': -1, 'H': -1}
        }
        fragment_types = fragment_types if fragment_types is not None else list('abcxyz')
        fragments = {}
        for fragment_type in fragment_types:
            if fragment_type in ['a', 'b', 'c']:
                if self.idx_cn > 1:  # skip if at the C-terminus
                    # N-terminal fragments
                    lbl = '{}{}'.format(fragment_type, self.idx_nc)
                    fragments[lbl] = self._formula_with_nterm() + ft[fragment_type]
            elif fragment_type in ['x', 'y', 'z']:
                if self.idx_nc > 1:  # skip if at the N-terminus
                    # C-terminal fragments
                    lbl = '{}{}'.format(fragment_type, self.idx_cn)
                    fragments[lbl] = self._formula_with_cterm() + ft[fragment_type]
            else:
                msg = '_AminoAcid: get_fragments: fragment_type "{}" not recognized'
                raise ValueError(msg.format(fragment_type))
        return fragments


class Peptide():
    """
    object encapsulating a peptide

    Attributes
    ----------
    seq1 : ``str``
        amino acid sequence, N->C, 1-letter codes
    seq3 : ``str``
        amino acid sequence, N->C, 3-letter codes
    n : ``int``
        count of amino acids in peptide
    formula : ``mzapy.isotopes.MolecularFormula``
        molecular formula for the peptide
    mass : ``float``
        monoisotopic mass (neutral) for the peptide
    """

    def __init__(self, sequence):
        """
        inits a new Peptide instance from an amino acid sequence (1- or 3-letter AA codes)

        Parameters
        ----------
        sequence : ``str``
            amino acid sequence as 3-letter or 1-letter codes
        """
        self._amino_acids = []
        parsed_seq = self._parse_sequence(sequence)
        self.n = len(parsed_seq)
        for i, aa in enumerate(parsed_seq):
            if i == 0:
                self._amino_acids.append(_AminoAcid(aa, 1, self.n - i))
            else:
                self._amino_acids.append(_AminoAcid(aa, i + 1, self.n - i, n_aa=self._amino_acids[i - 1]))

    def _parse_sequence(self, sequence):
        """ parses the amino acid sequence (``str``) and returns the split sequence (``list(str)``) """
        # determine if the codes are 1-letter or 3-letter
        pat_1letter = re.compile('^[A-Z]+$')
        pat_3letter = re.compile('^([A-Z][a-z]{2})+$')
        pat_3letter_h = re.compile('^([A-Z][a-z]{2}-)+([A-Z][a-z]{2})$')  # hyphen-separated 3-letter codes
        if pat_1letter.match(sequence):
            return list(sequence)
        elif pat_3letter.match(sequence):
            return re.findall('...', sequence)
        elif pat_3letter_h.match(sequence):
            return sequence.split('-')
        else:
            msg = 'Peptide: _parse_sequence: unable to parse peptide sequence "{}"'
            raise ValueError(msg.format(sequence))

    def __len__(self):
        return len(self._amino_acids)

    @property
    def seq1(self):
        return ''.join([_.name1 for _ in self._amino_acids])

    @property
    def seq3(self):
        return '-'.join([_.name3 for _ in self._amino_acids])

    @property
    def formula(self):
        formula = MolecularFormula(H=2, O=1)
        for aa in self._amino_acids:
            formula += aa.formula
        return formula

    @property
    def mass(self):
        return monoiso_mass(self.formula)

    def __repr__(self):
        return 'Peptide("{}")'.format(self.seq1)

    def __str__(self):
        return self.seq1

    def get_fragment_formulas(self, fragment_types=None):
        """
        produces molecular formulas for typical peptide fragments (a, b, c, x, y, z) from 
        all possible cleavages in the peptide. The fragments are produced from the neutral peptide 
        formula and have +1 charge
        
        Parameters
        ----------
        fragment_types : ``list(str)``, optional
            list of fragment types to produce formulas for, if not provided defaults to:
                ['a', 'b', 'c', 'x', 'y', 'z']

        Returns
        -------
        fragment_formulas : ``dict(str:mzapy.isotopes.MolecularFormula)``
            acts like dict mapping fragment type and index (e.g., "b2" or "y3") to fragment formula
        """
        fragment_formulas = {}
        for aa in self._amino_acids:
            fragments = aa.get_fragments(fragment_types)
            for frag, form in fragments.items():
                fragment_formulas[frag] = form.copy()
        return fragment_formulas

    def get_fragment_masses(self, fragment_types=None):
        """
        produces monoisotopic masses for typical peptide fragments (a, b, c, x, y, z) from 
        all possible cleavages in the peptide. The fragments are produced from the neutral peptide 
        formula and have +1 charge
        
        Parameters
        ----------
        fragment_types : ``list(str)``, optional
            list of fragment types to produce masses for, if not provided defaults to:
                ['a', 'b', 'c', 'x', 'y', 'z']

        Returns
        -------
        fragment_masses : ``dict(str:float)``
            dict mapping fragment type and index (e.g., "b2" or "y3") to fragment monoisotopic mass
        """
        # get the formulas then convert to masses and return
        fragment_formulas = self.get_fragment_formulas(fragment_types)
        return {k: monoiso_mass(v) for k, v in fragment_formulas.items()}

