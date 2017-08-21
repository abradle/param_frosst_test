import json

def return_atom_types():
  electroneg      = '[N,n,O,o,F,S,s,Cl,Br]'
  ring_electroneg = '[N,n,O,o,s,S]' # are the non-arom needed?
  branched_electroneg = '(~' + electroneg + ')'
  atom_types = [

        # Hydrogen
        # initally from http://www.chem.cmu.edu/courses/09-560/docs/msi/ffbsim/B_AtomTypes.html
        #('H2', '[H][C^3;H1][O]', 0), # O is an example electroneg, more specific than H1
        ('HP',  '[H][C^3;X4][+]',    0), # H on a CT next to atom with formal positive charge
        ('H3a', '[H][C^3;H1]'+branched_electroneg+branched_electroneg+'~'+electroneg, 0),
        ('H2a', '[H][C^3;H1]'+branched_electroneg+electroneg, 0),
        ('H2a', '[H][C^3;H2]'+branched_electroneg+electroneg, 0),
        # ('H5',  '[H][c^2;H1]([n])[n]', 0),
        ('H5',  '[H][c^2;H1]'+branched_electroneg+'~'+electroneg, 0),
        ('H5',  '[H][C^2;H1]'+branched_electroneg+'~'+electroneg, 0),
        ('H1a', '[H][C^3;H1]'+electroneg, 0),
        ('H1b', '[H][C^3][N,n]',  0), # {aliph H; 1 elneg grp}
        ('H1e', '[H][C^3]' + electroneg,  0), # {aliph H; 1 elneg grp}
        ('H4a', '[H][C^2,c^2]'+electroneg, 0), # H4 is H-C&sp2~elneg
        ('H4b', '[H][C^2,c^2]~'+electroneg, 0),
        ('H4c', '[H][C^2,c^2;H1][n]', 0), # was '[H][c^2;H1][n][H]'
        ('H4d', '[H][c^2;H1][s]', 0), # not electroneg! parm_Frosst.pcp inconsistent!
        ('HA', '[H][C^2]',       0),
        ('HA', '[H][cr;^2]',     0), # aromatic H, HA also bonds to [C^2]
        ('H',  '[H][NX3;H1;^2]', 0), # both this  and the one below? Checkme
        ('H',  '[H][nX3;H1;^2]', 0),
        ('HC', '[H][CH3;^3]',    0),
        ('HC', '[H][c,CH1]',     0),
        ('HC', '[H][c,CH2]',     0),
        ('HO', '[H][O;H1]',      0),
        ('HS', '[H][S]',         0),
        ('HW', '[H][O;H2]',      0),
        ('H2b', '[H][C^3]'+branched_electroneg+electroneg, 0), # not sure this works
        # fallback
        ('H', '[H]',    0),

        # -------------------------------pathological specialities ------------------------------
	#
	# fused 5,6 unsat systems
	# the v and ws go together (the same SMART) for normal 5,6 ring and "his-like" sets:
	# normal 5,6 ring
        ('CBv', 'c12aaaac1ccn2',  (0,5)), # indole, hit the C[5,6] carbons
        ('CCv', 'c12aaaac1ccn2',  6),     # indole, hits a non-fused carbon
        ('CWv', 'c12aaaac1ccn2',  7),     # indole, hits a non-fused carbon

        ('CBw', 'c12aaaac1aaa2',  (0,5)), # hit the C[5,6] carbons # which way round do we go?
        ('CCw', 'c12aaaac1caa2',  6),     # hits a non-fused carbon
        ('CWw', 'c12aaaac1aca2',  7),     # hits a non-fused carbon
        # 5,6 ring with "HIS-like" 5-ring
	('CBx', '[cr5]12aaaa[c5]2-'+electroneg+'=[cr5]-'+ring_electroneg+'1', (0,5)),
	('CRx', '[cr5]12aaaa[c5]2-'+electroneg+'=[cr5]-'+ring_electroneg+'1', 7),

	# real histidine-like
        ('CWy', '[cr5]1[cr5]'+ring_electroneg+'[cr5]'+ring_electroneg+'1', 0), # HIS
        ('CCy', '[cr5]1[cr5]'+ring_electroneg+'[cr5]'+ring_electroneg+'1', 1), # HIS
        ('CRy', '[cr5]1[cr5]'+ring_electroneg+'[cr5]'+ring_electroneg+'1', 3), # HIS # success - all hits!
	#
	# maybe we should allow CWy to match CC and CCy to match CW - 2 way rounds an imidazole

        ('CWb',  '[cr5;R1]1[cr5][cr5][ar5][nr5;H1]1', 0), # correct - all hits
        ('CCa',  '[cr5]1[cr5][cr5][ar5][ar5]1', 1),
        ('CWd',  '[cr5;H0]1@[cr5;H1]@[sr5]@[ar5]@[ar5]1', 1), # heuristic C5 in PF5
        ('CWe',  '[cr5;H0]1@[cr5;H0]([I,F,Cl,Br])@[sr5]@[ar5]@[ar5]1', 1), # heuristic C6 in PF5
        ('CC!',  '[cr5;H0]1[sr5][ar5][ar5][cr5]1', 0),

        # 5-ring unsat C&ar5=N&ar5-g&ar5~g&ar5~g&ar5-@1	-> CR NB * * *
        ('CR-5ring-a',  '[C,c;^2]1~[N]~[r5;^3]~[r5]~[C,c]1', 0),
        ('NB-5ring-a',  '[C,c;^2]1~[N]~[r5;^3]~[r5]~[C,c]1', 1),

	# 5-ring unsat C&ar5=C&ar5-C&ar5=C&ar5-g&ar5-@1 -> CW CC CC CW
        ('CW-5ring-b',  '[cr5]1[cr5][cr5][cr5][ar5]1', (0,3)),
        ('CC-5ring-b',  '[cr5]1[cr5][cr5][cr5][ar5]1', (1,2)),

	# 5-ring unsat
        ('NBz',  '[nr5;R1]1[cr5;R1][cr5R1][nr5R1][ar5R1]1', (0,3)), # should this be here?
        ('CRz',  '[nr5;R1]1[cr5;R1][cr5R1][nr5R1][ar5R1]1', (1,2)),


	# ------------------------------ back to sanity ---------------------------

        # Carbon

	# carbon order: CP {rings/fused-rings} C CR CB C* CA CM C2 CJ CT
	#

        ('CPa', '[cr6]-[cr6]', 0), # biphenyl bridge
        ('CPb', '[cr6]-[CR6]', 0), # biphenyl bridge, doesn't hit anything
        ('CPc', '[cr6]=[cr6]', 0), # biphenyl bridge, doesn't hit anything
        ('CPd', '[CR6]-[cr6]', 0), # biphenyl bridge, doesn't hit anything
        ('CPe', 'c1c(cccccc)aacc1', 1), # biphenyl bridge, hits C7 in PF23 (CP[abcd] does not).

        # ('CCb',  '[R]1:[c]:[R]:[R]:[R]:1', 1), # hits things that should be CR or CW
        ('CCc',   '[cr5][c;H0;X3;r5][cr5][nr5][cr5]', 1),

        # ('CCd', '[cr5;X3;^2]', 0), # 5 ring unsaturated

        ('CRa',   '[cr5]1[nr5;H0][ar5][ar5][ar5]1', 0), # this is useful
        ('CRb',   'n[cr5;R2]n', 1),
        ('CRc',  'n[cr5;R1;H0]n', 1),
        ('CRd', 'n1[cH1]ncc1', 1), # NE2 HIS

        # ('Ca',  '[c]:[cX3]=[N]', 1), # was '[c]~[cX3]=[N]' for C10 in PF7, but hits CRs
        ('Ca',  '[C]-[CX3]=[N]', 1),
        ('Cb',  '[CX3]=[O]',  0),
        ('Cc',  '[cX3]=[O]',  0),
        ('Cd',  '[H][CR0]=[N]', 1),
        ('Ce',  '[H][CR0]=[N]', 1), # also conyl
        ('Cf',  '[!C][CR0]=[N]', 1), # also conyl

	# CB is fused aromatic bridge-head
        ('CBa', 'c12aaaac1aaa2', 0), # doesn't hit anything
        ('CBb', 'c12aaaac1aaaa2', (0,5)),
        ('CBc', 'c12aaaac1aan2',  (0,5)), # indole

        # atom CG of a TYR has type CA
        #
        ('CAa', '[a][cr6;H0;R1]([C,c,N,n])[a]', 1),
        ('CAb', 'c', 0),  # CA on fusion atoms in PF6
        ('CAc', '[cr6;X3]', 0),
        ('CAd', '[cr6;H1]', 0),
        ('CAe', '[cr6]',    0), # Does this work? (PF6)
        ('CAf', '[Cr6]',    0), #  [CR6] means in 6 rings :-)
        ('CAg', 'c1ccccc1',  (0,1,2,3,4,5)),
        ('CAh', 'c1ccccn1',  (0,1,2,3,4)),
        ('CAg', '[C^2]1~[C^2]~[C^2]~[C^2]~[C^2]~[C^2]1',  (0,1,2,3,4,5)),
        ('CAh', '[C^2]1~[C^2]~[C^2]~[C^2]~[C^2]~N1',  (0,1,2,3,4)),
        #('CAi', '[C^2]1~[N^2]~[C^2]~[A]~[N^2]~N1',  0,), No hits
        #('CAj', '[C^2]1~[A^2]~[A^2]~[A]~[A]~N1',  0,), # going backwards
        #('CAk', '[C^2]1~[A^2]~[A^2]~[A]~[A]~A1',  0,), # going backwards

	# SMARTS put an atom in a 5 ring before a 6 ring. but CA is 6-ring
        ('C*', '[cr5;^2;X3]', 0),

        ('CMa', 'n1cnccc1', (3,4,5)), # positions 5 or 6 in pyrimidine # or 4 presumably!
        ('CMb', '[C^2]=[C^2]', (0,1)),  # by validation
        ('CMc', '[C]=A', 0),   # by validation
        ('CMd', '[C]=a', 0),   # by validation

        ('C2', '[CX2]#[CX2]', (0,1)),   # by validation
        ('C2', '[CX2]#A', 0),   # by validation
        ('C2', 'A#[CX2]', 1),   # by validation [1]

        ('CWd', '[cr5;X3;^2][c]=O', 0),
        ('CWa', 'n[cr5;R2][cr5;R2]c', 1), # C4 in PF4

        # ('CK', '[cr5;H1;N2]', 0), CK is not a thing in parm@Frosst

        # ('CJ', 'n1n[cX3]cc1', (2,3,4)), # positions 5 or 6 in pyrimidine # or 4 presumably! - not Amber

	# cyclopropyl and epoxide
	('CJa',   '[CX4]1[CX4][CX4]1',  (0,1,2)),
	('CJb',   '[CX4]1[CX4]A1', (0,1)), # remove unindex symmetry

        # ('CN', '[C]'), # how is this different to CB? CN is not a thing in parm@Frosst

        # ('CQ', 'n1[cH1]nccc1', 1),  # CQ is not a thing in parm@Frosst
        # ('CV', '[cr5;^2;H1]~n', 0), # CV is not a thing in parm@Frosst

        ('CT', '[CX4]', 0), # bonded to 4 things
        # Carbon fallback
        ('C',  '[C,c]', 0), # sp hybrizided. Hmmm.

        # [1] needed because above [CX2]#A doesn't find the bond both ways.  I suppose
        #     that's the way SMARTS work.


        # Oxygen
        ('OS',  "[OX2;H0]", 0), # ester
        ('OS',  "[oX2;H0]", 0), # aromatic, should I add [n]?
        ('OH',  "[OH1]", 0), # alcohol
        ('OW',  "[OH2]", 0), # water
        # O2 should match sulphate (and phosphates, I think) but not
        # O=S-{non-oxygen}
        ('O2',  'OP(~O)~O',   0), # O on a phosphate
        ('O2',  'OS(~O)~O',   0), # guess based on above
        ('O2',  'O=P',   0), # O on a phosphate
        ('O2',  'O=S=O',   (0,2)), # guess based on above
        ('O',   'O=*',   0), # carbonyl oxygen
        # fallback
        ('Ou',  'O',   0),
        ('Ou',  'o',   0),


        # Nitrogen
	#
	# Order: NJ NL N3 NC N2 NB N N2 N NA N3 N2 N3 ND N2 N N* N3 Nu
        #
        # Amber pathological cases first:
        # C&ar5=N&ar5-g&ar5~g&ar5~g&ar5-@1        > CR NB * * * ; { 5-ring unsat }
        # N&ar5=C&ar5-C&ar5=N&ar5~g&ar5-@1        > NB CR CR NB * ; { 5-ring unsat }
        ('NB', '[CR5]1[NR5][AR5][AR5][AR5]1', 1),
        ('NB', '[NR5]1[AR5][AR5][NR5][AR5]1', (0,3)),

        ('N2',   '[NX3;!R]=[R]',     0), # {sp2 amino N} - sigh PF7 - N1 is not N!
        ('N',    '[NX3;H1;^2]C=O',   0), # amide
        ('N',    '[NX3;H0;^2]C=O',   0), # PF5
        ('NA',   '[nr5;H1]',      0), # both this and the one below? Checkme
        ('NC',   '[nr6;X2;H0]',   0), # pyridine
        ('NB',   '[c][nr;X2;H0]',  0), # {e.g. HIS C=N-C} # does this catch anything? checkme
        ('NB',   '[c,s][nr5;R1;H0]c',  1),
        ('NB',   '[c;H1,s][nr5;R1;H0]n',  1), # hits N4, N5 in PF4.
        ('NB',   '[cr5;R2][nr5;R1;H0]n',  1), # hits N4, N5 in PF4.
        ('N*',   '[nr5;X2;H0]',   0),
        ('N*',   '[nr6;X2;H0]',   0),
        ('N*',   '[nr5;R2;X3;H0]', 0), # N3 in PF4
        ('NT',   '[N^3;X4]',      0),
        ('N3',   '[N^3;X3]',      0),
        ('N2',   '[NX3;H2^2]', 0),     # N of sp2 NH2 (as in ARG).
        ('N2',   'aN', 1),     # book
        ('N2',   '[NX3][C][NX3]', (0,2)),     # book
        ('N2',   '[H][NX3]CO', 1),     # book
        ('N2',   '[H][NX3]C=O', 1),     # book
        ('Nc',    '[NX3;H1;^2]',   0), # amide
        ('NC',   '[NR6;X2;H0]',   0),
        ('NL',   '[NX1]#A',   0), # triple bond N
        # fall-back nitrogen
        ('Nu',    '[N,n]',      0),

        ('SO', '[S](=O)[N]', 0),  # hypervalent sulfur
        ('S',  '[S,s][S,s]', (0,1)), # sulfide
        ('S',  '[s]', 0), # "sulfide" - hmm. PF5
        # sulfur
        ('Su', 'S', 0),
        ('Su', 's', 0), # yikes

        # P
        ('P',    'P', 0),
        # Cl
        ('CL',   '[Cl]', 0),
        # F
        ('F',    '[F]',  0),
        # Br
        ('BR',   '[Br]',  0),

        ]
  return atom_types