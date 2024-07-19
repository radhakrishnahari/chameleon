"""
 @file local_subs.py

 Python Chameleon specific substitution rules for the Precision Generator script.

 @copyright 2019-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
                      Univ. Bordeaux. All rights reserved.

 @version 1.3.0
 @author Mathieu Faverge
 @author Florent Pruvost
 @author Nathalie Furmento
 @author Alycia Lisito
 @date 2024-11-14

"""
_extra_blas = [
    # ----- Additional BLAS
    ('',                     'dsgesv',               'dsgesv',               'zcgesv',               'zcgesv'              ),
    ('',                     'sgesplit',             'dgesplit',             'cgesplit',             'zgesplit'            ),
    ('',                     'slascal',              'dlascal',              'clascal',              'zlascal'             ),
    ('',                     'slapack',              'dlapack',              'clapack',              'zlapack'             ),
    ('',                     'stile',                'dtile',                'ctile',                'ztile'               ),
    ('',                     'sgecon',               'dgecon',               'cgecon',               'zgecon'              ),
    ('',                     'spocon',               'dpocon',               'cpocon',               'zpocon'              ),
    ('',                     'strasm',               'dtrasm',               'ctrasm',               'ztrasm'              ),
    ('',                     'sgecfi',               'dgecfi',               'cgecfi',               'zgecfi'              ),
    ('',                     'splgtr',               'dplgtr',               'cplgtr',               'zplgtr'              ),
    ('',                     'splrnk',               'dplrnk',               'cplrnk',               'zplrnk'              ),
    ('',                     'splssq',               'dplssq',               'cplssq',               'zplssq'              ),
    ('',                     'sy2sb',                'sy2sb' ,               'he2hb',                'he2hb'               ),
    ('',                     'she2ge',               'dhe2ge',               'che2ge',               'zhe2ge'              ),
    ('',                     'slatro',               'dlatro',               'clatro',               'zlatro'              ), #=> Replace by getmo/gecmo as in essl
    ('',                     'sbuild',               'dbuild',               'cbuild',               'zbuild'              ), #=> Replace by map function
    ('',                     'sgram',                'dgram',                'cgram',                'zgram'               ),
    ('',                     'slaran',               'dlaran',               'claran',               'zlaran'              ),
    ('',                     'slaran',               'dlaran',               'slaran',               'dlaran'              ),
    ('',                     'slatm1',               'dlatm1',               'clatm1',               'zlatm1'              ),
    ('',                     'slatm1',               'dlatm1',               'slatm1',               'dlatm1'              ),
    ('',                     'sgenm2',               'dgenm2',               'cgenm2',               'zgenm2'              ),
    ('',                     'slag2c_fake',          'dlag2z_fake',          'slag2c',               'dlag2z'              ),
    ('',                     'slag2d',               'slag2d',               'clag2z',               'clag2z'              ),
    ('',                     'slag2h',               'dlag2h',               'slag2h',               'dlag2h'              ),
    ('',                     'hlag2s',               'hlag2d',               'hlag2s',               'hlag2d'              ),
    ('',                     'slag2h',               'dlag2h',               'clag2x',               'zlag2x'              ),
    ('',                     'hlag2s',               'hlag2d',               'xlag2c',               'xlag2z'              ),
    ('',                     'sgepdf',               'dgepdf',               'cgepdf',               'zgepdf'              ),
    ('',                     'scesca',               'dcesca',               'ccesca',               'zcesca'              ),
    ('',                     'sgesum',               'dgesum',               'cgesum',               'zgesum'              ),
    ('',                     'sgersum',              'dgersum',              'cgersum',              'zgersum'             ),
    ('',                     'sprint',               'dprint',               'cprint',               'zprint'              ),
    ('',                     'sgered',               'dgered',               'cgered',               'zgered'              ),
    ('',                     'sgerst',               'dgerst',               'cgerst',               'zgerst'              ),
    ('',                     'sipiv_allreduce',      'dipiv_allreduce',      'cipiv_allreduce',      'zipiv_allreduce'     ),
    ('',                     'sperm_allreduce',      'dperm_allreduce',      'cperm_allreduce',      'zperm_allreduce'     ),
]

_extra_BLAS = [ [ x.upper() for x in row ] for row in _extra_blas ]

subs = {
    # ------------------------------------------------------------
    # replacements applied to mixed precision files.
    'mixed' : [
        ('dstile',                'zctile'                ),
        # double/single,          double/single-complex
        #'12345678901234567890', '12345678901234567890')
        (r'\bdouble',           r'\bCHAMELEON_Complex64_t'),
        (r'\bChamRealDouble',   r'\bChamComplexDouble'    ),
        (r'\bfloat',            r'\bCHAMELEON_Complex32_t'),
        (r'\bChamRealFloat',    r'\bChamComplexFloat'     ),
        (r'\breal\b',           r'\bcomplex\b'            ),

        ('dsgels',               'zcgels'                 ),
        ('dsorgesv',             'zcungesv'               ),
        ('codelet_ds',           'codelet_zc'             ),
    ],
    # ------------------------------------------------------------
    # replacements applied to mixed precision files.
    'normal': [
        # pattern                single                  double                  single-complex          double-complex
        #'12345678901234567890', '12345678901234567890', '12345678901234567890', '12345678901234567890', '12345678901234567890')


        ('int',                  'float',                'double',               'rpk_complex32_t',     r'\brpk_complex64_t'   ),
        ('int',                  'float',                'double',               'CHAMELEON_Complex32_t', r'\bCHAMELEON_Complex64_t'),
        ('ChamPattern',          'ChamRealFloat',        'ChamRealDouble',       'ChamComplexFloat',    r'\bChamComplexDouble' ),
        ('ChamPattern',          'ChamRealFloat',        'ChamRealDouble',       'ChamRealFloat',       r'\bChamRealDouble'    ),
        ('ChamPattern',          'ChamRealFloat',        'ChamRealFloat',        'ChamComplexFloat',    r'\bChamComplexFloat'  ),
        ('ChamPattern',          'ChamRealFloat',        'ChamRealFloat',        'ChamRealFloat',       r'\bChamRealFloat'     ),
        ('ChamPattern',          'ChamRealHalf',         'ChamRealHalf',         'ChamComplexHalf',     r'\bChamComplexHalf'   ),
        ('int',                  'float',                'double',               'complex32',             'complex64'          ),
        ('Int',                  'Float',                'Double',               'Complex32',             'Complex64'          ),
        ('Int',                  'HMAT_SIMPLE_PRECISION','HMAT_DOUBLE_PRECISION','HMAT_SIMPLE_COMPLEX',   'HMAT_DOUBLE_COMPLEX'),

        # ----- Additional BLAS
        ('',                     'sTile',                'dTile',                'cTile',                'zTile'               ),
        ('',                     'sLapack',              'dLapack',              'cLapack',              'zLapack'             ),
        ('',                     'ORMQR',                'ORMQR',                'UNMQR',                'UNMQR'               ),
        ('',                     'ORMLQ',                'ORMLQ',                'UNMLQ',                'UNMLQ'               ),
        ('',                     'ORGQR',                'ORGQR',                'UNGQR',                'UNGQR'               ),
        ('',                     'ORGLQ',                'ORGLQ',                'UNGLQ',                'UNGLQ'               ),
        ('',                     'SYEV',                 'SYEV',                 'HEEV',                 'HEEV'                ),
        ('',                     'SYG',                  'SYG',                  'HEG',                  'HEG'                 ),
        ('',                     'const float ',         'const double ',       r'const void \*',       r'const void \*'       ), # Line that is extremely dangerous.
    ]
    + _extra_blas
    + _extra_BLAS
    + [

        # ----- For norms: compute result in Real or Double
        ('',                     'slange',               'dlange',               'slange',               'dlange'              ),
        ('',                     'slaset',               'dlaset',               'slaset',               'dlaset'              ),
        ('',                     'splssq',               'dplssq',               'splssq',               'dplssq'              ),
        ('',                     'slacpy',               'dlacpy',               'slacpy',               'dlacpy'              ),
        ('',                     'saxpy',                'daxpy',                'saxpy',                'daxpy'               ),

        (r'\b',                 r'szero\b',             r'dzero\b',             r'czero\b',             r'zzero\b'             ),
#        (r'\b',                 r'sone\b',              r'done\b',              r'cone\b',              r'zone\b'              ),

        # ----- Chameleon Prefixes
        ('CHAMELEON_P',          'CHAMELEON_S',          'CHAMELEON_D',          'CHAMELEON_C',          'CHAMELEON_Z'         ),
        ('CHAMELEON_p',          'CHAMELEON_s',          'CHAMELEON_d',          'CHAMELEON_c',          'CHAMELEON_z'         ),
        ('RUNTIME_P',            'RUNTIME_s',            'RUNTIME_d',            'RUNTIME_c',            'RUNTIME_z'           ),
        ('chameleon_p',          'chameleon_s',          'chameleon_d',          'chameleon_c',          'chameleon_z'         ),
        ('codelet_p',            'codelet_ds',           'codelet_ds',           'codelet_zc',          r'codelet_zc\b'        ),
        ('codelet_p',            'codelet_s',            'codelet_d',            'codelet_c',            'codelet_z'           ),
        ('runtime_p',            'runtime_s',            'runtime_d',            'runtime_c',            'runtime_z'           ),
        ('testing_p',            'testing_s',            'testing_d',            'testing_c',            'testing_z'           ),
        ('workspace_p',          'workspace_s',          'workspace_d',          'workspace_c',          'workspace_z'         ),
        ('check_p',              'check_s',              'check_d',              'check_c',              'check_z'             ),
#        ('CORE_P',               'CORE_S',               'CORE_D',               'CORE_C',               'CORE_Z'              ),
#        ('vec_p',                'vec_s',                'vec_d',                'vec_c',                'vec_z'               ),

        # ----- Vendors Prefixes
        ('vendor_p',             'vendor_s',             'vendor_d',             'vendor_c',             'vendor_z'            ),
        ('',                     'cublasS',              'cublasD',              'cublasC',              'cublasZ'             ),
      # ('',                     'starpu_s',             'starpu_d',             'starpu_c',             'starpu_z'            ),
      # ('',                     'STARPU_S',             'STARPU_D',             'STARPU_C',             'STARPU_Z'            ),
      # ('',                     's_',                   'd_',                   'c_',                   'z_'                  ),
      # ('',                     'S_',                   'D_',                   'C_',                   'Z_'                  ),
      # ('',                     'FLT_EPSILON',          'DBL_EPSILON',          'FLT_EPSILON',          'DBL_EPSILON'         ),
      # ('',                     's_RAFF_FLOAT',         'd_RAFF_FLOAT',         'c_RAFF_FLOAT',         'z_RAFF_FLOAT'        ),
      # # ----- unused?
      # ('',                     's_check',              'd_check',              'c_check',              'z_check'             ),
      # ('',                     'stesting',             'dtesting',             'ctesting',             'ztesting'            ),
      # ('',                     'SAUXILIARY',           'DAUXILIARY',           'CAUXILIARY',           'ZAUXILIARY'          ),
      # ('',                     'sbuild',               'dbuild',               'cbuild',               'zbuild'              ),
        # Hmat-OSS kernels
        ('hmat_p',               'hmat_s',               'hmat_d',               'hmat_c',               'hmat_z'              ),
        ('HMatrix<P_t>',         'HMatrix<S_t>',         'HMatrix<D_t>',         'HMatrix<C_t>',         'HMatrix<Z_t>'        ),
        ('RkMatrix<P_t>',        'RkMatrix<S_t>',        'RkMatrix<D_t>',        'RkMatrix<C_t>',        'RkMatrix<Z_t>'       ),
        ('FullMatrix<P_t>',      'FullMatrix<S_t>',      'FullMatrix<D_t>',      'FullMatrix<C_t>',      'FullMatrix<Z_t>'     ),
        ('ScalarArray<P_t>',     'ScalarArray<S_t>',     'ScalarArray<D_t>',     'ScalarArray<C_t>',     'ScalarArray<Z_t>'    ),
        # ----- RAPACK Prefixes
        ('RapackPattern',        'RapackFloat',          'RapackDouble',         'RapackComplex32',     r'\bRapackComplex64'   ),
        ('RapackPattern',        'RapackFloat',          'RapackDouble',         'RapackFloat',         r'\bRapackDouble'      ),
        ('rpk_p',                'rpk_s',                'rpk_d',                'rpk_c',                'rpk_z'               ),
        ('rpkx_p',               'rpkx_s',               'rpkx_d',               'rpkx_c',               'rpkx_z'              ),
        ('rpk_precgentmp_diag_t','rpk_precgentmp_diag_t','rpk_precgentmp_diag_t','rpk_precgentmp_diag_t','rpk_diag_t'          ),
        ('rpk_p',                'rpk_s',                'rpk_d',                'rpk_s',                'rpk_d'               ),
        ('rpkx_p',               'rpkx_s',               'rpkx_d',               'rpkx_s',               'rpkx_d'              ),
        ('RPK_P',                'RPK_S',                'RPK_D',                'RPK_C',                'RPK_Z'               ),
        ('RPK_P',                'RPK_S',                'RPK_D',                'RPK_S',                'RPK_D'               ),
        ('prapack',              'srapack',              'drapack',              'crapack',              'zrapack'             ),
        ('rpk_',                 'rpk_',                 'rpk_',                 'rpk_',                 'rpk_precgentmp_'     ),
        ('pge2lr',               'sge2lr',               'dge2lr',               'cge2lr',               'zge2lr'              ),
        ('CHAMELEON_plr',               'sge2lr',               'dge2lr',               'cge2lr',               'zge2lr'              ),
        ('pge2ps',               'sge2ps',               'dge2ps',               'cge2ps',               'zge2ps'              ),
        ('phmat2rpk',            'shmat2rpk',            'dhmat2rpk',            'chmat2rpk',            'zhmat2rpk'           ),
        ('ptile2rpk',            'stile2rpk',            'dtile2rpk',            'ctile2rpk',            'ztile2rpk'           ),
        ('prpk2fmt',             'srpk2fmt',             'drpk2fmt',             'rpk2fmt',              'zrpk2fmt'            ),
        ('plrcpr',             'slrcpr',             'dlrcpr',             'lrcpr',              'zlrcpr'            ),
        ('pfmt2fmt',             'sfmt2fmt',             'dfmt2fmt',             'fmt2fmt',              'zfmt2fmt'            ),
        ('plrcrk',             'slrcrk',             'dlrcrk',             'lrcrk',              'zlrcrk'            ),
    ]
}

exceptfrom = [ r'\brpk_diag_t' ]
