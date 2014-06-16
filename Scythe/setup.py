from setuptools import setup, find_packages
setup(
        name='scythe',
        version='0.1a1',
        author='Janina Mass',
        author_email='janina.mass@hhu.de',
        packages=find_packages(),
        scripts=['scythe/convert/scythe_ensembl2grp.py',
            'scythe/convert/scythe_ensembl_fasta.py',
            'scythe/convert/scythe_ensembl_ortho_mysql.py',
            'scythe/convert/scythe_grp_orthomcl.py',
            'scythe/convert/scythe_grp_proteinortho.py',
            'scythe/convert/scythe_grp_tsv.py',
            'scythe/convert/scythe_loc_ensemblfasta.py',
            'scythe/convert/scythe_loc_gff.py',
            'scythe/convert/scythe_loc_tsv.py',
            'scythe/convert/scythe_mergeSubsets.py',
            'scythe/scythe.py',
            'scythe/scythe-gui.py'],
        license='GPLv3',
        description='Find best matching set of transcripts for one-to-one orthologous genes from two or more species',
        long_description=open('README.txt').read(),
        classifiers=[
            'Topic :: Scientific/Engineering :: Bio-Informatics'
            ],
        )
