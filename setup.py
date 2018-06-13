from distutils.core import setup

setup(
    name='comb_chamber',
    version='0.0.1',
    package_dir={'comb_chamber': 'core', 'comb_chamber.templates': 'core/templates'},
    packages=['comb_chamber', 'comb_chamber.templates'],
    package_data={'comb_chamber': ['templates/comb_geom.tex']},
    url='',
    license='',
    author='Alexander Zhigalkin',
    author_email='',
    description=''
)
