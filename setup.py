from setuptools import setup

setup(
    name='qsogen',
    version='1.1.0',    
    description='Model Quasar SEDs',
    url='https://github.com/MJTemple/qsogen',
    author='Matthew Temple',
    author_email='Matthew.Temple@mail.udp.cl',
    license='MIT',
    package_dir={'qsogen': 'qsogen'},
    package_data={'qsogen': ['S0_template_norm.sed',
                             'pl_ext_comp_03.sph',
                             'qsosed_emlines_20210625.dat',
                             'vega_2007.lis',
                             'filters/*.filter',]
                   },
    packages=['qsogen'],
    install_requires=['numpy',
                      'scipy',
                      'astropy',
                      ],
    include_package_data=True,
)
