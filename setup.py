from setuptools import setup

setup(
    name='qsogen',
    version='1.1.0',    
    description='Model Quasar SEDs from Temple, Hewett & Banerji (2021) MNRAS 508, 737',
    url='https://github.com/MJTemple/qsogen',
    author='Matthew Temple',
    author_email='Matthew.Temple@mail.udp.cl',
    license='MIT',
    # package_dir={'qsogen': ['qsosed', 'config', 'model_colours'],
    #              'filters': 'filters'},
    package_data={'qsogen': ['S0_template_norm.sed',
                             'pl_ext_comp_03.sph',
                             'qsosed_emlines_20210625.dat',
                             'vega_2007.lis',],
                  'filters': ['*.filter']
                   },
    packages=['qsogen'],
    install_requires=['numpy',
                      'scipy',
                      'astropy',
                      ],
    include_package_data=True,
)
