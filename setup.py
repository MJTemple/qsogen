from setuptools import setup

setup(
    name='qsogen',
    version='1.1.0',    
    description='Model Quasar SEDs',
    url='https://github.com/MJTemple/qsogen',
    author='Matthew Temple',
    author_email='Matthew.Temple@mail.udp.cl',
    license='MIT',
    packages=['qsogen'],
    install_requires=['numpy',
                      'scipy',
                      'astropy',
                      ],
)
