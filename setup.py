import os
from setuptools import setup, find_packages

def read(fname):
    with open(os.path.join(os.path.dirname(__file__), fname)) as f:
        return f.read()

package_name = 'pipepermcalc'
setup(
    name=package_name,
    version='{{VERSION_PLACEHOLDER}}',
    packages=find_packages(exclude=['tests*', 'tests.*', 'research*',]),
    license='All rights reserved',
    description='A python package to ',
    long_description=read('README.rst'),
    long_description_content_type="text/x-rst",
    
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3.7',
    ],
    python_requires='>=3.7',
    project_urls={
    'Source': f'https://github.com/KWR-Water/{package_name}',
    'Documentation': f'http://{package_name}.readthedocs.io/en/latest/',
    'Tracker': f'https://github.com/KWR-Water/{package_name}/issues',
    'Help': f'https://github.com/KWR-Water/{package_name}/issues',
    },
    install_requires=[
        'pandas',
        'openpyxl',
        'pytest',
        'fuzzywuzzy',
        ],
    url=f'https://github.com/KWR-Water/{package_name}',
    author='KWR Water Research Institute',
    author_email='alex.hockin@kwrwater.nl, bram.hillebrand@kwrwater.nl',
    package_data={'': ['database/*.csv']},
)
