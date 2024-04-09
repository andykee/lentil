import re

from setuptools import setup

with open('lentil/__init__.py') as f:
    version = re.search("__version__ = '(.*?)'", f.read()).group(1)

setup(
    name='Lentil',
    version=version,
    install_requires=[
        'numpy>=1.17,<2.0',
        'scipy'
    ]
)
