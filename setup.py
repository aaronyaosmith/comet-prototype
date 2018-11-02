from setuptools import setup
from docs.source import conf

setup(
    name='hgmd',
    version=conf.version,
    packages=['hgmd'],
    entry_points={
        'console_scripts': [
            'hgmd = hgmd.__main__:main'
        ]
    }
)
