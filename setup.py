from setuptools import setup

setup(
    name = 'hgmd',
    version = '0.1.0',
    packages = ['hgmd'],
    entry_points = {
        'console_scripts': [
            'hgmd = hgmd.__main__:main'
        ]
    }
)
