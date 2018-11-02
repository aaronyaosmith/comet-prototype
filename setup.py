import re

from setuptools import setup


verstr = "unknown"
try:
    verstrline = open('yourpackage/_version.py', "rt").read()
except EnvironmentError:
    pass  # Okay, there is no version file.
else:
    VSRE = r"^__version__ = ['\"]([^'\"]*)['\"]"
    mo = re.search(VSRE, verstrline, re.M)
    if mo:
        verstr = mo.group(1)
    else:
        raise RuntimeError("unable to find version in yourpackage/_version.py")


setup(
    name='hgmd',
    version=verstr,
    packages=['hgmd'],
    entry_points={
        'console_scripts': [
            'hgmd = hgmd.__main__:main'
        ]
    }
)
