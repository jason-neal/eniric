try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

config = {
    'description': 'Extended NIR Information Content',
    'author': 'Jason Neal',
    'url': 'https://github.com/jason-neal/eniric.git',
    'download_url': 'https://github.com/jason-neal/eniric.git',
    'author_email': 'jason.neal@astro.up.pt',
    'version': '0.1',
    'install_requires': ['pytest'],
    'packages': ['eniric'],
    'scripts': [],
    'name': 'eniric'
}

setup(**config)
