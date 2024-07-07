from setuptools import setup, find_packages

setup(
    name='bolcan',
    version='0.0.1',
    packages=find_packages('tests*', 'docs*'),
    install_requires=[
        'pandas',
    ],
    author='Justice Nyame',
    author_email='nyamejustice2000@gmail.com',
    description='A Python package offering essential tools for petroleum and other engineering calculations',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    url='https://github.com/Jnyame21/bolcan',
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.6',
)


