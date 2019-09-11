from setuptools import setup

setup(name='qtb_plot',
      version='0.0.1',
      description='Institute standard plotting styles',
      url='https://gitlab.com/marvin.vanaalst/qtb_plot',
      author='Marvin van Aalst',
      author_email='marvin.van.aalst@hhu.de',
      license='GPL4',
      packages=['qtb_plot'],
      install_requires=[
          'numpy',
          'matplotlib'
      ],
      zip_safe=False)
