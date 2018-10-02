from distutils.core import setup

setup(name='Contouring',
      version='0.1.0',
      author='Pavlos Papaconstadopoulos',
      author_email='pavpapac@gmail.com',
      url='https://github.com/pavpapac/Contouring',
      license='MIT',
      packages=['cmp_contours', ],
      description='A study of tumour contour structures outlined on CT images',
      long_description=open('README.md').read(),
      install_requires=[
          'numpy', 'shapely', 'pydicom'
      ],
      entry_points={
          'console_scripts': [
              'contouring = cmp_contours.study:main',
          ]
      })
