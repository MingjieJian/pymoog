import setuptools

setuptools.setup(
      name='pymoog',
      version='0.0.1',
      description='The python wrapper to run LTE spectra synthesis code MOOG.',
      long_description="",
      long_description_content_type="text/markdown",
      url='https://github.com/MingjieJian/ir_ldr',
      author='Mingjie Jian, Pranav Satheesh and Kruthi Krishna',
      author_email='ssaajianmingjie@gmail.com',
      classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Framework :: IPython",
        "Operating System :: OS Independent",
        "Development Status :: 3 - Alpha",
        "Topic :: Scientific/Engineering :: Astronomy"
      ],
      packages=["pymoog"],
      install_requires=[
          'numpy >= 1.18.0',
          'pandas >= 1.0.0',
          'PyAstronomy >= 0.15.0',
          'matplotlib >= 3.1.0',
          'mendeleev >= 0.6.0'
      ],
      include_package_data=True,
      zip_safe=False)