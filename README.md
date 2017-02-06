# fir - additions to the Scipy signals toolbox.

This is a package to compliment the the SciPy signals toolbox.  Some functions are redundant (the remez function) from the SciPy package.  Indeed, the remez function is a direct copy.  However, some are new.

## Installation

This package does require some external packages, particularly the latest version of boost (v1.63).  At the moment this is installed by source.  Information can be found [here](http://www.boost.org/doc/libs/1_63_0/more/getting_started/unix-variants.html).

To use, you will first need to build the program.  Assuming you are on a Linux machine
```shell
cd fir
mkdir obj # build objects will be placed here
make
make install
```
The `install` step is determined where you place the install path.  The install location is established by 
```shell
INSTALL_PATH = /location/of/python/path
```
## Usage
Each of the functions should be documented, but the usage is similar to the Scipy function calls.  

The documentation can be found by launching python in the command line and then
```python
>>>import fir
>>>help(fir)
```
to see that documentation.

