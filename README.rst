FGH
===

About
-----

This code is a Fortran/Python implementation of Marston's "Fourier Grid
Hamiltonian" as described in D. Tannor's book "Introduction to Quantum
Mechanics: a time-dependent perspective". This method allows the computation
of quantum eigenstates of bound 1d systems.

Usage
-----

The Hamiltonian matrix to be diagonalized is computed by a Fortran routine
(double loop!). The routine can be compiled to a Python extension using
``f2py`` (command given in the source). The potential definition, sampling,
diagonalisation and plotting is performed in the Python script ``FGH.py``.

License
-------

Copyright (c) 2012 Alexander Eberspächer

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
