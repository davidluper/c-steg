c-steg
======
 
jpeg steganography library using libjpeg

this repo is a result of playing around with using libjpeg for steganography.  
it is not meant to be a cryptographic tool and should not be for anything requiring
any level of security as neither cryptography or security are addressed at all in this 
code base.  there is an android make file in the src/jpeg folder that should build it 
via the NDK.

to build from the console:
  1. go into src/jpeg and run configure
  2. in the src/jpeg folder run make
  3. from the src directory run make
  
that will make a binary called c-steg in the src directory.  
usage instructions for the binary can be obtained by running the binary with --help (c-steg --help)
