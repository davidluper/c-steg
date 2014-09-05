c-steg
======

<h2>jpeg steganography library using libjpeg</h2>

<p>this repo is a result of playing around with using libjpeg for steganography.  
it is not meant to be a cryptographic tool and should not be for anything requiring
any level of security as neither cryptography or security are addressed at all in this
code base.  there is an android make file in the src/jpeg folder that should build it
via the NDK.

<b>to build from the console:</b>
  1. go into src/jpeg and run configure
  2. in the src/jpeg folder run make
  3. from the src directory run make

that will make a binary called c-steg in the src directory.  
usage instructions for the binary can be obtained by running the binary with --help (c-steg --help)

<p>there are two steg modes.  the first one is called lossless and it assumes you get the
jpeg you encode back without it being recompressed.  the second mode is called lossy and
it tries to accomodate cropping or recompression of the jpeg after encoding it with information<br /><br />

<b>to run a lossless encode from the command line run this command</b><br />
c-steg --encode-lossless --input-file	./my-photo.jpg --output-file ./my-photo-encoded.jpg \
--text "this is the text i want to encode in the jpg" --quality	95 --bits-to-steal 1

<b>to decode the lossless payload from the image</b><br />
c-steg --decode-lossless --input-file ./my-photo-encoded.jpg

<b>to encode a lossy payload into a jpg use the following command</b><br />
c-steg --encode-lossy --input-file	./my-photo.jpg --output-file ./my-photo-encoded.jpg \
--text "this is the text i want to encode in the jpg"  --quality 92	--compressions 3

<b>to decode the lossy payload from the image</b><br />
c-steg --decode-lossless --input-file ./my-photo-encoded.jpg
