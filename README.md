# libbinary

This is the library of functions used by my various close binary star analysis and visualisation codes from my PhD.

I am resurrecting this initially as an exercise to test various coding tools. Sibling repos depending on this will follow.

Early imports may be non-functional or have out-of-date docs - this code is from 2004.

It should be working when I merge it into master.

# Build and install libbinary

Install dependencies, e.g. In Ubuntu 14.04

    sudo apt-get install build-essential

Set the install location environment variable <tt>LIBBINARY</tt> and add to the library path, e.g.

    export LIBBINARY=/usr/local/libbinary

to your <tt>.bashrc</tt>

Compile and install

    make
    sudo -E make install

That should be it. Won't do much until I add something which depends on it.

One of my aims with creating this is to put together a modern CMake makefile and proper build configuration to remove the need for the ancient methods above. But first I will get it working again the old way.

# Original uncorrected documentation

Contents:
* libbinary:
* idl binary routines:
* stream_example

1). Setup IDL so it can find everything.

Add the following lines to your .cshrc or equivalent:

# IDL
setenv IDL_STARTUP ${HOME}/.idlrc
setenv IDL_DIR /usr/local/rsi/idl
setenv IDL_PATH \+${HOME}/idl
setenv IDL_PATH ${IDL_PATH}:\+$IDL_DIR/lib:\+$IDL_DIR/examples
alias idl "xterm -fg yellow -bg black -tn xterm -T 'IDL' -e $IDL_DIR/bin/idl &"

The IDL_DIR should point to the idl directory in your local installation of
IDL. You can leave out the alias command if you like - it just makes IDL run
in a separate window with a helpful title.

2) Installing libbinary.

Add these lines to you .cshrc (after the IDL stuff)
where path_to/libbinary should be changed to the appropriate path to
the libbinary directory.

# Binary library
setenv BINARY path_to/libbinary/
if ${?LD_LIBRARY_PATH} then
    setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:${BINARY}
else
    setenv LD_LIBRARY_PATH ${BINARY}
endif
setenv IDL_PATH ${IDL_PATH}:\+${BINARY}/idl

Run your .cshrc script or equivalent: source ~/.cshrc

cd to the libbinary directory.
type 'make install'.

This should compile and setup the binary library.

To your .idlrc file in your home directory (create it if it's not there), add
the following line:

DEFSYSV, '!LIBBINARY', '/local/dj/code/libbinary/libbinary.so'

where the path should be to the full path to libbinary.so file.

3) Testing it.

run IDL then '.run stream_example' (the included file).
This should calculate and plot the trajectory then velocity
of a stream. stream_example is commented to explain how to use it.
