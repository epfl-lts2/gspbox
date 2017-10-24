Graph Signal Processing Toolbox (GSPBox)
======

The official Graph Signal Processing Toolbox (GSPBox)

Running the toolbox
==

In Matlab type

        gsp_start
        
as the first command from the installation directory. This will set up the correct paths.

The gsp_start command will add all the necessary subdirectories (so
please don't add these manually), and it will print a statement
telling you which backend you are currently using.

Compiling the toolbox
===

In Matlab type

        gsp_make
        
as the first command from the installation directory. This will try to compile the mex interface for 
the third party components: AMD and LDL. 

You need the library stdc++ to compile AMD and LDL. On Ubuntu you can type:

        sudo apt-get install libstdc++6


Installing third party software
====

In Matlab type 

        gsp_install

as the first command from the installation directory. This will try to download and install third party 
software.


Documentation
===

You can find the complete online documentation here : http://lts2research.epfl.ch/gsp/doc/start/
  
