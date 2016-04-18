Accelerated Computational Engine (ACE) is a GPU-enabled framework to simplify 
creation of GPU-capable appliations. .... <need more details>

ACE does not provide any data or analytical plugins.  ACE is a framework that
should be forked for other projects to develop their own data and analytical 
plugins, and release those projects as stand-alone software for use by their 
respective communities.

Dependencies
------------
1) Make 4.1
2) GCC 4.9
3) GPU drivers for your hardware
4) OpenCL


Prerequisite Installation for Ubuntu 16.04 LTS
----------------------------------------------
Ubuntu 16.05 has a compatible Make and GCC version. Installation of GPU 

Prerequisite Installation for Ubuntu 16.04 LTS
----------------------------------------------

1) To upgrade the make utility use the following:

  wget http://ftp.us.debian.org/debian/pool/main/m/make-dfsg/make_4.1-9_amd64.deb
  sudo dpkg --install make_4.1-9_amd64.deb 
  
2) Ubuntu 14.04 provides an older version of GCC. To update the following commands
should be executed to update to GCC v4.9:

  sudo add-apt-repository ppa:ubuntu-toolchain-r/test
  sudo apt-get update
  sudo apt-get install gcc-4.9 g++-4.9
  sudo update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-4.9 60 \
    --slave /usr/bin/g++ g++ /usr/bin/g++-4.9
    
3) GPU or graphics card drivers should be installed as per the manufacture's 
recommendations.  For NVidia it is recommended to download and install the
drivers from the manufacturers site.  
 
4) To install OpenCL:
  
  sudo apt-get install ocl-icd-opencl-dev

Compile ACE
------------
Inside of the 'src' directory. Exeucte the following

  make

Test ACE
--------
After compilation, a new 'run' directory is created at the same level
as the 'src' directory.  Test ACE by executing the 'unit' binary. For example,
after compiling in the 'src' directory:

  ../run/unit

If all tests passed you should see a message such as:
  
  144 unit test(s) passed <<<
  
After tests are executed you will be at the console prompt of ACE. You can exit
by typing 'quit'.

Running ACE
-----------
ACE does not provide any data or analytical plugins but it does have a 
functioning console and can identify existing accelerated computing 
infrastructure on the machine.  The ACE binary is in the 'run' directory. To
run ACE, execute the following

  ./run/ace
  
 
  