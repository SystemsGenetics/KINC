
# Coding for KINC on Ubuntu 16.04
To begin development work in KINC on Ubuntu 16.04 the following steps can be used to setup your coding environment:

## Step 1: Install Dependencies
First we need to instal the OpenCL development libraries.  

```bash
sudo apt install ocl-icd-libopencl1
sudo apt install opencl-headers
sudo apt install clinfo
sudo apt install ocl-icd-opencl-dev
```
For Intel or NVidia drivers see the respective drivers.

Next, we need to install the OpenMPI library

```bash
sudo apt-get install openmpi-bin libopenmpi-dev openmpi-common
```

## Step 2: Download ACE and KINC
Use Git to clone the ACE and ACE repositories and switch to the develop branch

```bash
git clone git@github.com:SystemsGenetics/ACE.git
cd ACE
git checkout develop 

cd ../
git clone git@github.com:SystemsGenetics/KINC.git
cd KINC
git checkout develop
```
### Step 3: Compile ACE
