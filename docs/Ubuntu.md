
# Installation on Ubuntu

Use the following steps to install KINC from source on Ubuntu 16.04:

## Step 1: Install Dependencies

Most of the dependencies are available as packages:
```bash
sudo apt install build-essential libgsl-dev libopenblas-dev libopenmpi-dev ocl-icd-opencl-dev liblapacke-dev
```

For device drivers (AMD, Intel, NVIDIA, etc), refer to the manufacturer's website.

### Qt (>=5.7)

Select a suitable version of Qt at http://download.qt.io/official_releases/qt and install Qt:

```bash
wget http://download.qt.io/official_releases/qt/5.7/5.7.1/qt-opensource-linux-x64-5.7.1.run
sh ./qt-opensource-linux-x64-5.7.1.run
```

If you install Qt locally then you must add Qt to the executable path:

```bash
# append to ~/.bashrc
export QTDIR="$HOME/Qt/5.7.1/gcc_64"
export PATH="$QTDIR/bin:$PATH"
```


## Step 2: Download and Build ACE
Determine which [version of ACE](https://github.com/SystemsGenetics/ACE/releases) you want to build and set the following environment variable:
```
export ACE_VERSION=v3.0.2
```

Next, clone the ACE source code:
```
git clone https://github.com/SystemsGenetics/ACE.git
cd ACE/build
git checkout $ACE_VERSION
```

By default, ACE will try to install itself into `/usr/local`. To ease installation, set this to be the INSTALL_PREFIX environemtn variable:
```
export INSTALL_PREFIX=/local/software
```

Suppose you desire to install ACE in a different directory (e.g. `/local/software`) You can adjust the installation directory by changing this envrionement variable:
```
export INSTALL_PREFIX=/local/software
```

Now, within the `ACE/build` directory run the following to build the ACE libraries:
```
qmake ../src/ACE.pro PREFIX=$INSTALL_PREFIX/ACE-$ACE_VERSION
make qmake_all
make
make qmake_all
make install
```
This will install ACE into the directory specified by the `INSTALL_PREFIX` variable in a directory named with the ACE version.

## Step 3: Build KINC
To ease installation will set the default version for both ACE and KINC. The version of ACE used here should be the same specified above when compiling ACE.  For KINC, determine which [version of KINC](https://github.com/SystemsGenetics/KINC/releases) you want to build and set the environment variable:
```
export ACE_VERSION=v3.0.2
export KINC_VERSION=v3.2.2
```

Next, clone the KINC repository.
```bash
git clone https://github.com/SystemsGenetics/KINC.git
cd KINC/build
git checkout $KINC_VERSION
```

By default, KINC will try to install itself into `/usr/local`. To ease installation, set this to be the INSTALL_PREFIX environemtn variable:
```
export INSTALL_PREFIX=/usr/local
```

Suppose you desire to install ACE in a different directory (e.g. `/local/software`) You can adjust the installation directory by changing this envrionement variable:
```
export INSTALL_PREFIX=/local/software
```

Before you can build KINC, the compiler must be able to find the ACE libraries.  Several environment variables help with this:
```
export PATH="$INSTALL_PREFIX/ACE-$ACE_VERSION/bin:$PATH"
export LD_LIBRARY_PATH="$INSTALL_PREFIX/ACE-$ACE_VERSION/lib:$LD_LIBRARY_PATH"
export LIBRARY_PATH="$INSTALL_PREFIX/ACE-$ACE_VERSION/lib:$LIBRARY_PATH"
export CPATH="$INSTALL_PREFIX/ACE-$ACE_VERSION/include:$CPATH"
export C_INCLUDE_PATH="$INSTALL_PREFIX/ACE-$ACE_VERSION/include:$C_INCLUDE_PATH"
export CPLUS_INCLUDE_PATH="$INSTALL_PREFIX/ACE-$ACE_VERSION/include:$CPLUS_INCLUDE_PATH"
export OBJC_INCLUDE_PATH="$INSTALL_PREFIX/ACE-$ACE_VERSION/include:$OBJC_INCLUDE_PATH"
```

Now build & install KINC:

```bash
cd KINC/build
qmake ../src/KINC.pro PREFIX=$INSTALL_PREFIX/KINC-$KINC_VERSION
make qmake_all
make
make qmake_all
make install
```
To run KINC you must update the `LD_LIBRARY_PATH` in your `~/.bashrc` file.  Use the following command to get the exact text you need to add.

```
echo "export LD_LIBRARY_PATH=\"$INSTALL_PREFIX/ACE-$ACE_VERSION/lib:$INSTALL_PREFIX/KINC-$KINC_VERSION/lib:\$LD_LIBRARY_PATH\""
echo "export PATH=\"$INSTALL_PREFIX/ACE-$ACE_VERSION/bin:$INSTALL_PREFIX/KINC-$KINC_VERSION/bin:\$PATH\""
```
Add the resulting text to the end of your ~/.bashrc file.
You should now be able to run KINC
