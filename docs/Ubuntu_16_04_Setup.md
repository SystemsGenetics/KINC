
# Development on Ubuntu 16.04

Use the following steps to setup KINC for development on Ubuntu 16.04:

## Step 1: Install Dependencies

Most of the dependencies are available as packages:
```bash
sudo apt install build-essential libgsl-dev libopenblas-dev libopenmpi-dev ocl-icd-opencl-dev
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

## Step 2: Download ACE and KINC

Clone the ACE and KINC repositories from Github.

```bash
git clone git@github.com:SystemsGenetics/ACE.git
git clone git@github.com:SystemsGenetics/KINC.git
```

## Step 3: Build ACE and KINC

Follow the ACE instructions to build ACE. If you install ACE locally then you must add ACE to the linker path:

```bash
# append to ~/.bashrc
export INSTALL_PREFIX="$HOME/software"
export CPLUS_INCLUDE_PATH="$INSTALL_PREFIX/include:$CPLUS_INCLUDE_PATH"
export LIBRARY_PATH="$INSTALL_PREFIX/lib:$LIBRARY_PATH"
export LD_LIBRARY_PATH="$INSTALL_PREFIX/lib:$LD_LIBRARY_PATH"
```

Build & install KINC:

```bash
cd build
qmake ../src/KINC.pro PREFIX=$INSTALL_PREFIX
make qmake_all
make
make qmake_all
make install
```

You should now be able to run KINC.

## (Optional) Use QtCreator

Select **File** > **Open File or Project** and then navigate in the file browser to the ACE directory and select the ACE.pro file. Navigate through configure setup. Repeat for KINC.
