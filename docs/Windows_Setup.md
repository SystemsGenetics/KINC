# WARNING
Windows is temporarily not supported because there is no OpenMPI library for the windows platform. This will be changed in the future when MPI becomes optional at compile time.

# Coding for KINC in Windows
To begin development work in KINC on Windows the following steps can be used to setup your coding environment:

## Step 1:  Download QT Creator
KINC is build with the [ACE library](https://github.com/systemsgenetics/ACE) which uses QT libraries for display and data types. To develop with KINC, the QTCreator programming should be installed. You can obtain the free, community version here:  https://www.qt.io/download.

When installing QTCreater, be sure to install the QT 5.7 (or later) and MinGW 5.3.0 (or later) compiler.  You will be able to do this using the installation dialogue.

## Step 2: Download ACE and KINC
Use the Git Bash tool to clone the ACE and ACE repositories and switch to the develop branch

```bash
git clone git@github.com:SystemsGenetics/ACE.git
git clone git@github.com:SystemsGenetics/KINC.git
```

## Step 3: Setup ACE and KINC in QTCreator
First add ACE. Click **File** > **Open File or Project** and then navigate in the file browser to the ACE directory and click the ACE.pro file. Navigate through configure setup. Do the same for KINC.
