# Coding for KINC in Windows
To begin development work in KINC on Windows the following steps can be used to setup your coding environment:

## Step 1:  Download QT Creator
KINC is build with the [ACE library](https://github.com/systemsgenetics/ACE) which uses QT libraries for display and data types. To develop with KINC, the QTCreator programming should be installed. You can obtain the free, community version here:  https://www.qt.io/download

## Step 2: Download ACE and KINC
Use the Git Bash tool to clone the ACE and ACE repositories and switch to the develop branch

```bash
git clone git@github.com:SystemsGenetics/ACE.git
cd ACE
git checkout develop 

cd ../
git clone git@github.com:SystemsGenetics/KINC.git
cd KINC
git checkout develop
```

## Step 3: Setup ACE and KINC in QTCreator
Within QTCreater click **File** > **New File or Project** and click **Import Existring Project** and in the box that appears, select the path to the directory where ACE was cloned.  Be sure to select the ACE/src directory.  Follow the dialog prompts until the project is imported.  Do the same with KINC.
