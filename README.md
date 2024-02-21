## CCD Package for [BEAST 2](beast2.org/)

Implementation of the conditional clade distribution (CCD) with two parametrisations,
namely, based on clade split frequencies (CCD1) and clade frequencies (CCD0).
Furthermore, point estimators based on the CCDs are implemented,
which allows TreeAnnotator to produce better summary trees than via MCC trees (which is restricted to the sample).

## Install 

* Download the package from [here](https://github.com/CompEvol/CCD/releases/download/v0.0.1/CCD.package.v0.0.1.zip)
* Create CCD directory inside BEAST package directory
  * for Windows in Users\<YourName>\BEAST\2.X\VSS
  * for Mac in /Users/<YourName>\/Library/Application Support/BEAST/2.X/VSS
  * for Linux /home/<YourName>/.beast/2.X/VSS
  Here <YourName> is the username you use, and in “2.X” the X refers to the major version of BEAST, so 2.X=2.1 for version 2.1.3.
* Unzip the file `CCD.package.v0.0.1.zip` inside the CCD directory

## Build from code

* Get code for beast2, BeastFX and CCD repositories:
  * git clone https://github.com/CompEvol/beast2.git
  * git clone https://github.com/CompEvol/BeastFX.git
  * git clone https://github.com/CompEvol/CCD/releases/download/v0.0.1/CCD.package.v0.0.1.zip
* Run `ant install` from the CCD directory
  
## Usage

Start TreeAnnotator, and select `Conditional Clade Distribtion0` from the drop down box next to `Target tree type`:

![tree annotator](doc/treeannotator.png)


From a terminal, to use CCD0 for Linux and OS X run

```
/path/to/beast/bin/treeannotator -topology CCD0 input.trees output.tree
```

where `/path/to` the path to where BEAST is installed. For Windows, use

```
\path\to\BEAST\bat\treeannotator.bat -topology CCD0 input.trees output.tree
```


