## CCD Package for [BEAST 2](https://www.beast2.org/)

Implementation of the conditional clade distribution (CCD), which offers an estimates of the posterior tree (topology) distribution learned from the sample, as well as associated tools and algorithms.
Three parametrisations of a CCD are implemented, namely, based on clade frequencies (CCD0), clade split frequencies (CCD1), and pairs of clade split frequencies (CCD2).
Furthermore, point estimators based on the CCDs are implemented,
which allows TreeAnnotator to produce better summary trees than via MCC trees (which is restricted to the sample).
The package further incldues tools for entropy-based rogue analysis and rogue removal.

See the [CCD-Research repository](https://github.com/CompEvol/CCD-Research) for information, scripts, and data of the associated research papers,
and see the instructions below on how to use the tools.

## Install the package

### Installing through BEAUti

CCD is a [BEAST2](http://beast2.org) package that requires BEAST 2 v2.7.6.
If you have not already done so, you can get BEAST 2 from [here](http://beast2.org).

To install CCD, it is easiest to start BEAUti (a program that is part of BEAST), and select the menu `File/Manage packages`. A package manager dialog pops up, that looks something like this:

![Package Manager](https://github.com/CompEvol/CCD/raw/master/doc/package_repos.png)

If the CCD package is listed, just click on it to select it, and hit the `Install/Upgrade` button.

If the CCD package is not listed, you may need to add a package repository by clicking the `Package repositories` button. A window pops up where you can click `Add URL` and add `https://raw.githubusercontent.com/CompEvol/CBAN/master/packages-extra-2.7.xml` in the entry. After clicking OK, the dialog should look something like this:

![Package Repositories](https://github.com/CompEvol/CCD/raw/master/doc/package_repos0.png)

Click OK and now CCD should be listed in the package manager (as in the first dialog above). Select and click Install/Upgrade to install.

### Install by hand

* Download the package from [here](https://github.com/CompEvol/CCD/releases/download/v0.0.1/CCD.package.v0.0.4.zip)
* Create CCD directory inside BEAST package directory
  * for Windows in Users\<YourName>\BEAST\2.X\VSS
  * for Mac in /Users/<YourName>\/Library/Application Support/BEAST/2.X/VSS
  * for Linux /home/<YourName>/.beast/2.X/VSS
  Here <YourName> is the username you use, and in “2.X” the X refers to the major version of BEAST, so 2.X=2.7 for version 2.7.6.
* Unzip the file `CCD.package.v0.0.1.zip` inside the CCD directory

## Build from code

* Get code for beast2, BeastFX and CCD repositories:
  * git clone https://github.com/CompEvol/beast2.git
  * git clone https://github.com/CompEvol/BeastFX.git
  * git clone https://github.com/CompEvol/CCD.git
* Run `ant install` from the CCD directory
  
## Usage

### Point Estimates

The tree topolgy point estimate provided by CCD's is called the CCD MAP tree and is best accessed through TreeAnnotator.
See also the general [tutorial](https://beast2.blogs.auckland.ac.nz/treeannotator/);
the steps are to start TreeAnnotator, and select `MAP (CCD0)` from the drop down box next to `Target tree type`:

![tree annotator](doc/treeannotator.png)

TreeAnnotator and thus the CCD point estimates can also be accessed from the terminal.
Simply but `CCD0` as option for the topology (`-topology`).
In full, the call from a terminal for Linux and OS X would be

```
/path/to/beast/bin/treeannotator -topology CCD0 input.trees output.tree
```

where `/path/to` the path to where BEAST is installed. For Windows, use

```
\path\to\BEAST\bat\treeannotator.bat -topology CCD0 input.trees output.tree
```

For CCD1 based point estimates select `MAP (CCD1)` from the drop down box in the GUI, or use `CCD1` instead of `CCD0` for the command line version.

### Rogue Analysis

A rogue analysis computes a rogue score for each clade, including each taxon, based on a given posterior sample of *binary* trees and some parameters and prints it to a csv file.
Furthermore, it annotates each clade in a CCD MAP tree with the respective rogue score.
For more information, see the [paper](https://www.biorxiv.org/content/10.1101/2024.09.25.615070v1) 
and for further information and example data see the [research paper repository](https://github.com/CompEvol/CCD-Research/tree/main/skeletonsAndRogues).

The tool can be accessed via the BEAST2 AppLauncher (see this [blog post](http://www.beast2.org/2019/07/23/better-apps-for-the-beast-appstore.html) for more information).
To start a rogue analysis on your sample of trees, you can use the following command in the terminal; when you omit any necessary parameter, a GUI will pop up and give you extra information.
```
/path/to/applauncher RogueAnalysis -trees treeInputFile.trees -burnin 10 -minProbability 0.1 -heightSettingStrategy CA -out outputFileWithoutEnding
```
In particular, note that the `out` filename will be used with `.csv` to print the full rogue score information and with `.trees` for the annotated MAP tree.
To get more information on the parameters, you can use:
```
/path/to/applauncher RogueAnalysis -help
```

### Skeleton Analysis

TBA
