# MorphodynamicAtlas: an approach for integrating dynamic datasets

[![Documentation](https://img.shields.io/badge/Documentation-Link-blue.svg)](https://npmitchell.github.io/dynamicAtlas/)

> MATLAB toolkit for interfacing a library of fixed and dynamic datasets, with methods to aligning data in master timeline(s)

## Overview

**dynamicAtlas** is a MATLAB toolkit for interfacing with a library of fixed and live datasets of a morphological process, here focusing on the stage of development called 'germ band extension' in the fruit fly.

Suppose we want to create a lookup table of ALL data in the library. This is simple. First where is the data? Call that path ``atlasPath``

```matlab
	atlasPath = '/path/to/the/data/dynamicAtlas/'
	da = dynamicAtlas.dynamicAtlas(atlasPath, {'WT'}) ;
```

Now ``da`` is an instance of the dynamicAtlas class with access to all the data through its methods.

What if we just want a subset of the WT library? We can look for just certain geonotypes (parent directories) and/or for certain labels (subdirectories in the library). Let's grab all the data stained for the transcription factor Runt in embryos that have a wild-type genotype. We can create an atlas of this subset of the data like:

```matlab
	options = struct() ;
	options.labels = { 'Runt'} ;
	da = dynamicAtlas.dynamicAtlas(atlasPath, {'WT'}, options) ;
```

Then we could slice all of this data with a given timestamp range
(timestamp-delta, timestamp+delta) via

```matlab
    qs = da.findTime(timestamp, delta)
```

Then we have a queriedSample object called ``qs``. It has some interesting methods like ``getData()`` and ``getPIV()`` which load the pullback images of the tissue or the velocities for us to analyze it. To execute these methods, simply run
``qs.getData()`` or ``qs.getPIV()``.
We can also ask for things like the mean velocity across all the data in the queried sample using ``qs.getMeanPIV()``.

Alternatively, we could slice the atlas data for only frames of live datasets of a given genotype and a given label with a given timestamp range:

```matlab
    qs = da.findDynamicGenotypeLabelTime('WT', 'Runt', timestamp, delta)
```

Similar methods to check out are ``da.findEmbryo('202006261115')``, ``da.findStaticGenotypeLabel()`` and similar.

## System and software requirements

- OS: Tested on Mac OS Sequoia 15.1, Windows 10, and Linux Ubuntu 18.04.6 LTS, but will also work on other operating systems, as long as these systems support MATLAB installation.

- Software: Tested on MATLAB R2024a and R2024b, but will work on any recent MATLAB release. No other software dependencies are required.
  
- Hardware: No non-standard hardware is required. Performance is optimized when RAM is at least 18 GB.

## Installation

A current MATLAB version is all you need to use dynamicAtlas. 
Simply clone the repository as usual. Installation should occur within 1 minute.
```bash
    git clone https://github.com/npmitchell/dynamicAtlas.git && cd dynamicAtlas
```

- Demo Dataset Access: available for download on Zenodo at: https://doi.org/10.5281/zenodo.14792464
  
- Full Atlas Data Access: available for download on Dryad at: https://datadryad.org/stash/dataset/doi:10.25349/D9WW43 (works best on the Firefox and Google Chrome
browsers)

- Download the data from the repository using the link(s) provided above, and unzip it.  

- Open MATLAB, set the variable ‘atlasPath’ as the path to the directory containing the unzipped data, and set the variable 'codePath' as the path of the dynamicAtlas code folder, as described in the Demo Script provided in GitHub repository: “demo_dynamicAtlas_functionality.m”. See demo tutorial PDF "DynamicAtlas_Demo_Tutorial.pdf" in the repository, which is also included in the Supplementary Information of the manuscript, for a tutorial walkthrough that runs this script on a demo dataset.

## Demo
A demo script is provided on the GitHub repository as “demo_dynamicAtlas_functionality.m”, and also included as part of the Supplementary Information. This script is the basis of the walkthrough included in the tutorial mentioned above, "DynamicAtlas_Demo_Tutorial.pdf". The tutorial includes instructions for initializing the code and running the software on the data, as well as screenshots of the expected outputs that correspond to each section of the code. Demo run time on the demo dataset should be less than 1 hour on a standard computer.

## Citing
N. P. Mitchell*, M. F. Levebvre*, V. Jain-Sharma*, N. Claussen, M. K. Raich, H. J. Gustafson, A. R. Bausch, S. J. Streichan. “Morphodynamic atlas of Drosophila development.” bioRxiv 10.1101/2022.05.26.493584 (2022). 

## License

[MIT License](LICENSE)
