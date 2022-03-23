# dynamicAtlas

[![Documentation](https://img.shields.io/badge/Documentation-Link-blue.svg)](https://npmitchell.github.io/dynamicAtlas/)

> MATLAB toolkit for interfacing a library of fixed and dynamic datasets, with methods to aligning data in master timeline(s)

## Overview

**dynamicAtlas** is a MATLAB toolkit for interfacing with a library of fixed and live datasets of a morphological process, here focusing on the stage of development called 'germ band extension' in the fruit fly.

Suppose we want to create an lookup table of ALL data in the library. This is simple. First where is the data? Call that path ``atlasPath``

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



## Installation

A current MATLAB version is all you need to use dynamicAtlas.

Simply clone the repository as usual.
```bash
    git clone https://github.com/npmitchell/dynamicAtlas.git && cd dynamicAtlas
```

## Citing
TBA

## License

[MIT License](LICENSE)