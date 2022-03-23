.. dynamicAtlas documentation master file, created by
   sphinx-quickstart on Mon Mar 21 12:01:45 2022.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to dynamicAtlas's documentation!
========================================

DynamicAtlas is a MATLAB toolkit for interfacing with a library of fixed and live datasets of a morphological process, here focusing on the stage of development called 'germ band extension' in the fruit fly.

Suppose we want to create an lookup table of ALL data in the library. This is simple. First where is the data? Call that path ``atlasPath``


.. code-block:: matlab

	atlasPath = '/path/to/the/data/dynamicAtlas/'
	da = dynamicAtlas.dynamicAtlas(atlasPath, {'WT'}) ;``

Now ``da`` is an instance of the dynamicAtlas class with access to all the data through its methods.

What if we just want a subset of the WT library? We can look for just certain geonotypes (parent directories) and/or for certain labels (subdirectories in the library). Let's grab all the data stained for the transcription factor Runt in embryos that have a wild-type genotype. We can create an atlas of this subset of the data like:


.. code-block:: matlab

	options = struct() ;
	options.labels = { 'Runt'} ;
	da = dynamicAtlas.dynamicAtlas(atlasPath, {'WT'}, options) ;``

Then we could slice all of this data with a given timestamp range
(timestamp-delta, timestamp+delta) via

.. code-block:: matlab
    qs = da.findTime(timestamp, delta)``

Then we have a queriedSample object called ``qs``. It has some interesting methods like ``getData()`` and ``getPIV()`` which load the pullback images of the tissue or the velocities for us to analyze it. To execute these methods, simply run
``qs.getData()`` or ``qs.getPIV()``.
We can also ask for things like the mean velocity across all the data in the queried sample using ``qs.getMeanPIV()``.

Alternatively, we could slice the atlas data for only frames of live datasets of a given genotype and a given label with a given timestamp range:

.. code-block:: matlab
    qs = da.findDynamicGenotypeLabelTime('WT', 'Runt', timestamp, delta)``

Similar methods to check out are ``da.findEmbryo('202006261115')``, ``da.findStaticGenotypeLabel()`` and similar.


Examples of using dynamicAtlas
------------------------------

.. toctree::
   :maxdepth: 1
   :caption: Example scripts:

   example_usage_genes
   example_usage_piv


Classes and Folders
--------------------------------------

.. toctree::
   :maxdepth: 1
   :caption: Classes in +dynamicAtlas:
   
   dynamicAtlas
   lookupMap
   queriedSample
   

.. toctree::
   :maxdepth: 1
   :caption: Folders of functions:
   
   basics
   curve_handling
   data_handling
   external		
   matchTime		
   piv_handling		
   plotting		
   scripts			
   statistics	
   stripeExtraction
   tiff_handling

Some details on how dynamicAtlas works
--------------------------------------

We saw how to make a dynamicAtlas class instance and a queriedSample class instance.
There is one more class in this toolkit that helps it all run smoothly called a lookupMap. You don't have to use the lookupMap if you want to keep things simple, but in case you want to look under the hood, let's take a look. 
You can instantiate the map using the ``lookup`` method of dynamicAtlas:

.. code-block:: matlab
    mapWT = da.lookup('WT') ;

Now mapWT has methods:

.. code-block:: matlab
    methods(mapWT)

To find embryos with a Runt stain, you can say:

.. code-block:: matlab
    runts = findGenotypeLabel('WT', 'Runt')

Alternatively, we can take a peek within the lookup map of da to get
all Runt labels, regardless of genotype:

.. code-block:: matlab
    runts = mapWT.findLabel('Runt')

or even index the map stored in da's lookup map:

.. code-block:: matlab
    runts = mapWT.map('Runt')

To find embryos with a t=10 +/- 2 min

.. code-block:: matlab
    snaps = mapWT.findTime(10, 2)

To find Runt stains with a t=10 +/- 2 min

.. code-block:: matlab
    runtsnaps = mapWT.findLabelTime('Runt', 10, 2)


Master timeline generation
--------------------------

Fixed samples are timestamped according to their morphology against live datasets. 
If we want to build such a master timeline from scratch, we can do so as follows.

Use dynamic datasets within the lookupMap to build master timeline.
To control how this is performed, toggle da.timeLineMethod:

.. code-block:: matlab

    da.makeMasterTimeline('WT', 'Runt')

Timestamp other data against the master timeline.
To control how this is performed, toggle da.timeStampMethod.

.. code-block:: matlab
	da.timeStampStripe7('WT', 'Runt')


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
