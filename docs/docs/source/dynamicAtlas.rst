dynamicAtlas
============


.. toctree::
   :maxdepth: 4



Example Usage
-------------

Suppose we want to grab all the Wild-type Runt data:

``>> da = dynamicAtlas.dynamicAtlas('/path/to/data/', {'WT', 'Runt'})``

Then we could slice all of this data with a given timestamp range
(timestamp-delta, timestamp+delta)

``>> qs = da.findTime(timestamp, delta)``

Then we have a queriedSample object

``mapWT = da.lookup('WT') ;``

Now mapWT has methods:

``>> methods(mapWT)``

To find embryos with a Runt stain, you can say:

``>> runts = findGenotypeLabel('WT', 'Runt')``

Alternatively, we can take a peek within the lookup map of da to get
all Runt labels, regardless of genotype:

``>> runts = mapWT.findLabel('Runt')``

or even index the map stored in da's lookup map:

``>> runts = mapWT.map('Runt')``

To find embryos with a t=10 +/- 2 min

``>> snaps = mapWT.findTime(10, 2)``

To find Runt stains with a t=10 +/- 2 min

``>> runtsnaps = mapWT.findLabelTime('Runt', 10, 2)``

Use dynamic datasets within the lookupMap to build master timeline.
To control how this is performed, toggle da.timeLineMethod:

``>> da.makeMasterTimeline('WT', 'Runt')``

Timestamp other data against the master timeline.
To control how this is performed, toggle da.timeStampMethod:

``>> da.timeStampStripe7('WT', 'Runt')``

To Do:
------
handle PIV-based timeline creation
	
Full Contents:
--------------

.. automodule:: @dynamicAtlas
    :show-inheritance:
    :members:

	
	
Indices and tables
------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
