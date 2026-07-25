Example: 1
----------

Generate a molecule object from .mol2 file and shift atoms:

.. code:: python

    import mooonpy
	
    molecule = mooonpy.molspace('detda.mol2')
    molecule.atoms.shift(x=5, y=10, z=30)

