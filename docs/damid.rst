DamID principle
---------------

DamID is a method that identifies genomic binding sites of a protein of interest (POI). In contrast to ChIP-seq, no antibodies against the POI are required.

The method is based on fusing the POI to a DNA adenine methyltransferase (Dam) from *E. coli*. Dam methylates adenine bases in GATC motifs, which are not present in the eukaryotic genome. 

To identify the genomic binding sites of the POI, the Dam-fusion protein is expressed in cells. Over time, the Dam-fusion protein will bind to the chromatin and methylate the GATC motifs in the vicinity of the binding sites. These methylation marks effectively function as a recording of the history of the chromatin binding of the POI.

For quantification, DNA is isolated and methylated DNA is cloned and subsequently sequenced.


.. figure:: images/damid.png
   :align: center
   :width: 1000
   
   Figure adapted from Van den Ameele et al. 2019 Current Opinion in Neurobiology


Experimental considerations
---------------------------

#. Make sure that you have sequence validated your Dam expressing plasmids prior to your experiment. Cryptic expression of Dam can cause toxicity in *E. coli* and can cause selection for Dam inactivating mutations.