Internal storage
================

Default top path (controlled by t) is ~/.opentreeoflife/taxalotl


```
|=== normalized
|             \=== ${resource}-${version}
|                           |=== about.json
|                           |=== details.json
|                           |=== synonyms.tsv
|                           \=== taxonomy.tsv
|=== partitioned
|             |=== Life
|             |             |=== Archaea
|             |             |             \=== __inputs__
|             |             |                           |=== ${resource}-${version}
|             |             |                           |             |=== __roots__.json
|             |             |                           |             |=== taxonomy.tsv
|             |             |                           \             \=== synonyms.tsv
|             |             |=== Bacteria
|             |             |             \=== __inputs__
|             |             |                           \=== ${resource}-${version}
|             |             |                                         |=== __roots__.json
|             |             |                                         |=== taxonomy.tsv
|             |             |                                         \=== synonyms.tsv
|             |             |=== Eukaryota
|             |             |             |=== __inputs__
|             |             |                           \=== ${resource}-${version}
|             |             |                                         |=== __roots__.json
|             |             |                                         |=== taxonomy.tsv
|             |             |                                         \=== synonyms.tsv
|             |             \=== __misc__
|             |                           \=== __inputs__
|             |                                         \=== ${resource}-${version}
|             |                                                       |=== __accum_des__.json
|             |                                                       |=== taxonomy.tsv
|             |                                                       \=== synonyms.tsv
|             |=== __mapping__.json
|             |=== __separator_names__.json
|             \=== __separator_names_to_dir__.json
|=== raw
|             \=== ${resource}-${version}
|                           \===  variable filenames depending on resource
|=== resources
|             |=== ${resource}.json
|             \=== ${resource}.json
|=== taxalotl.conf

```

## Partition command
Partitioning a resource breaks it into some hard-coded groups (see below) to make the `taxonomy.tsv` files more maneagable. 
Partitioning will result in many of the taxa being moved into one of the clades that is daughter of the current level.
The (typically paraphyletic) assemblage of other taxa not included in one of the name clades will be placed into the appropriate
  folder (for the resource) in the `__misc__` diretory at the parent level.
For each level of the hierarchy, the resources are placed into an `__inputs__` directory inside a subdirectory
  named with the full, versioned name of the resource.

For every level except the "Life" level, a `__roots__.json` file will contain a JSON version of any taxon/taxa that root that subtree. 
Multiple roots can exist for an input resource that does not formally name the level of the taxon.
For example if an input had all of the families for some suborder, but lacked the suborder itself, if the resource were partitioned at the level of that suborder, then all of the family roots would be listed in the roots file.
The roots are also in the `taxonomy.tsv` file.

Similarly, the `__accum_des__.json` contains and object with JSON records of the roots of the subtree that the resource included and have been partitioned below this level.
So, a parent level's `__accum_des__.json` is essentially a concatenation of it's daughters' `__roots__.json` file.
TODO: This file is probably not necessary because it is just redundant info.



### Current hard-coded partition breaks
```
Life
  Archaea
  Bacteria
  Eukaryota
    Haptophyta
    SAR
    Fungi
    Archaeplastida
      Chloroplastida
      Glaucophyta
      Rhodophyta
    Metazoa
      Bryozoa
      Chordata
      Cnidaria
      Ctenophora
      Mollusca
      Nematoda
      Platyhelminthes
      Porifera
      Annelida
      Arthropoda
        Arachnida
        Malacostraca
        Insecta
          Coleoptera
          Diptera
          Hymenoptera
          Lepidoptera
```


