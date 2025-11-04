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



