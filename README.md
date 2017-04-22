# Taxalotl
This is python package for merging taxonomies to continue
the releases in the OTT series of taxonomies.
This package might end up just being a sub-package of peyotl,
or die a young death.

The code is based on parts of the
    [reference-taxonomy](https://github.com/OpenTreeOfLife/reference-taxonomy)
    repository, which was written primarily by 
    [Jonathan A. Rees](https://github.com/jar398) as 
    a part of the [NSF](https://www.nsf.gov/)-funded 
    [Open Tree of Life](https://tree.opentreeoflife.org)
    project.

## installation
`virtualenv` is always recommended for experimental
    Python packages.

Depends on:
  * `requests`
  * some utilities that are only on a dev branch of `peyotl`

So, the easiest way to install right now is:

    ./install.sh

which will:
   1. create a virtualenv called `venv`,
   2. install prerequisites in it,
   3. install the correct version of `peyotl` using the "develop"
    command to pip (to make a symlink), and
   4. install the taxalotl package (also using the "develop" 
   command)

See [peyotl docs](http://opentreeoflife.github.io/peyotl/) for
info about the config files of peyotl.
These affect the logging message handling of Taxalotl.




## Thanks and stuff
Thanks to the US NSF
<a href="https://www.nsf.gov/"><img src="./doc/nsf1.jpg" alt="NSF logo" /></a>

The package is named after [_Ambystoma mexicanum_](https://en.wikipedia.org/wiki/Axolotl) ...
and the fact that there are a lot o' taxa...
and the whole Open Tree of Life thing...

