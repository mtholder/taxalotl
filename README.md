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

Taxalotl depends on:
  * `requests`
  * `beautifulsoup4`
  * some utilities that are only on the `taxalotl` branch of `peyotl`

So, the easiest way to install right now is:

    ./install.sh

which will:
   1. create a virtualenv called `env` sister to this (taxalotl) dir,
   2. install prerequisites in it,
   3. clone (as sister to this dir) and install the correct version of `peyutil` and `peyotl` 
   4. install the taxalotl package (also using the "develop" 
   command)

The script ends with some comments about the actions you need
to take to configure the package.

See [peyotl docs](http://opentreeoflife.github.io/peyotl/) for
info about the config files of peyotl.
These affect the logging message handling of Taxalotl.


## Usage
The `taxalotlcli` script provides the command-line interface which
is broken up into several commands.
The command-level documentation is below.
See the [Tutorial.md](./Tutorial.md) for an overview of usage.

The syntax `${x}` in the documentation below refers to
the value of some variable (`x` in this case),
that should be one of the configuration variables specified in
the `taxolotl.conf` file.
    

### status command
`taxalotlcli status` reports on the status of each "resource".

`taxalotlcli status ncbi` reports just on the status of the
    ncbi resource.

### download command
`taxalotlcli download ID` downloads the archive for
    the `ID` resource into the `${raw}` directory if that
    archive is not present.

### unpack command
`taxalotlcli unpack ID` unpacks the archive for
    the `ID` resource to the `${raw}/ID`.
Downloads the archive if necessary.

### normalize command
`taxalotlcli normalize ID` unpacks the raw archive
for `ID` from `${raw}/ID` into the 
[OTT Interim Taxonomy](https://github.com/OpenTreeOfLife/reference-taxonomy/wiki/Interim-taxonomy-file-format)
format in `${normalized}/ID`
Unpacks the raw archive if necessary.

An "extra" `details.json` file may also be written with 
more information about the normalization process.
This information is "extra" in the sense that it was not
emitted by the reference-taxonomy repo's version of the code.


## Structure
### Resources directory
That dir should hold descriptions of the taxonomies that
    are sources of information.
The file `.merged.json` in that directory is automatically
    generated as the union of the fields in all of the
    other files in the directory;
so you should not edit that file by hand.
Most info should probably just go in `resources.json`, but
    you can also put info in a file with the name
    `<resource-id>.json` to make the resources list more
    manageable.

### Resources and id aliasing
Tags used to describe the resources are still in flux.
Every resource should have either:
  * a `resource_type` 
    property (with value "external taxonomy",
    "open tree taxonomy idlist", or "open tree taxonomy"), or
  * an `inherits_from` property with the value
    corresponding to the ID of another resource.

So, there is a resource with `id="ncbi"` that is intended
    to hold info that applies to any version of the NCBI
    Taxonomy.
Specific dated snapshots of NCBI inherit from either that
    "base" resource or the previous snapshot.
This will create a directed linear graph. 
If you use the base id (`ncbi`) in a command that operates
    on a concrete resource, `taxalotl` will assume that you
    mean the latest version of that resource.

## Extras
If you `source completion.sh` in your `bash` session and you have
    the top-level directory on your `PATH` then you'll have some
    cool tab completion of commands and options.

## Thanks and stuff
Thanks to the US NSF
<a href="https://www.nsf.gov/"><img src="./doc/nsf1.jpg" alt="NSF logo" /></a>

The package is named after [_Ambystoma mexicanum_](https://en.wikipedia.org/wiki/Axolotl) ...
and the fact that there are a lot o' taxa...
and the whole Open Tree of Life thing...

However, see https://www.youtube.com/watch?v=Ka0Fj6P3T-w on the correct pronunciation
  of "Axolotl" (and some other Nahuatl words).

It might become `Taxolotl` instead of `Taxalotl`.
