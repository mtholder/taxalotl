# Taxalotl tutorial
See [README.md](./README.md) for general documentation.

## Installation and configuration
The following assumes that you are running bash (or some bash-compatible) shell:

    git clone https://github.com/mtholder/taxalotl.git
    cd taxalotl
    bash install.sh
    source env/bin/activate
    export PEYOTL_CONFIG_FILE="$PWD/recommended-peyotl-conf.ini"
    export PATH="${PATH}:${PWD}"
    source completion.sh

This should set up a python 3 virtualenv inside an `env` directory
and set up a reasonable configuration.

## Default configuration
The default scratch directory for is `~/.opentreeoflife/taxalotl` (this
will be referred to as TAXALOTL_SCRATCH but that is just a name used in documentation,
not an environmental variable that you need).
The default configuration file location is `~/.opentreeoflife/taxalotl/taxalotl.conf`.
When `taxalotl` becomes more stable, you may want to switch 
the `crash_with_stacktraces` setting in `~/.opentreeoflife/taxalotl/taxalotl.conf`
to `false`, but you probably want to leave it set to true now so that you can
submit more useful bug reports.

## Grokking taxalotl
The script `taxalotlcli` should be on your PATH, this is the front-end to taxalotl.
Sourcing the `completion.sh` script (as was done above) should enable tab-completion of
    commands and options in taxalotl; this makes the front-end **much easier** to work with.


taxalotl uses the [OTifacts](https://github.com/mtholder/OTifacts) repository for information
about the taxonomic resources that Open Tree uses.
taxalotl uses its scratch directory to manage you local copies (and manipulations) of these
resources.
In the context of taxalotl, a "resource" is almost always a taxonomy.

    taxalotlcli status

reports on the state of "your" local copies of the taxonomic resources. If you run this
immediately after installation, your taxalotl resource directory will be empty.
The resources directory is where taxalotl manages notes about what resources could be relevant.
It is a copy of the relevant information in the OTifacts repo.  

If there is any chance that your local resources view lags behind the OTifacts repo on the web,
    you'll want to run:

    taxalotlcli pull-otifacts

to populate your local clone of the OTifacts repo and process for access via taxalotl.
This action fills the TAXALOTL_SCRATCH/resources directory with JSON files about the
resources.
After running that

    taxalotlcli status

should report on the status of your local copies of each of the taxonomy-oriented artifacts
that are listed in OTifacts.
You can use:

    taxalotlcli status ott

to see info just on OTT, or you can substitute a series of other resource IDs in place
of "ott" in that command to a more targeted status report.

The goal of taxalotl is to manage these taxonomic resources for you in a series of flat
files for ease of exploration.
On each machine that you are using taxalotl on, you'll need to execute operations:
  1. download,
  2. unpack,
  3. normalize, and 
  4. partition
for each of the resources to download the archive, unpack it, normalize it (convert it to
Open Tree's taxonomy format), and partition it into major "separations" of the taxonomy.

The first time that you go through this process it may take a while, 
because the taxonomies are large and Python is not a fast language.
The motivation of taxalotl is to work on updates to smaller pieces of the taxonomies, assuming
that the current separation that has been used in OTT correctly partitions the input taxonomies.
So, when comparing new information on subsets of the taxonomies, the fact that Python is slow
should not be a major concern.
If you are working on multiple machines, you can just sync the TAXALOTL_SCRATCH directories 
because that is stage that taxalotl-cli is working.
Thus after partitioning all in the inputs, you have done most of the slow work, and you might
want to sync.

While you can execute each of the aforementioned operations as a separate command, taxalotl
does know that partitioning requires a normalized file, which requires an unpacked archive, 
which requires downloading...

So,

    taxalotlcli partition ott

is sufficient to trigger download, unpack, normalize, and partitions actions on the latest
version of OTT.

**NOTE** taxalotl will try to show you the license information that it knows about 
when you trigger a download, but it is your responsibility to make sure that you are abiding
by the terms of any license.
So you will usually have to answer a prompt before the download
starts.

Currently taxalotl has a somewhat quirky mix of information that is hard-coded based on
OTT version 3.0 (the start point for the taxalotl updating of taxonomy) and
information that it gleans from the taxonomy.
Before you parition other resources, you'll need to run:

    taxalotlcli build-partition-maps

