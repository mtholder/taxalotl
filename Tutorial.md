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
The default scratch directory for is `~/.opentreeoflife/taxalotl`, and
the default configuration file location is `~/.opentreeoflife/taxalotl/taxalotl.conf`.
When `taxalotl` becomes more stable, you may want to switch 
the `crash_with_stacktraces` setting in `~/.opentreeoflife/taxalotl/taxalotl.conf`
to `false`, but you probably want to leave it set to true now so that you can
submit more useful bug reports.

## Grokking taxalotl
The script `taxalotl-cli.py` should be on your PATH, this is the front-end to taxalotl.
Sourcing the `completion.sh` script (as was done above) should enable tab-completion of
    commands and options in taxalotl; this makes the front-end **much easier** to work with.


taxalotl uses the [OTifacts](https://github.com/mtholder/OTifacts) repository for information
about the taxonomic resources that Open Tree uses.
taxalotl uses its scratch directory to manage you local copies (and manipulations) of these
resources.
In the context of taxalotl, a "resource" is almost always a taxonomy.

    taxalotl-cli.py status

reports on the state of "your" local copies of the taxonomic resources. If you run this
immediately after installation, your taxalotl resource directory will be empty.
The resources directory is where taxalotl manages notes about what resources could be relevant.
It is a copy of the relevant information in the OTifacts repo.  

If there is any chance that your local resources view lags behind the OTifacts repo on the web,
    you'll want to run:

    taxalotl-cli.py pull-otifacts