# Taxalotl tutorial
See [README.md](./README.md) for general documentation.

## Installation and configuration
The following assumes that you are running bash (or some bash-compatible) shell:

    git clone https://github.com/mtholder/taxalotl.git
    cd taxalotl
    bash install.sh
    source env/bin/activate
    export PEYOTL_CONFIG_FILE="$PWD/recommended-peyotl-conf.ini"
    source completion.sh

This should set up a python 3 virtualenv inside an `env` directory
and set up a reasonable configuration.
The configuration file is at `~/.opentreeoflife/taxalotl/taxalotl.conf` and
the scratch directory for operations will be in `~/.opentreeoflife/taxalotl`
When `taxalotl` becomes more stable, you may want to switch 
the `crash_with_stacktraces` setting in `~/.opentreeoflife/taxalotl/taxalotl.conf`
to `false`, but you probably want to leave it set to true now so that you can
submit more useful bug reports.

