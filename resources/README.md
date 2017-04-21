# Resources subdirectory
This is a small datastore of information about the taxonomy resources.
Eventually, we can put these JSON in JSON-LD syntax.

Currently the idea is that every taxonomy-processing script would know
    that it should be able to look in this directory for a file called
    `resources.json`
Optional files could be added with the convention that their names should
    be `<id>.json` where `<id>` is the identifier used to refer to them
    in the `resources.json` or scripts.
Ideally, the scripts that use these files should be able to read any subset
    and take unions of the base-level attributes.
So, the 