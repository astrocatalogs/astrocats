# Open Supernova Catalog #

[![Build Status](https://travis-ci.org/astrotransients/sne.svg?branch=master)](https://travis-ci.org/astrotransients/sne)

This is a repository of supernova light curves, currently maintained by James Guillochon. The catalog is displayed in tabular form on [sne.space](https://sne.space). Individual event data can be downloaded by click the "Source" link on the left and downloading the .dat file associated with a given event; however it's probably easier to use the search function on the table displayed on [sne.space](https://sne.space) to find the event you're looking for.

## Contributing Data ##

If you have a historical archive of data that likely won't be updated, feel free to [e-mail your data directly to me](mailto:jguillochon@cfa.harvard.edu) and we'll do the hard work of converting it to the standard format. However, if you are an active observer and will add data all the time, we prefer to have you [fork the repository](https://bitbucket.org/Guillochon/sne/fork) and update/add events that way. This also allows observers to keep event data "private" in their fork until they are ready to add it to the repository.

## Format of Data Files ##

The data files are in JSON format. Below is an example datafile for 2002ap:

```json
{
    "SN2002ap": {
        "name": "SN2002ap",
        "aliases": [
            "SN2002ap"
        ],
        "sources": [
            {
                "name": "<a href='https://www.nhn.ou.edu/~suspect/'>SUSPECT</a>",
                "alias": "1",
                "secondary": true
            },
            {
                "name": "<a href='http://adsabs.harvard.edu/cgi-bin/nph-bib_query?bibcode=astro-ph/0307136&amp;db_key=AST'>astro-ph/0307136</a>",
                "alias": "2"
            },
            {
                "name": "<a href='http://adsabs.harvard.edu/cgi-bin/nph-bib_query?bibcode=2002MNRAS.332L..73G&amp;db_key=AST'>2002MNRAS.332L..73G</a>",
                "alias": "3"
            }
        ],
        "photometry": [
            {
                "timeunit": "MJD",
                "time": "52304.169999999925",
                "band": "I",
                "abmag": "13.6990",
                "aberr": "0.0150",
                "source": "1,2"
            },
            {
                "timeunit": "MJD",
                "time": "52305.12999999989",
                "band": "I",
                "abmag": "13.2810",
                "aberr": "0.0180",
                "source": "1,2"
            },
            {
                "timeunit": "MJD",
                "time": "52306.12000000011",
                "band": "I",
                "abmag": "13.0190",
                "aberr": "0.0290",
                "source": "1,2"
            }
        ],
        "discoveryear": 2002,
        "host": "NGC 628",
        "redshift": "0.0021",
        "hvel": "632",
        "claimedtype": "Ic pec",
        "maxday": 15,
        "maxmonth": 2,
        "galra": "013642",
        "galdec": "+154701",
        "snra": "013623.85",
        "sndec": "+154513.2",
        "discoverer": "Hirose",
        "maxyear": 2002,
        "discovermonth": 1,
        "discoverday": 29,
        "z": 0.00211
    }
}
```

At the moment, the only mandatory field for a given supernova is its name.
