# Supernova Light Curve Catalog #

This is a repository of supernova light curves, currently maintained by James Guillochon. The catalog is displayed in tabular form on [sne.space](https://sne.space). Individual event data can be downloaded by click the "Source" link on the left and downloading the .dat file associated with a given event; however it's probably easier to use the search function on the table displayed on [sne.space](https://sne.space) to find the event you're looking for.

## Contributing Data ##

If you have a historical archive of data that likely won't be updated, feel free to [e-mail your data directly to me](mailto:jguillochon@cfa.harvard.edu) and we'll do the hard work of converting it to the standard format. However, if you are an active observer and will add data all the time, we prefer to have you [fork the repository](https://bitbucket.org/Guillochon/sne/fork) and update/add events that way. This also allows observers to keep event data "private" in their fork until they are ready to add it to the repository.

## Format of Data Files ##

The data files are always in plain text format, although in the future we will preferring large datasets to be added in compressed format (probably bzip2). The data is expected to be in key-value form, where a field name is presented, followed by an associated value. For many fields the value is just a single number/string, but the value can itself be a list of key-value pairs, as is the case for photometry entries. Below is an example datafile for 1993J:


```
#!text

"name"	"SN1993J"
"claimedtype"	"IIb"
"discoverer"	"Garcia"
"galdec"	"+690355"
"galra"	"095533"
"host"	"NGC3031"
"hvel"	"-35"
"redshift"	"-0.000117"
"sndec"	"+690113.38"
"snra"	"095524.95"
"photometry"	"MJD"	"49076.0"	"band"	"B"	"instrument"	""	"abmag"	"10.8000"	"aberr"	""	"upperlimit"	"0"
"photometry"	"MJD"	"49076.0"	"band"	"V"	"instrument"	""	"abmag"	"10.7500"	"aberr"	""	"upperlimit"	"0"
"photometry"	"MJD"	"49087.9"	"band"	"V"	"instrument"	""	"abmag"	"11.3900"	"aberr"	""	"upperlimit"	"0"
```