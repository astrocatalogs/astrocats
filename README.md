# Supernova Light Curve Catalog #

This is a repository of supernova light curves, currently maintained by James Guillochon. The catalog is displayed in tabular form on [sne.space](https://sne.space). Individual event data can be downloaded by click the "Source" link on the left and downloading the .dat file associated with a given event; however it's probably easier to use the search function on the table displayed on [sne.space](https://sne.space) to find the event you're looking for.

## Contributing Data ##

If you have a historical archive of data that likely won't be updated, feel free to [e-mail your data directly to me](mailto:jguillochon@cfa.harvard.edu) and we'll do the hard work of converting it to the standard format. However, if you are an active observer and will add data all the time, we prefer to have you [fork the repository](https://bitbucket.org/Guillochon/sne/fork) and update/add events that way. This also allows observers to keep event data "private" in their fork until they are ready to add it to the repository.

## Format of Data Files ##

The data files are always inherently plain text format, although it is preferred to have all datafiles compressed in bzip2 format (.bz2) to minimize disk space. The data is expected to be in key-value form, where a field name is presented, followed by an associated value. For many fields the value is just a single number/string, but the value can itself be a list of key-value pairs, as is the case for photometry entries. Both the key and the value should be enclosed in quotes,  Below is an example datafile for 2002ap:

```
#!text

"name"	"SN2002ap"
"claimedtype"	"Ic"
"host"	"NGC628"
"hvel"	"632.0000"
"redshift"	"0.0021"
"year"	"2002"
"source"	"name"	"<a href='https://www.nhn.ou.edu/~suspect/'>SUSPECT</a>"	"alias"	"1"	"secondary"	"1"
"source"	"name"	"<a href='http://adsabs.harvard.edu/cgi-bin/nph-bib_query?bibcode=astro-ph/0307136&amp;db_key=AST'>astro-ph/0307136</a>"	"alias"	"2"
"source"	"name"	"<a href='http://adsabs.harvard.edu/cgi-bin/nph-bib_query?bibcode=2002MNRAS.332L..73G&amp;db_key=AST'>2002MNRAS.332L..73G</a>"	"alias"	"3"
"source"	"name"	"<a href='http://adsabs.harvard.edu/cgi-bin/nph-bib_query?bibcode=astro-ph/0304010&amp;db_key=AST'>astro-ph/0304010</a>"	"alias"	"4"
"source"	"name"	"<a href='http://adsabs.harvard.edu/cgi-bin/nph-bib_query?bibcode=astro-ph/0209507&amp;db_key=AST'>astro-ph/0209507</a>"	"alias"	"5"
"source"	"name"	"<a href='http://dau.itep.ru/sn/node/72'>Sternberg Astronomical Institute Supernova Light Curve Catalogue</a>"	"alias"	"6"	"secondary"	"1"
"source"	"name"	"VSNET 2003,"	"alias"	"7"
"source"	"name"	"R.J.Foley et al., 2003, PASP, 115, 1220; astro-ph/0307136"	"alias"	"8"
"source"	"name"	"A.Gal-Yam, E.O.Ofek, O.Shemmer, 2002, MNRAS, 332, L73; astro-ph/0204008"	"alias"	"9"
"source"	"name"	"Y.Yoshii et al., 2003, ApJ, 592, 467; astro-ph/0304010"	"alias"	"10"
"source"	"name"	"Unknown,?"	"alias"	"11"
"source"	"name"	"S.B.Pandey et al., 2003, MNRAS, 340, 375; astro-ph/0209507"	"alias"	"12"
"photometry"	"timeunit"	"MJD"	"time"	"52303.04"	"band"	"R"	"abmag"	"15.2"	"source"	"6,7"
"photometry"	"timeunit"	"MJD"	"time"	"52303.2"	"band"	"C"	"abmag"	"14.4"	"source"	"6,7"
"photometry"	"timeunit"	"MJD"	"time"	"52303.39"	"band"	"C"	"abmag"	"14.5"	"source"	"6,7"
```

At the moment, the only mandatory field for a given supernova is its name.