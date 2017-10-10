"""
"""


class Biblio:

    def __init__(self):
        self.all_sources = {}


    def add_source_to_entry(self, entry_name, warn=True, **kwargs):
        """Add a source both to a particular entry and to this catalogs `all_sources` dict.
        """
        if entry_name not in self.entries:
            err = "Entry '{}' not in entries!".format(entry_name)
            utils.log_raise(self.log, err, KeyError)

        derive = kwargs.setdefault('derive_parameters', True)
        if not derive and warn:
            self.log.warning("`derive_parameters` should probably be set to True!")

        source = self.entries[entry_name].add_source(**kwargs)
        source_data = self.entries[entry_name][self.proto._KEYS.SOURCES][source]
        print(source)
        print(source_data)
        bibcode = source_data.get(SOURCE.BIBCODE)
        if bibcode is None:
            self.log.warning("no bibcode in source, cannot add to `all_sources`!\n{}".format(
                source_data))
        else:
            self.all_sources[bibcode] = source_data
        return
