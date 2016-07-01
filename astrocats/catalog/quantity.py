"""Class for representing spectra.
"""
from astrocats.catalog.catdict import CatDict
from astrocats.catalog.key import KEY_TYPES, Key, KeyCollection
from astrocats.catalog.utils import is_number


class QUANTITY(KeyCollection):
    # Any
    VALUE = Key('value')
    # Numeric
    ERROR = Key('error', KEY_TYPES.NUMERIC)
    PROB = Key('prob', KEY_TYPES.NUMERIC)
    # Booleans
    UPPER_LIMIT = Key('upperlimit', KEY_TYPES.BOOL)
    # Strings
    UNIT = Key('unit', KEY_TYPES.STRING)
    KIND = Key('kind', KEY_TYPES.STRING)
    SOURCE = Key('source', KEY_TYPES.STRING, compare=False)


class Quantity(CatDict):
    """
    """

    _KEYS = QUANTITY

    def __init__(self, parent, **kwargs):
        self.REQ_KEY_TYPES = [
            [QUANTITY.VALUE],
            [QUANTITY.SOURCE]
        ]
        super().__init__(parent, **kwargs)

        value = self.get(QUANTITY.VALUE, '')
        serror = self.get(QUANTITY.ERROR, '')
        sunit = self.get(QUANTITY.UNIT, '')
        skind = self.get(QUANTITY.KIND, '')
        sprob = self.get(QUANTITY.PROB, '')
        quantity = self._name

        if not self[QUANTITY.VALUE] or value == '--' or value == '-':
            return
        if serror and (not is_number(serror) or float(serror) < 0):
            raise ValueError(parent[parent._KEYS.NAME] + "'s quanta " +
                             quantity +
                             ' error value must be a number and positive.')

        # Set default units
        if not unit and quantity == 'velocity':
            unit = 'KM/s'
        if not unit and quantity == 'ra':
            unit = 'hours'
        if not unit and quantity == 'dec':
            unit = 'degrees'
        if not unit and quantity in ['lumdist', 'comovingdist']:
            unit = 'Mpc'

        # Handle certain quantity
        if quantity == 'alias':
            value = name_clean(value)
            for df in self.get(KEYS.DISTINCT_FROM, []):
                if value == df['value']:
                    return

        if quantity in ['velocity', 'redshift', 'ebv', 'lumdist',
                        'comovingdist']:
            if not is_number(value):
                return
        if quantity == 'host':
            if is_number(value):
                return
            if value.lower() in ['anonymous', 'anon.', 'anon',
                                  'intergalactic']:
                return
            value = host_clean(value)
            if ((not skind and ((value.lower().startswith('abell') and
                                 is_number(value[5:].strip())) or
                                'cluster' in value.lower()))):
                skind = 'cluster'
        elif quantity == KEYS.CLAIMED_TYPE:
            isq = False
            value = value.replace('young', '')
            if value.lower() in ['unknown', 'unk', '?', '-']:
                return
            if '?' in value:
                isq = True
                value = value.strip(' ?')
            for rep in self._source_syns:
                if value in self._source_syns[rep]:
                    value = rep
                    break
            if isq:
                value = value + '?'

        elif quantity in ['ra', 'dec', 'hostra', 'hostdec']:
            (value, sunit) = radec_clean(value, quantity, unit=unit)
        elif quantity == 'maxdate' or quantity == 'discoverdate':
            # Make sure month and day have leading zeroes
            sparts = value.split('/')
            if len(sparts[0]) > 4 and int(sparts[0]) > 0:
                raise ValueError('Date years limited to four digits.')
            if len(sparts) >= 2:
                value = sparts[0] + '/' + sparts[1].zfill(2)
            if len(sparts) == 3:
                value = value + '/' + sparts[2].zfill(2)

            for ii, ct in enumerate(my_quantity_list):
                # Only add dates if they have more information
                if len(ct['value'].split('/')) > len(value.split('/')):
                    return

        if is_number(value):
            value = '%g' % Decimal(value)
        if serror:
            serror = '%g' % Decimal(serror)

        for ii, ct in enumerate(my_quantity_list):
            if ct['value'] == value and sources:
                if 'kind' in ct and skind and ct['kind'] != skind:
                    return
                for source in sources.split(','):
                    if source not in my_quantity_list[ii]['source'].split(','):
                        my_quantity_list[ii]['source'] += ',' + source
                        if serror and 'error' not in my_quantity_list[ii]:
                            my_quantity_list[ii]['error'] = serror
                        if sprob and 'probability' not in my_quantity_list[ii]:
                            my_quantity_list[ii]['probability'] = sprob
                return

        if not sunit:
            sunit = unit
