'''Module for gel-related simulations and data structures.'''
try:
    from matplotlib import pyplot
except ImportError:
    # Warning message is handled during library import
    pass


class Band(object):
    '''Represent a gel band - contains length (bp) and optional mass
    properties.

    :param length: Length of the band in base pairs (bp).
    :type length: int
    :param mass: Mass of the band in mg.
    :type mass: float

    '''
    def __init__(self, length, mass=None):
        self.length = length
        self.mass = mass

    def __repr__(self):
        return 'A {} mg, {} bp band.'.format(self.mass, self.length)


class Column(object):
    '''Represent a gel column - contains several Bands and is potentially
    named.

    :param bands: A list of bands.
    :type bands: pymbt.simulation.gel.Band
    :param name: The name of the column.
    :type name: str

    '''
    def __init__(self, band_list, name=None):
        self.bands = band_list
        self.name = name

    def __repr__(self):
        if self.name is not None:
            return '{}: A column of {} bands.'.format(self.name,
                                                      len(self.bands))
        else:
            return 'A column of {} bands.'.format(len(self.bands))


_neb_1kb_lengths = [10002, 8001, 6001, 5001, 4001, 3001, 2000, 1500, 1000, 517,
                    500]
_neb_1kb_masses = [42, 42, 50, 42, 33, 125, 48, 36, 42, 21, 21]
neb_1kb = Column([Band(length, mass) for length, mass in
                  zip(_neb_1kb_lengths, _neb_1kb_masses)])


class Gel(object):
    '''Represent a gel - a set of columns of bands.

    :param

    '''
    def __init__(self, column_list, add_ladder=False,
                 ladder=neb_1kb):
        self.columns = column_list

        if add_ladder:
            column_list.insert(0, ladder)

    def plot(self):
        '''Render gel to plot.'''
        fig = pyplot.figure(figsize=(12, 9), dpi=90)
        sub1 = fig.add_subplot(111)
        max_band = max([band.length for column in self.columns for band in
                        column.bands])
        print max_band
        for i, column in enumerate(self.columns):
            for j, band in enumerate(column.bands):
                if band.mass is not None:
                    alpha = band.mass / 100.0
                    color = 'black'
                else:
                    alpha = 1.0
                    color = 'blue'
                sub1.broken_barh([(i + 0.25, 0.5)],
                                 (band.length, 60), alpha=alpha,
                                 facecolors=color, edgecolors='none')
        sub1.set_xlim([0, len(self.columns)])
        pyplot.show()
