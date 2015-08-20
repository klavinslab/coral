'''Retrieve restriction enzymes from rebase.'''
import shutil
import tempfile
import urllib2
import pymbt


class Rebase(object):
    '''Retrieve restriction enzymes from rebase database.'''
    def __init__(self):
        self._tmpdir = None
        self.update()

    def update(self):
        '''Update definitions.'''
        # Download http://rebase.neb.com/rebase/link_withref to tmp
        self._tmpdir = tempfile.mkdtemp()
        try:
            self._rebase_file = self._tmpdir + '/rebase_file'
            print 'Downloading latest enzyme definitions'
            url = 'http://rebase.neb.com/rebase/link_withref'
            header = {'User-Agent': 'Mozilla/5.0'}
            req = urllib2.Request(url, headers=header)
            con = urllib2.urlopen(req)
            with open(self._rebase_file, 'wb') as rebase_file:
                rebase_file.write(con.read())
            # Process into self._enzyme_dict
            self._process_file()
        except urllib2.HTTPError, e:
            print 'HTTP Error: {} {}'.format(e.code, url)
            print 'Falling back on default enzyme list'
            self._enzyme_dict = pymbt.constants.fallback_enzymes
        except urllib2.URLError, e:
            print 'URL Error: {} {}'.format(e.reason, url)
            print 'Falling back on default enzyme list'
            self._enzyme_dict = pymbt.constants.fallback_enzymes
        # Process into RestrictionSite objects? (depends on speed)
        print 'Processing into RestrictionSite instances.'
        self.restriction_sites = {}
        # TODO: make sure all names are unique
        for key, (site, cuts) in self._enzyme_dict.iteritems():
            # Make a site
            try:
                r = pymbt.RestrictionSite(pymbt.DNA(site), cuts, name=key)
                # Add it to dict with name as key
                self.restriction_sites[key] = r
            except ValueError:
                # Encountered ambiguous sequence, have to ignore it until
                # pymbt.DNA can handle ambiguous DNA
                pass

    def get(self, name):
        '''Retrieve enzyme by name.

        :param name: Name of the restriction enzyme, e.g. EcoRV.
        :type name: str
        :returns: Restriction site matching the input name.
        :rtype: pymbt.RestrictionSite
        :raises: Exception when enzyme is not found in the database.

        '''
        # Looks for restriction enzyme by name
        try:
            return self.restriction_sites[name]
        except KeyError:
            raise Exception('Enzyme not found.')

    def _process_file(self):
        '''Process rebase file into dict with name and cut site information.'''
        print 'Processing file'
        with open(self._rebase_file, 'r') as f:
            raw = f.readlines()
        names = [line.strip()[3:] for line in raw if line.startswith('<1>')]
        seqs = [line.strip()[3:] for line in raw if line.startswith('<5>')]
        if len(names) != len(seqs):
            raise Exception('Found different number of enzyme names and '
                            'sequences.')
        self._enzyme_dict = {}
        for name, seq in zip(names, seqs):
            if '?' in seq:
                # Is unknown sequence, don't keep it
                pass
            elif seq.startswith('(') and seq.endswith(')'):
                # Has four+ cut sites, don't keep it
                pass
            elif '^' in seq:
                # Has reasonable internal cut sites, keep it
                top_cut = seq.index('^')
                bottom_cut = len(seq) - top_cut - 1
                site = seq.replace('^', '')
                self._enzyme_dict[name] = (site, (top_cut, bottom_cut))
            elif seq.endswith(')'):
                # Has reasonable external cut sites, keep it
                # (4-cutter also starts with '(')
                # separate site and cut locations
                site, cuts = seq.split('(')
                cuts = cuts.replace(')', '')
                top_cut, bottom_cut = [int(x) + len(site) for x in
                                       cuts.split('/')]
                self._enzyme_dict[name] = (site, (top_cut, bottom_cut))
        shutil.rmtree(self._tmpdir)
