# -*- coding: utf-8 -*-

######################################################################
######################################################################
######################################################################
#                                                                    #
# maps.py                                                            #
#                                                                    #
# Tanaus\'u del Pino Alem\'an                                        #
#   Instituto de Astrof\'isica de Canarias                           #
#                                                                    #
######################################################################
######################################################################
#                                                                    #
# Class for map generation and visualization                         #
#                                                                    #
######################################################################
######################################################################
#                                                                    #
######################################################################
######################################################################
#                                                                    #
#  31/10/2018 - V1.1.0 - Changed the whole error handling (TdPA)     #
#                      - Updated help and parameter listing (TdPA)   #
#                      - Added proper management of the wind         #
#                        nodes (TdPA)                                #
#                      - Introduced the shortcuts for a more         #
#                        convenient generation of multiple variables #
#                        of a map (TdPA)                             #
#                      - Added options to generate_map that are to   #
#                        be used by internal routines. May hide      #
#                        them in self.__ variables (TdPA)            #
#                      - Introduced detailed river tracer (TdPA)     #
#                      - Updated saving and loading of maps (TdPA)   #
#                      - Introduced routine for proper refining/     #
#                        zooming of maps (TdPA)                      #
#                      - Fixed height scale, which showed wrong      #
#                        coast, and ensured that scales are kept     #
#                        when using refine (TdPA)                    #
#                                                                    #
#  23/10/2018 - V1.0.0 - First real version. (TdPA)                  #
#                                                                    #
#  19/04/2018 - V0.1.0 - First version. (TdPA)                       #
#                                                                    #
######################################################################
######################################################################
######################################################################

_maps_class__cursor_up_one = '\x1b[1A'
_maps_class__erase_line = '\x1b[2K'
_maps_class__tnormal = '\033[0m'
_maps_class__tnormalgreen = '\033[92m'
_maps_class__tnormalbold = '\033[1m'
_maps_class__tnormalblue = '\033[94m'
_maps_class__tnormalred = '\033[91m'
_maps_class__tnormalseparator = ''
_maps_class__twarning = _maps_class__tnormalred + \
                        '##Warning## ' + \
                        _maps_class__tnormal
_maps_class__terror = _maps_class__tnormalred + \
                      _maps_class__tnormalbold + \
                      '##Error## ' + \
                      _maps_class__tnormal

import sys,copy,struct
try:
    from simplex import *
except ImportError:
    msg = 'Missing simplex.py'
    print _maps_class__terror + msg
    for err in sys.exc_info()[:2]:
        print err
    sys.exit()
except:
    msg = 'Unexpected error importing simplex'
    print _maps_class__terror + msg
    for err in sys.exc_info()[:2]:
        print err
    sys.exit()
try:
    import numpy as np
except ImportError:
    msg = 'Missing numpy'
    print _maps_class__terror + msg
    for err in sys.exc_info()[:2]:
        print err
    sys.exit()
except:
    msg = 'Unexpected error importing numpy'
    print _maps_class__terror + msg
    for err in sys.exc_info()[:2]:
        print err
    sys.exit()
try:
    import random
    _maps_class__random = True
except:
    _maps_class__random = False
try:
    import matplotlib.pyplot as plt
    from matplotlib import cm
    from mpl_toolkits.basemap import Basemap
    from mpl_toolkits.basemap import addcyclic
    from matplotlib.colors import LinearSegmentedColormap
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    import matplotlib.gridspec as gridspec
    _maps_class__2Dsup = True
except:
    _maps_class__2Dsup = False
try:
    from scipy import interpolate
    _maps_class__interp = True
except:
    _maps_class__interp = False
try:
    from scipy.ndimage import filters
    _maps_class__smooth = True
except:
    _maps_class__smooth = False

######################################################################
######################################################################
######################################################################
######################################################################

class maps_class():
    ''' Class with functions to generate/visualize maps
    '''

######################################################################
######################################################################

    def __init__(self, nth=None, nch=None, thrange=None, \
                 chrange=None, seed=None, octaves=None, \
                 frequency=None, persistence=None, water=None, \
                 maxdepth=None, maxheight=None, projection=None, \
                 pthrange=None, pchrange=None, windnnodes=None, \
                 windnodes=None, maxwindspeed=None, \
                 mintemperature=None, maxtemperature=None, name=None):
        ''' Initialize class
        '''

        self.__random = __random
        if not self.__random:
            msg = 'Missing random. No automatic random seeds posible'
            self.__warning(msg)

        self.__can_plot_2D = __2Dsup
        if not self.__can_plot_2D:
            msg = 'Missing pyplot/basemap'
            self.__warning(msg)

        self.__interp = __interp
        if not self.__interp:
            msg = 'Missing scypy/interpolate'
            self.__warning(msg)

        self.__smooth = __smooth
        if not self.__smooth:
            msg = 'Missing scypy/ndimage'
            self.__warning(msg)


        #
        # Parameters
        self.__rade = 180./np.pi
        self.__dera = np.pi/180.

        # Supported maps
        self.__modes = { \
             'cea':'Cylindrical Equal Area', \
             'mbtfpq':'McBryde-Thomas Flat-Polar Quartic', \
             'aeqd':'Azimuthal Equidistant', \
             'sinu':'Sinusoidal', \
             'poly':'Polyconic', \
             'omerc':'Oblique Mercator', \
             'gnom':'Gnomonic', \
             'moll':'Mollweide', \
             'lcc':'Lambert Conformal', \
             'tmerc':'Transverse Mercator', \
             'nplaea':'North-Polar Lambert Azimuthal', \
             'gall':'Gall Stereographic Cylindrical',
             'npaeqd':'North-Polar Azimuthal Equidistant', \
             'mill':'Miller Cylindrical', \
             'merc':'Mercator', \
             'stere':'Stereographic', \
             'eqdc':'Equidistant Conic', \
             'cyl':'Cylindrical Equidistant', \
             'npstere':'North-Polar Stereographic', \
             'spstere':'South-Polar Stereographic', \
             'hammer':'Hammer', \
             'geos':'Geostationary', \
             'nsper':'Near-Sided Perspective', \
             'eck4':'Eckert IV', \
             'aea':'Albers Equal Area', \
             'kav7':'Kavrayskiy VII', \
             'spaeqd':'South-Polar Azimuthal Equidistant', \
             'ortho':'Orthographic', \
             'cass':'Cassini-Soldner', \
             'vandg':'van der Grinten', \
             'laea':'Lambert Azimuthal Equal Area', \
             'splaea':'South-Polar Lambert Azimuthal', \
             'robin':'Robinson'}

        self.__default = {'nth': 360, \
                          'nch': 360, \
                          'thrange': [-90.,90.], \
                          'chrange': [-180.,180.], \
                          'seed': 26894, \
                          'octaves': 20, \
                          'frequency': 0.8, \
                          'persistence': 1.65, \
                          'water': -1, \
                          'maxdepth': 20., \
                          'maxheight': 20., \
                          'windnnodes': -1, \
                          'windnodes': None, \
                          'maxwindspeed': 200., \
                          'mintemperature': -10., \
                          'maxtemperature': 30., \
                          'projection': 'cea'}

        #
        # Control inputs

        # Polar nodes
        if nth is None:
            self.__nth = 360
        else:
            try:
                self.__nth = int(nth)
                if self.__nth < 3:
                    msg = 'nth must be at least 3. Set default'
                    self.__warning(msg)
                    self.__nth = self.__default['nth']
            except ValueError:
                msg = 'ValueError in nth. Set default'
                self.__warning(msg)
                self.__nth = self.__default['nth']
            except:
                msg = 'Unexpected error setting nth'
                error = sys.exc_info()[:2]
                self.__error(msg, error)
                sys.exit()

        # Azimuthal nodes
        if nch is None:
            self.__nch = self.__default['nch']
        else:
            try:
                self.__nch = int(nch)
                if self.__nch < 3:
                    msg = 'nch must be at least 3. Set default'
                    self.__warning(msg)
                    self.__nch = self.__default['nch']
            except ValueError:
                msg = 'ValueError in nch'
                self.__warning(msg)
                self.__nch = self.__default['nch']
            except:
                msg = 'Unexpected error setting nch'
                error = sys.exc_info()[:2]
                self.__error(msg, error)
                sys.exit()

        # Polar range
        if thrange is None:
            self.__thrange = self.__default['thrange']
        else:
            if len(thrange) == 2:
                try:
                    self.__thrange = []
                    for th in thrange:
                        self.__thrange.append(float(th))
                    if self.__thrange[0] >= self.__thrange[1] or \
                       self.__thrange[0] < -90. or \
                       self.__thrange[0] > 90. or \
                       self.__thrange[1] < -90. or \
                       self.__thrange[1] > 90.:
                        msg = 'Bad range in thrange. Set default'
                        self.__warning(msg)
                        self.__thrange = self.__default['thrange']
                except ValueError:
                    msg = 'ValueError in thrange. Set default'
                    self.__warning(msg)
                    self.__thrange = self.__default['thrange']
                except:
                    msg = 'Unexpected error setting thrange'
                    error = sys.exc_info()[:2]
                    self.__error(msg, error)
                    sys.exit()
            else:
                msg = 'Wrong format thrange. Set default'
                self.__warning(msg)
                self.__thrange = self.__default['thrange']

        # Polar range for plot
        if thrange is None:
            self.__pthrange = copy.deepcopy(self.__thrange)
        else:
            if len(pthrange) == 2:
                try:
                    self.__pthrange = []
                    for th in pthrange:
                        self.__pthrange.append(float(th))
                    if self.__pthrange[0] >= self.__pthrange[1] or \
                       self.__pthrange[0] < -90. or \
                       self.__pthrange[0] > 90. or \
                       self.__pthrange[1] < -90. or \
                       self.__pthrange[1] > 90.:
                        msg = 'Bad range in pthrange. Set default'
                        self.__warning(msg)
                        self.__pthrange = \
                                         copy.deepcopy(self.__thrange)
                except ValueError:
                    msg = 'ValueError in pthrange. Set default'
                    self.__warning(msg)
                    self.__pthrange = copy.deepcopy(self.__thrange)
                except:
                    msg = 'Unexpected error checking pthrange'
                    error = sys.exc_info()[:2]
                    self.__error(msg, error)
                    sys.exit()
            else:
                msg = 'Wrong format pthrange. Set default'
                self.__warning(msg)
                self.__pthrange = copy.deepcopy(self.__thrange)

        # Azimuthal range
        if chrange is None:
            self.__chrange = self.__default['chrange']
        else:
            if len(chrange) == 2:
                try:
                    self.__chrange = []
                    for ch in chrange:
                        self.__chrange.append(float(ch))
                    if self.__chrange[0] >= self.__chrange[1] or \
                       self.__chrange[0] < -180. or \
                       self.__chrange[0] > 180. or \
                       self.__chrange[1] < -180. or \
                       self.__chrange[1] > 180.:
                        msg = 'Bad range in chrange. Set default'
                        self.__warning(msg)
                        self.__chrange = self.__default['chrange']
                except ValueError:
                    msg = 'ValueError in chrange. Set default'
                    self.__warning(msg)
                    self.__chrange = self.__default['chrange']
                except:
                    msg = 'Unexpected error setting chrange'
                    error = sys.exc_info()[:2]
                    self.__error(msg, error)
                    sys.exit()
            else:
                msg = 'Wrong format chrange. Set default'
                self.__warning(msg)
                self.__chrange = self.__default['chrange']

        # Azimuthal range for plot
        if chrange is None:
            self.__pchrange = copy.deepcopy(self.__chrange)
        else:
            if len(pchrange) == 2:
                try:
                    self.__pchrange = []
                    for ch in pchrange:
                        self.__pchrange.append(float(ch))
                    if self.__pchrange[0] >= self.__pchrange[1] or \
                       self.__pchrange[0] < -180. or \
                       self.__pchrange[0] > 180. or \
                       self.__pchrange[1] < -180. or \
                       self.__pchrange[1] > 180.:
                        msg = 'Bad range in pchrange. Set default'
                        self.__warning(msg)
                        self.__pchrange = \
                                         copy.deepcopy(self.__chrange)
                except ValueError:
                    msg = 'ValueError in pchrange'
                    self.__warning(msg)
                    self.__pchrange = copy.deepcopy(self.__chrange)
                except:
                    msg = 'Unexpected error setting pchrange'
                    error = sys.exc_info()[:2]
                    self.__error(msg, error)
                    sys.exit()
            else:
                msg = 'Wrong format pchrange'
                self.__warning(msg)
                self.__pchrange = copy.deepcopy(self.__chrange)

        # Seed
        if seed is None:
            if self.__random:
                self.__seed = random.randint(1,100000)
            else:
                self.__seed = self.__default['seed']
        else:
            try:
                self.__seed = int(seed)
            except ValueError:
                msg = 'ValueError in seed. Set default'
                self.__warning(msg)
                if self.__random:
                    self.__seed = random.randint(1,100000)
                else:
                    self.__seed = self.__default['seed']
            except:
                msg = 'Unexpected error setting seed'
                error = sys.exc_info()[:2]
                self.__error(msg, error)
                sys.exit()

        # Octaves
        if octaves is None:
            self.__octaves = self.__default['octaves']
        else:
            try:
                self.__octaves = int(octaves)
                if self.__optaves < 1:
                    msg = 'Octaves must be positive. Set default'
                    self.__warning(msg)
                    self.__octaves = self.__default['octaves']
            except ValueError:
                msg = 'ValueError in octaves. Set default'
                self.__warning(msg)
                self.__octaves = self.__default['octaves']
            except:
                msg = 'Unexpected error setting octaves'
                error = sys.exc_info()[:2]
                self.__error(msg, error)
                sys.exit()

        # Frequency
        if frequency is None:
            self.__frequency = self.__default['frequency']
        else:
            try:
                self.__frequency = float(frequency)
                if self.__frequency < 0. or self.__frequency > 1.:
                    msg = 'Frequency must be between ' + \
                          ' 0 and 1. Set default'
                    self.__warning(msg)
                    self.__frequency = self.__default['frequency']
            except ValueError:
                msg = 'ValueError in frequency. Set default'
                self.__warning(msg)
                self.__frequency = self.__default['frequency']
            except:
                msg = 'Unexpected error setting frequency'
                error = sys.exc_info()[:2]
                self.__error(msg, error)
                sys.exit()

        # Persistence
        if persistence is None:
            self.__persistence = self.__default['persistence']
        else:
            try:
                self.__persistence = int(persistence)
                if self.__persistence < 0.:
                    msg = 'Persistence must be positive. Set default'
                    self.__warning(msg)
                    self.__persistence = self.__default['persistence']
            except ValueError:
                msg = 'ValueError in persistence. Set default'
                self.__warning(msg)
                self.__persistence = self.__default['persistence']
            except:
                msg = 'Unexpected error setting persistence'
                error = sys.exc_info()[:2]
                self.__error(msg, error)
                sys.exit()

        # Water
        if water is None:
            self.__water = self.__default['water']
        else:
            try:
                self.__water = float(water)
                if self.__water < 0. or self.__water > 100.:
                    msg = 'Water must be between ' + \
                          '0 and 100. Set default'
                    self.__warning(msg)
                    self.__water = self.__default['water']
            except ValueError:
                msg = 'ValueError in water. Set default'
                self.__warning(msg)
                self.__water = self.__default['water']
            except:
                msg = 'Unexpected error setting water'
                error = sys.exc_info()[:2]
                self.__error(msg, error)
                sys.exit()

        # Maxdepth
        if maxdepth is None:
            self.__maxdepth = self.__default['maxdepth']
        else:
            try:
                self.__maxdepth = float(maxdepth)
                if self.__maxdepth <= 0.:
                    msg = 'Maxdepth must be ' + \
                          'positive. Set default'
                    self.__warning(msg)
                    self.__maxdepth = self.__default['maxdepth']
            except ValueError:
                msg = 'ValueError in maxdepth. Set default'
                self.__warning(msg)
                self.__maxdepth = self.__default['maxdepth']
            except:
                msg = 'Unexpected error setting maxdepth'
                error = sys.exc_info()[:2]
                self.__error(msg, error)
                sys.exit()

        # Maxheight
        if maxheight is None:
            self.__maxheight = self.__default['maxheight']
        else:
            try:
                self.__maxheight = float(maxheight)
                if self.__maxheight <= 0.:
                    msg = 'Maxheight must be ' + \
                          'positive. Set default'
                    self.__maxheight = self.__default['maxheight']
            except ValueError:
                msg = 'ValueError in maxheight. Set default'
                self.__maxheight = self.__default['maxheight']
            except:
                msg = 'Unexpected error setting maxheight'
                error = sys.exc_info()[:2]
                self.__error(msg, error)
                sys.exit()

        # Windnnodes
        if windnnodes is None:
            self.__wind_nnodes = self.__default['windnnodes']
        else:
            try:
                self.__wind_nnodes = int(windnnodes)
            except ValueError:
                msg = 'ValueError in maxheight. Set default'
                self.__wind_nnodes = self.__default['windnnodes']
            except:
                msg = 'Unexpected error setting windnnodes'
                error = sys.exc_info()[:2]
                self.__error(msg, error)
                sys.exit()

        # Windnodes
        if windnodes is None:
            self.__wind_nodes = self.__default['windnodes']
        else:
            try:
                if isinstance(windnodes, list):
                    self.__wind_nodes = []
                    for nod in windnodes:
                        if self.__wind_nodes is None:
                            break
                        if isinstance(nod, list):
                            if len(nod) == 4:
                                if nod[0] < -90. or nod[0] > 90.:
                                    msg = 'Latitude of node must ' + \
                                          'be between -90 and 90.'
                                    self.__warning(msg)
                                    self.__wind_nodes = \
                                           self.__default['windnodes']
                                    continue
                                if nod[1] < -180. or nod[1] > 180.:
                                    msg = 'Longitude of node ' + \
                                           'must be between ' + \
                                           '-180. and 180.'
                                    self.__warning(msg)
                                    self.__wind_nodes = \
                                           self.__default['windnodes']
                                    continue
                                if int(np.absolute(nod[2])) != 1:
                                    msg = 'Sign of node must be ' + \
                                          '-1. or 1.'
                                    self.__warning(msg)
                                    self.__wind_nodes = \
                                           self.__default['windnodes']
                                    continue
                                if nod[3] < 0.:
                                    msg = 'Weight of node must ' + \
                                          'be positive'
                                    self.__warning(msg)
                                    self.__wind_nodes = \
                                           self.__default['windnodes']
                                    continue
                                self.__wind_nodes.append(nod)
                            else:
                                msg = 'Windnodes must be a list ' + \
                                      'of 4 element lists'
                                self.__warning(msg)
                                self.__wind_nodes = \
                                           self.__default['windnodes']
                                continue
                        else:
                            msg = 'Windnodes must be a list of 4 ' + \
                                  'element lists'
                            self.__warning(msg)
                            return
                else:
                    msg = 'Windnodes must be a list of 4 element ' + \
                          'lists'
                    self.__warning(msg)
                    self.__wind_nodes = self.__default['windnodes']
            except:
                msg = 'Unexpected error setting windnodes'
                error = sys.exc_info()[:2]
                self.__error(msg, error)
                sys.exit()

        # Maxwindspeed
        if maxwindspeed is None:
            self.__maxwindspeed = self.__default['maxwindspeed']
        else:
            try:
                self.__maxwindspeed = float(maxwindspeed)
                if self.__maxwindspeed <= 0.:
                    msg = 'Maxwindspeed must be ' + \
                          'positive. Set default'
                    self.__warning(msg)
                    self.__maxwindspeed = \
                                        self.__default['maxwindspeed']
            except ValueError:
                msg = 'ValueError in maxwindspeed. ' + \
                      'Set default'
                self.__warning(msg)
                self.__maxwindspeed = self.__default['maxwindspeed']
            except:
                msg = 'Unexpected error setting maxwindspeed'
                error = sys.exc_info()[:2]
                self.__error(msg, error)
                sys.exit()

        # Mintemperature
        if mintemperature is None:
            self.__mintemperature = self.__default['mintemperature']
        else:
            try:
                self.__mintemperature = float(mintemperature)
                if self.__mintemperature >= 0.:
                    msg = 'mintemperature must be ' + \
                          'negative. Set default'
                    self.__warning(msg)
                    self.__mintemperature = \
                                      self.__default['mintemperature']
            except ValueError:
                msg = 'ValueError in mintemperature. ' + \
                      'Set default'
                self.__warning(msg)
                self.__mintemperature = \
                                      self.__default['mintemperature']
            except:
                msg = 'Unexpected error setting mintemperature'
                error = sys.exc_info()[:2]
                self.__error(msg, error)
                sys.exit()

        # Maxtemperature
        if maxtemperature is None:
            self.__maxtemperature = self.__default['maxtemperature']
        else:
            try:
                self.__maxtemperature = float(maxtemperature)
                if self.__maxtemperature <= 0.:
                    msg = 'Maxtemperature must be ' + \
                          'positive. Set default'
                    self.__warning(msg)
                    self.__maxtemperature = \
                                      self.__default['maxtemperature']
            except ValueError:
                msg = 'ValueError in maxtemperature. ' + \
                      'Set default'
                self.__warning(msg)
                self.__maxtemperature = \
                                      self.__default['maxtemperature']
            except:
                msg = 'Unexpected error setting maxtemperature'
                error = sys.exc_info()[:2]
                self.__error(msg, error)
                sys.exit()

        # Projection
        if projection is None:
            self.__projection = self.__default['projection']
        else:
            if projection in self.__modes.keys():
                self.__projection = projection
            else:
                msg = 'Unknown projection. Set default'
                self.__warning(msg)
                self.__projection = self.__default['projection']

        # Name
        if name is None:
            self.__name = 'map'
        else:
            self.__name = str(name)

        # Generated
        self.__exist = False
        self.__wind_exist = False
        self.__temperature_exist = False
        self.__moist_exist = False
        self.__biome_exist = False
        self.__river_exist = False
        self.__lake_exist = False

######################################################################
######################################################################

    def __print(self, msg=None):
        ''' Normal message
        '''

        # Argument handling
        if msg is None:
            msg = ''

        # Set message
        msg = __tnormal + msg

        # Print error
        print(msg)

######################################################################
######################################################################

    def __error(self, msg=None, error=None):
        ''' Error message
        '''

        # Argument handling
        if msg is None:
            msg = ''

        # Set message
        msg = __terror + msg

        # Print error
        print(msg)
        if error is not None:
            for err in error:
                print('  '+str(err))


######################################################################
######################################################################

    def __warning(self, msg=None, error=None):
        ''' Warning message
        '''

        # Argument handling
        if msg is None:
            msg = ''

        # Set message
        msg = __twarning + msg

        # Print error
        print(msg)
        if error is not None:
            for err in error:
                print('  '+str(err))

######################################################################
######################################################################

    def help(self):
        ''' Shows a short guide of the class
        '''

        msg = ' Welcome to my map generator.\n\n'

        msg += ' If you are reading this, you already figured ' + \
               "out how\n to call this class, so let's skip " + \
               'that part.\n\n'

        msg += ' To see a list of the parameters of the map ' + \
               'generation\n use the method maps_class.' + \
               'parameter_list().\n The list shows also the ' + \
               'methods to change each\n parameter and the valid ' + \
               'limits\n The same methods, changing "set" to ' + \
               '"get", can be\n used to recall the value of an ' + \
               'individual parameter\n The method maps_class.pa' + \
               'rameter_list() also shows\n the current values.\n\n'

        msg += ' To store the map parameters, use\n maps_class.' + \
               'save_parameters(). The default name is\n the ' + \
               'parameter "name" plus ".par"\n To load map ' + \
               'parameters, use\n maps_class.load_parameters(). ' + \
               'The default name is\n the parameter "name" plus ' + \
               '".par"\n\n'

        msg += ' To generate a map, use the method\n maps_class.' + \
               'generate_map()\n\n'

        msg += ' To generate a wind map, use the method\n ' + \
               'maps_class.generate_wind()\n\n'

        msg += ' To generate a temperature map, use the method\n ' + \
               'maps_class.generate_temperature()\n\n'

        msg += ' To generate a moisture map, use the method\n ' + \
               'maps_class.generate_moist()\n\n'

        msg += ' To generate a biome map, use the method\n ' + \
               'maps_class.generate_biome()\n\n'

        msg += ' To generate rivers and lakes, use the method\n ' + \
               'maps_class.generate_rivers()\n'

        msg += ' You can use the argument "detailed=True" for a ' + \
               'more precise, albeit slower, tracing of rivers.\n\n'

        msg += ' Each new element that can be generated ' + \
               'requires\n every previous element in the order ' + \
               'described above\n\n'

        msg += ' You can use maps_class.generate_weather() to ' + \
               'generate\n wind, temperature, moisture and ' + \
               'biomes once heights\n are done\n\n'

        msg += ' You can use maps_class.generate_all() to ' + \
               'generate\n every aspect of the map in one go. ' + \
               'The mode for the rivers is "defailed=False" ' + \
               'with this shortcut\n\n'

        msg += ' Everything except the heights requires a full ' + \
               'globe map\n\n'

        msg += ' The method maps_class.refine() can be used to' + \
               ' get a \nsubsection of a map with a much better' + \
               ' resolution.\n If the weather and rivers are' + \
               ' already calculated, they will\n be also refined,' + \
               ' if not, you will not be able to generate them\n' + \
               ' afterwards. Be careful, this process is not' + \
               ' reversible, \nit is recommended to save the' + \
               ' global map before refining.\n\n'

        msg += ' To visualize a map in 2D, use the method\n maps' + \
               '_class.draw(). You can see a list with the\n' + \
               ' available projections using the method\n maps_cl' + \
               'ass.help_projection(). To select a\n projection, ' + \
               'set the projection parameter to the\n desired ' + \
               'short keyword using the set_projection(val)\n' + \
               ' method.\n For some projections you may want to ' + \
               'rotate the map\n in the longitude axis. Use self.' + \
               'rotate(val) to\n rotate the longitude by val\n\n'

        msg += ' To store a generated map, use maps_class.save_' + \
               'map().\n The default name is the parameter ' + \
               '"name" plus ".map".\n To load a map, use maps_' + \
               'class.load_map(). The default\n name is the ' + \
               'parameter "name" plus ".map".\n\n'

        msg += ' To store the map in vtk format, use\n maps_class' + \
               '.save_vtk(). The default name is the\n parameter' + \
               ' "name" plus ".vtk"\n\n'
        
        msg += ' Enjoy the maps!'

        self.__print(msg)


######################################################################
######################################################################

    def parameter_list(self):
        ''' Shows a list of the parameters and current values
        '''

        msg =  '*************************************************' + \
               '*************************************************' + \
               '*************************************************' + \
               '***\n'
        msg += 'nth, number of latitude nodes. Change with ' + \
               '.set_nth(integer>3)\n'
        msg += '  nth = '+str(self.__nth)+'\n'
        msg += 'nch, number of longitude nodes. Change with ' + \
               '.set_nch(integer>3)\n'
        msg += '  nch = '+str(self.__nch)+'\n'
        msg += 'thrange, minimum and maximum latitudes (deg). ' + \
               'Change with .set_thrange(float,float), with the ' + \
               'arguments the minimum and the maximum, ' + \
               'respectively. -90 <= float <= 90\n'
        msg += '  thrange = '+str(self.__thrange[0])+':'+ \
                              str(self.__thrange[1])+'\n'
        msg += 'chrange, minimum and maximum longitudes (deg). ' + \
               'Change  with .set_chrange(float,float), with the ' + \
               'arguments the minimum and the maximum, ' + \
               'respectively. -180 <= float <= 180\n'
        msg += '  chrange = '+str(self.__chrange[0])+':'+ \
                              str(self.__chrange[1])+'\n'
        msg += 'pthrange, minimum and maximum latitudes for ' + \
               'plotting (deg). Change with ' + \
               '.set_pthrange(float,float), with the arguments ' + \
               'the minimum and the maximum, respectively. -90' + \
               ' <= float <= 90\n'
        msg += '  pthrange = '+str(self.__pthrange[0])+':'+ \
                               str(self.__pthrange[1])+'\n'
        msg += 'pchrange, minimum and maximum longitudes for ' + \
               'plotting (deg). Change with ' + \
               '.set_pchrange(float,float), with the arguments ' + \
               'the minimum and the maximum, respectively. -180' + \
               ' <= float <= 180\n'
        msg += '  pchrange = '+str(self.__pchrange[0])+':'+ \
                               str(self.__pchrange[1])+'\n'
        msg += 'seed, it defines the map geography. Change with ' + \
               '.set_seed(integer) or .gen_seed() for a random ' + \
               'value\n'
        msg += '  seed = '+str(self.__seed)+'\n'
        msg += 'octaves, number of noise layers to add. Change ' + \
               'with .set_octaves(0<=integer)\n'
        msg += '  octaves = '+str(self.__octaves)+'\n'
        msg += 'frequency, how much the dimensions are expanded ' + \
               'in each octave. Change with ' + \
               '.set_frequency(0<=float<=1)\n'
        msg += '  frequency = '+str(self.__frequency)+'\n'
        msg += 'persistence, factor that each layer is weaker ' + \
               'than the last one. Change with ' + \
               '.set_persistence(float>=0)\n'
        msg += '  persistence = '+str(self.__persistence)+'\n'
        msg += 'water, relative quantity of points with water in ' + \
               '%. -1 means automatic. Change with ' + \
               '.set_water(float<=100)\n'
        msg += '  water = '+str(self.__water)+'\n'
        msg += 'maxdepth, maximum depth in km. Change with ' + \
               '.set_maxdepth(float>0)\n'
        msg += '  maxdepth = '+str(self.__maxdepth)+'\n'
        msg += 'maxheight, maximum height in km. Change with ' + \
               '.set_maxheight(float>0)\n'
        msg += '  maxheight = '+str(self.__maxheight)+'\n'
        msg += 'windnnodes, number of nodes in wind generator. ' + \
               'Change with .set_windnnodes(integer)\n'
        msg += '  windnnodes = '+str(self.__wind_nnodes)+'\n'
        msg += 'windnodes, precise nodes for the wind generator. ' + \
               'has priotiry over windnnodes. Change with .set_' + \
               'windnodes([[latitude, longitude, sign, weight], ' + \
               '...])\n'
        if self.__wind_nodes is None:
            msg += '  windnodes = None\n'
        else:
            msg += '  windnodes = ',self.__wind_nodes[0],'\n'
            for inod in range(1,len(self.__wind_nodes)):
                msg += '              '+ \
                        str(self.__wind_nodes[inod][0]) + ', ' + \
                        str(self.__wind_nodes[inod][1]) + ', ' + \
                        str(self.__wind_nodes[inod][2]) + ', ' + \
                        str(self.__wind_nodes[inod][3]) + '\n'
        msg += 'maxwindspeed, maximum speed in km/h. Change with ' + \
               '.set_maxwindspeed(float>0)\n'
        msg += '  maxwindspeed = '+str(self.__maxwindspeed)+'\n'
        msg += 'mintemperature, minimum temperature in C. Change ' + \
               'with .set_mintemperature(float<0)\n'
        msg += '  mintemperature = '+str(self.__mintemperature)+'\n'
        msg += 'maxtemperature, maximum temperature in C. Change ' + \
               'with .set_maxtemperature(float>0)\n'
        msg += '  mintemperature = '+str(self.__mintemperature)+'\n'
        msg += 'projection, type of projection to plot in 2D. ' + \
               'Change with .set_projection(string). Check ' + \
               '.help_projection() for possible strings\n'
        msg += '  projection = '+self.__modes[self.__projection]+'\n'
        msg += 'name, path to the file to store the map. Change ' + \
               'with .set_name(string)\n'
        msg += '  name = '+self.__name+'\n'
        msg += '*************************************************' + \
               '*************************************************' + \
               '*************************************************' + \
               '***\n'

        self.__print(msg)

######################################################################
######################################################################

    def set_nth(self,nth):
        ''' Changes value of polar nodes
        '''

        try:

            nth = int(nth)

            if nth < 3:
                msg = 'nth must be at least 3'
                self.__warning(msg)
            else:
                self.__nth = nth

        except ValueError:
            msg = 'nth must be integer'
            self.__error(msg)
        except:
            msg = 'Unexpected error'
            error = sys.exc_info()[:2]
            self.__error(msg, error)

######################################################################
######################################################################

    def get_nth(self):
        ''' Get value of polar nodes
        '''

        return self.__nth

######################################################################
######################################################################

    def set_nch(self,nch):
        ''' Changes value of azimuthal nodes
        '''

        try:

            nch = int(nch)

            if nch < 3:
                msg = 'nch must be at least 3'
                self.__warning(msg)
            else:
                self.__nch = nch

        except ValueError:
            msg = 'nch must be integer'
            self.__error(msg)
        except:
            msg = 'Unexpected error'
            error = sys.exc_info()[:2]
            self.__error(msg, error)

######################################################################
######################################################################

    def get_nch(self):
        ''' Get value of azimuthal nodes
        '''

        return self.__nch

######################################################################
######################################################################

    def set_thrange(self,th0,th1):
        ''' Set value of polar range, specify bottom and upper limits
        '''

        try:
            th0 = float(th0)
            th1 = float(th1)
            if th0 >= th1 or \
               th0 < -90. or \
               th0 > 90. or \
               th1 < -90. or \
               th1 > 90.:
                msg = 'Bad range.'
                self.__warning(msg)
            else:
                self.__thrange = [th0,th1]
        except ValueError:
            msg = 'ValueError in thrange'
            self.__error(msg)
        except:
            msg = 'Unexpected error'
            error = sys.exc_info()[:2]
            self.__error(msg, error)

######################################################################
######################################################################

    def get_thrange(self):
        ''' Get value of polar range
        '''

        return self.__thrange

######################################################################
######################################################################

    def set_pthrange(self,th0,th1):
        ''' Set value of polar range for plotting, specify bottom and
            upper limits
        '''

        try:
            th0 = float(th0)
            th1 = float(th1)
            if th0 >= th1 or \
               th0 < -90. or \
               th0 > 90. or \
               th1 < -90. or \
               th1 > 90.:
                msg = '#Warning# Bad range.'
                self.__warning(msg)
            else:
                self.__pthrange = [th0,th1]
        except ValueError:
            msg = 'ValueError in pthrange'
            self.__error(msg)
        except:
            msg = 'Unexpected error'
            error = sys.exc_info()[:2]
            self.__error(msg, error)

######################################################################
######################################################################

    def get_pthrange(self):
        ''' Get value of polar range for plotting
        '''

        return self.__pthrange

######################################################################
######################################################################

    def set_chrange(self,ch0,ch1):
        ''' Set value of azimuthal range, specify bottom and upper
            limits
        '''

        try:
            ch0 = float(ch0)
            ch1 = float(ch1)
            if ch0 >= ch1 or \
               ch0 < -180. or \
               ch0 > 180. or \
               ch1 < -180. or \
               ch1 > 180.:
                msg = 'Bad range.'
                self.__warning(msg)
            else:
                self.__chrange = [ch0,ch1]
        except ValueError:
            msg = 'ValueError in chrange'
            self.__error(msg)
        except:
            msg = 'Unexpected error'
            error = sys.exc_info()[:2]
            self.__error(msg, error)

######################################################################
######################################################################

    def get_chrange(self):
        ''' Get value of azimuthal range
        '''

        return self.__chrange

######################################################################
######################################################################

    def set_pchrange(self,ch0,ch1):
        ''' Set value of azimuthal range for plotting, specify bottom
            and upper limits
        '''

        try:
            ch0 = float(ch0)
            ch1 = float(ch1)
            if ch0 >= ch1 or \
               ch0 < -180. or \
               ch0 > 180. or \
               ch1 < -180. or \
               ch1 > 180.:
                msg = 'Bad range.'
                self.__warning(msg)
            else:
                self.__pchrange = [ch0,ch1]
        except ValueError:
            msg = 'ValueError in pchrange'
            self.__error(msg)
        except:
            msg = 'Unexpected error'
            error = sys.exc_info()[:2]
            self.__error(msg, error)

######################################################################
######################################################################

    def get_pchrange(self):
        ''' Get value of azimuthal range for plotting
        '''

        return self.__pchrange

######################################################################
######################################################################

    def set_seed(self,seed):
        ''' Changes value of the seed
        '''

        self.__seed = seed

######################################################################
######################################################################

    def gen_seed(self):
        ''' Generates a random seed
        '''

        if self.__random:
            self.__seed = random.randint(1,100000)
        else:
            msg = 'No support of random generator'
            self.__error(msg)

######################################################################
######################################################################

    def get_seed(self):
        ''' Get value of seed
        '''

        return self.__seed

######################################################################
######################################################################

    def set_octaves(self,octav):
        ''' Set value of octaves
        '''

        try:
            octav = int(octav)
            if optav < 1:
                msg = 'octaves must be positive.'
                self.__warning(msg)
            else:
                self.__octaves = octav
        except ValueError:
            msg = 'ValueError in octaves'
            self.__error(msg)
        except:
            msg = 'Unexpected error'
            error = sys.exc_info()[:2]
            self.__error(msg, error)

######################################################################
######################################################################

    def get_octaves(self):
        ''' Get value of octaves
        '''

        return self.__octaves

######################################################################
######################################################################

    def set_frequency(self,freq):
        ''' Set value of frequency
        '''

        try:
            freq = float(freq)
            if freq < 0. or freq > 1.:
                msg = 'frequency must be between ' + \
                      ' 0 and 1.'
                self.__warning(msg)
            else:
                self.__frequency = freq
        except ValueError:
            msg = '##Error## ValueError in octaves'
            self.__error(msg)
        except:
            msg = 'Unexpected error'
            error = sys.exc_info()[:2]
            self.__error(msg, error)

######################################################################
######################################################################

    def get_frequency(self):
        ''' Get value of frequency
        '''

        return self.__frequency

######################################################################
######################################################################

    def set_persistence(self,pers):
        ''' Set value of persistence
        '''

        try:
            pers = float(pers)
            if pers < 0.:
                msg = 'persistence must be ' + \
                      'positive.'
                self.__warning(msg)
            else:
                self.__persistence = pers
        except ValueError:
            msg = 'ValueError in persistence'
            self.__error(msg)
        except:
            msg = 'Unexpected error'
            error = sys.exc_info()[:2]
            self.__error(msg, error)

######################################################################
######################################################################

    def get_persistence(self):
        ''' Get value of persistence
        '''

        return self.__persistence

######################################################################
######################################################################

    def set_water(self,wat):
        ''' Set value of water. Negative means automatic
        '''

        try:
            wat = float(wat)
            if wat >= 0.:
                if wat > 100.:
                    msg = 'Maximum water is 100.'
                    self.__warning(msg)
                else:
                    self.__water = wat
            else:
                self.__water = -1
        except ValueError:
            msg = 'ValueError in water'
            self.__error(msg)
        except:
            msg = 'Unexpected error'
            error = sys.exc_info()[:2]
            self.__error(msg, error)


######################################################################
######################################################################

    def get_water(self):
        ''' Get value of water
        '''

        return self.__water

######################################################################
######################################################################

    def set_maxdepth(self,depth):
        ''' Set value of maxdepth
        '''

        try:
            depth = float(depth)
            if depth <= 0.:
                msg = 'maxdepth must be ' + \
                      'positive.'
                self.__warning(msg)
            else:
                self.__maxdepth = depth
        except ValueError:
            msg = 'ValueError in maxdepth'
            self.__error(msg)
        except:
            msg = 'Unexpected error'
            error = sys.exc_info()[:2]
            self.__error(msg, error)

######################################################################
######################################################################

    def get_maxdepth(self):
        ''' Get value of maxdepth
        '''

        return self.__maxdepth

######################################################################
######################################################################

    def set_maxheight(self,height):
        ''' Set value of maxheight
        '''

        try:
            height = float(height)
            if height <= 0.:
                msg = 'maxheight must be ' + \
                      'positive.'
                self.__warning(msg)
            else:
                self.__maxheight = height
        except ValueError:
            msg = 'ValueError in maxheight'
            self.__error(msg)
        except:
            msg = 'Unexpected error'
            error = sys.exc_info()[:2]
            self.__error(msg, error)

######################################################################
######################################################################

    def get_maxheight(self):
        ''' Get value of maxheight
        '''

        return self.__maxheight

######################################################################
######################################################################

    def set_windnnodes(self, nodes):
        ''' Set value of windnnodes
        '''

        try:
            self.__windnnodes = int(nodes)
        except:
            msg = 'Unexpected error'
            error = sys.exc_info()[:2]
            self.__error(msg, error)
            return

######################################################################
######################################################################

    def get_winddnnodes(self):
        ''' Get value of windnnodes
        '''

        return self.__windnnodes

######################################################################
######################################################################

    def set_windnodes(self, nodes):
        ''' Set wind nodes
        '''

        try:
            if isinstance(nodes, list):
                windnodes = []
                for nod in nodes:
                    if isinstance(nod, list):
                        if len(nod) == 4:
                            if nod[0] < -90. or nod[0] > 90.:
                                msg = 'Latitude of node must ' + \
                                      'be between -90 and 90.'
                                self.__warning(msg)
                                return
                            if nod[1] < -180. or nod[1] > 180.:
                                msg = 'Longitude of node must ' + \
                                      'be between -180. and 180.'
                                self.__warning(msg)
                                return
                            if int(np.absolute(nod[2])) != 1:
                                msg = 'Sign of node must be ' + \
                                      '-1. or 1.'
                                self.__warning(msg)
                                return
                            if nod[3] < 0.:
                                msg = 'Weight of node must be ' + \
                                      'positive'
                                self.__warning(msg)
                                return
                            windnodes.append(nod)
                        else:
                            msg = 'Windnodes must be a list of 4 ' + \
                                  'element lists'
                            self.__warning(msg)
                            return
                    else:
                        msg = 'Windnodes must be a list of 4 ' + \
                              'element lists'
                        self.__warning(msg)
                        return
                self.__wind_nodes = copy.deepcopy(windnodes)
            else:
                msg = 'Windnodes must be a list of 4 element lists'
                self.__warning(msg)
                return
        except:
            msg = 'Unexpected error'
            error = sys.exc_info()[:2]
            self.__error(msg, error)
            return

######################################################################
######################################################################

    def generate_windnodes(self):
        ''' Random generation of wind nodes
        '''

        # If there is random number generator
        if self.__random:

            # Initialize seed to the seed of the map
            random.seed(self.__seed)

            try:
                if self.__wind_nnodes < 0:
                    __wind_nnodes = random.choice(range(4,10))
            except:
                __wind_nnodes = self.__wind_nnodes

            # Prepare wind nodes taking some random points
            __wind_nodes = []
            tweight = 0.

            # Locate nodes
            for kk in range(__wind_nnodes):

                lat = random.uniform(-89.,89.)*self.__dera
                lon = random.uniform(-180,179.)*self.__dera
                sign = random.choice([-1.,1.])
                weight = random.uniform(0.,2.)
                __wind_nodes.append([lat, lon, sign, weight])
                tweight += weight

            self.__wind_nnodes = __wind_nnodes
            self.__wind_nodes = copy.deepcopy(__wind_nodes)

        # If no random
        else:

            self.__wind_nnodes = 4

            # Prepare wind nodes taking some points
            self.__wind_nodes = [ \
                 [-67*self.__dera,-160*self.__dera, 1,.4], \
                 [-15*self.__dera, 120*self.__dera, 1,.3], \
                 [ 45*self.__dera,  20*self.__dera,-1,.2], \
                 [ 61*self.__dera, -80*self.__dera,-1,.1]]

######################################################################
######################################################################

    def get_windnodes(self):
        ''' Get values of wind nodes
        '''

        return self.__wind_nodes

######################################################################
######################################################################

    def clear_windnodes(self):
        ''' Reset wind nodesa
        '''

        self.__wind_nodes = None


######################################################################
######################################################################

    def set_maxwindspeed(self,wind):
        ''' Set value of maxwindspeed
        '''

        try:
            wind = float(wind)
            if wind <= 0.:
                msg = 'maxwindspeed must be ' + \
                      'positive.'
                self.__warning(msg)
            else:
                self.__maxwindspeed = wind
        except ValueError:
            msg = 'ValueError in maxwindspeed'
            self.__error(msg)
        except:
            msg = 'Unexpected error'
            error = sys.exc_info()[:2]
            self.__error(msg, error)

######################################################################
######################################################################

    def get_maxwindspeed(self):
        ''' Get value of maxwindspeed
        '''

        return self.__maxwindspeed

######################################################################
######################################################################

    def set_mintemperature(self,temp):
        ''' Set value of mintemperature
        '''

        try:
            temp = float(temp)
            if temp > 0.:
                msg = 'mintemperature must be ' + \
                      'negative.'
                self.__warning(msg)
            else:
                self.__mintemperature = temp
        except ValueError:
            msg = 'ValueError in mintemperature'
            self.__error(msg)
        except:
            msg = 'Unexpected error'
            error = sys.exc_info()[:2]
            self.__error(msg, error)

######################################################################
######################################################################

    def get_mintemperature(self):
        ''' Get value of mintemperature
        '''

        return self.__mintemperature

######################################################################
######################################################################

    def set_maxtemperature(self,temp):
        ''' Set value of maxtemperature
        '''

        try:
            temp = float(temp)
            if temp < 0.:
                msg = 'maxtemperature must be ' + \
                      'positive.'
                self.__warning(msg)
            else:
                self.__maxtemperature = temp
        except ValueError:
            msg = 'ValueError in maxtemperature'
            self.__error(msg)
        except:
            msg = 'Unexpected error'
            error = sys.exc_info()[:2]
            self.__error(msg, error)

######################################################################
######################################################################

    def get_maxtemperature(self):
        ''' Get value of maxtemperature
        '''

        return self.__maxtemperature

######################################################################
######################################################################

    def set_name(self,name):
        ''' Set value of name
        '''

        self.__name = str(name)

######################################################################
######################################################################

    def get_name(self):
        ''' Get value of name
        '''

        return self.__name

######################################################################
######################################################################

    def set_projection(self,proj):
        ''' Set map projection
        '''

        if proj in self.__modes.keys():
            self.__projection = proj
        else:
            msg = 'Unknown projection'
            self.__warning(msg)

######################################################################
######################################################################

    def get_projection(self):
        ''' Get value of projection
        '''

        return self.__projection

######################################################################
######################################################################

    def help_projection(self):
        ''' Show available projections
        '''

        for key in self.__modes.keys():
            msg = '{0}:{1}'.format(key,self.__modes[key])

######################################################################
######################################################################

    def get_exist(self):
        ''' Checks if there is a map
        '''

        return 'Yes' if self.__exists else 'No'

######################################################################
######################################################################
######################################################################
######################################################################

    def generate_all(self):
        ''' Generates the full weather maps
        '''

        self.generate_map(full=True)

        if not self.__exist:
            msg = 'Heights could not be generated, abort order'
            self.__error(msg)
            return

        if not self.__fullmapx and not self.__fullmapy:
            msg = 'Map must be full globe to do weather and rivers.'
            self.__error(msg)
            self.__exist = False
            return

        self.generate_wind()

        if not self.__wind_exist:
            msg = 'Wind could not be generated, abort order'
            self.__error(msg)
            self.__exist = False
            return

        self.generate_temperature()

        if not self.__temperature_exist:
            msg = 'Temperature could not be generated, abort ' + \
                  'order'
            self.__error(msg)
            self.__exist = False
            self.__wind_exist = False
            return

        self.generate_moist()

        if not self.__moist_exist:
            msg = 'Moisture could not be generated, abort ' + \
                  'order'
            self.__error(msg)
            self.__exist = False
            self.__wind_exist = False
            self.__temperature_exist = False
            return

        self.generate_biome()

        if not self.__biome_exist:
            msg = 'Biome could not be generated, abort ' + \
                  'order'
            self.__error(msg)
            self.__exist = False
            self.__wind_exist = False
            self.__temperature_exist = False
            self.__moist_exist = False
            return

        self.generate_rivers()

        if not self.__river_exist:
            msg = 'Rivers could not be generated, abort ' + \
                  'order'
            self.__error(msg)
            self.__exist = False
            self.__wind_exist = False
            self.__temperature_exist = False
            self.__moist_exist = False
            self.__biome_exist = False
            self.__river_exist = False
            self.__lake_exist = False
            return

######################################################################
######################################################################

    def generate_weather(self):
        ''' Generates the full weather maps
        '''

        if not self.__exist:
            msg = 'Must generate a map first'
            self.__error(msg)
            return

        if not self.__fullmapx and not self.__fullmapy:
            msg = 'Map must be full globe'
            self.__error(msg)
            return

        self.generate_wind()

        if not self.__wind_exist:
            msg = 'Wind could not be generated, abort weather order'
            self.__error(msg)
            return

        self.generate_temperature()

        if not self.__temperature_exist:
            msg = 'Temperature could not be generated, abort ' + \
                  'weather order'
            self.__error(msg)
            self.__wind_exist = False
            return

        self.generate_moist()

        if not self.__moist_exist:
            msg = 'Moisture could not be generated, abort ' + \
                  'weather order'
            self.__error(msg)
            self.__wind_exist = False
            self.__temperature_exist = False
            return

        self.generate_biome()

        if not self.__biome_exist:
            msg = 'Biome could not be generated, abort ' + \
                  'weather order'
            self.__error(msg)
            self.__wind_exist = False
            self.__temperature_exist = False
            self.__moist_exist = False
            return

######################################################################
######################################################################

    def __generate_map(self, Lat, Lon, seed, persistence, octaves, \
                       frequency, shift, depth, height):
        ''' Generates a map for the current parameters
        '''

        # Translate Lat and Lon
        iLat = (90.-Lat)*self.__dera
        nth = iLat.size
        iLon = (Lon+180.)*self.__dera
        nch = iLon.size

        # Translate persistence
        ipersistence = 1./persistence

        # Size of noise
        pnum = nch*nth

        # Compute noise
        simplex = simplex_class(octaves=octaves, \
                                persistence=ipersistence, \
                                scale=frequency)
        noise = []
        for lat in iLat:
            for lon in iLon:
                x = np.sin(lat)*np.cos(lon) + 1.
                y = np.sin(lat)*np.sin(lon) + 1.
                z = np.cos(lat) + 1.
                noise.append( \
                      simplex.scaled_octave_noise_35d(x,y,z, \
                                                      seed))

        # Water control
        for kk in range(pnum):
            noise[kk] += shift
            if noise[kk] < 0.:
                noise[kk] *= depth
            else:
                noise[kk] *= height

        # Return map
        return np.array(noise).reshape(nth,nch)

######################################################################
######################################################################

    def generate_map(self, full=None, silent=False, refine=False):
        ''' Generates a map for the current parameters
        '''

        # Check silent
        if not isinstance(silent, bool):
            silent = False

        # Check refine
        if not isinstance(refine, bool):
            refine = False

        # Check if there is a previous map
        if refine and not self.__exist:
            msg = 'cannot use refine flag if there is no ' + \
                  'previous map'
            self.__error(msg)
            return

        # Translate thran and chran
        thrange = [self.__thrange[0],self.__thrange[1]]
        chrange = [self.__chrange[0],self.__chrange[1]]
        if abs(thrange[0]+90.) < 1e-4:
            thrange[0] = -89.9999
        if abs(thrange[1]-90.) < 1e-4:
            thrange[1] =  89.9999
        for ii in range(2):
            chrange[ii] = (chrange[ii]+180.)*self.__dera
            thrange[ii] = (90.-thrange[ii])*self.__dera

        # Translate persistence
        persistence = 1./self.__persistence

        # Size of noise
        pnum = self.__nch*self.__nth

        # Longitude and latitude
        Lat = np.linspace(np.cos(thrange[0]), \
                          np.cos(thrange[1]), self.__nth, \
                          endpoint=True)

        # If full map
        if np.absolute(np.absolute(thrange[1] - thrange[0]) - \
                       np.pi) < 1e-3:
            self.__fullmapy = True
        else:
            self.__fullmapy = False
        if np.absolute(chrange[1] - chrange[0] - np.pi*2.) < 1e-3:
            step = 2.*np.pi/self.__nch
            Tch = chrange[1] - step
            Lon = np.linspace(chrange[0],Tch,self.__nch,endpoint=True)
            self.__fullmapx = True
        else:
            Lon = np.linspace(chrange[0],chrange[1],self.__nch)
            self.__fullmapx = False

        # Check if a full map was required here
        if full is not None:
            # Check input type
            if isinstance(full, bool):
                # If requiring full map
                if full:
                    if not self.__fullmapx or not self.__fullmapy:
                        msg = 'a full map is required, adjust ' + \
                              'thrange and chrange'
                        self.__error(msg)
                        return
            # Wrong input type
            else:
                msg = 'full must be bool'
                self.__error(msg)
                return

        # Massage latitude
        for ii in range(self.__nth):
            Lat[ii] = np.arccos(Lat[ii])
        if abs(Lat[0]-np.pi) < 1e-2:
            if Lat[1] < (np.pi-0.001):
                Lat[0] = np.pi - 0.001
            else:
                Lat[0] = 0.5*(Lat[0] + Lat[1])
        if abs(Lat[-1]) < 1e-2:
            if Lat[-2] > 0.001:
                Lat[-1] = 0.001
            else:
                Lat[-1] = 0.5*(Lat[-1] + Lat[-2])

        # Compute noise
        msg = 'Creating height map'
        if not silent:
            self.__print(msg)
        simplex = simplex_class(octaves=self.__octaves, \
                                persistence=persistence, \
                                scale=self.__frequency)
        noise = []
        for lat in Lat:
            for lon in Lon:
                x = np.sin(lat)*np.cos(lon) + 1.
                y = np.sin(lat)*np.sin(lon) + 1.
                z = np.cos(lat) + 1.
                noise.append( \
                      simplex.scaled_octave_noise_35d(x,y,z, \
                                                      self.__seed))

        # Water control

        msg = 'Filling the world with salty water'
        if not silent:
            self.__print(msg)

        # If refining
        if refine:

            # Adjust water
            for kk in range(pnum):
                noise[kk] += self.__shift
                if noise[kk] < 0.:
                    noise[kk] *= self.__maxdepth
                else:
                    noise[kk] *= self.__maxheight

        else:

            if self.__water >= 0.:

                # If all water, shift maximum
                if abs(self.__water-100.) < 1e-2:
                    shift = -1.

                # If no water, no shift
                elif abs(self.__water) < 1e-2:
                    shift = 0.

                # If no extremes, compute shift
                else:

                    ints = 11
                    actual = 0.
                    pack = 0.
                    minv = 0.
                    maxv = 1.
                    histo = [0]*(ints - 1)
                    kk = 0

                    while abs(self.__water-actual) > 1e-2:

                        inter = np.linspace(minv,maxv,ints)

                        for nois in noise:
                            for jj in range(1,ints):

                                if nois > inter[jj-1] and \
                                   nois < inter[jj]:
                                    histo[jj-1] += 1
                                    break

                        cumold = 0.

                        for jj in range(0,ints-1):

                            cum = cumold + histo[jj]
                            actual = pack + cum*100./pnum

                            if (actual > self.__water):
                                minv = inter[jj]
                                maxv = inter[jj+1]
                                shift = -inter[jj]
                                pack = pack + cumold*100./pnum
                                actual = pack
                                break
                            else:
                                cumold = cum

                        kk += 1

                        if kk > 100: break

            else:

                shift = -0.5

            # Measure the water level 
            actual = 0.
            for kk in range(pnum):
                noise[kk] += shift
                if noise[kk] < 0.:
                    noise[kk] *= self.__maxdepth
                    actual += 1.
                else:
                    noise[kk] *= self.__maxheight

            actual = actual*100./(pnum*1.)

            # Store shift
            self.__shift = shift

            # Update water level to the real one
            if self.__water >= 0:
                self.__water = actual

        # Transform into real latutude and longitude
        for ii in range(self.__nth):
            Lat[ii] = Lat[ii]*self.__rade*(-1.) + 90.
        for ii in range(self.__nch):
            Lon[ii] = Lon[ii]*self.__rade - 180.

        # Store
        self.__lon = Lon
        self.__lat = Lat
        self.__height = np.array(noise).reshape(self.__nth,self.__nch)
        self.__exist = True

        # Update extremes
        if not refine:
            self.__minz = np.min(self.__height)
            self.__maxz = np.max(self.__height)


######################################################################
######################################################################
######################################################################
######################################################################

    def generate_wind(self, silent=False):
        ''' Generates a wind map for the current parameters.
        '''

        # Check silent
        if not isinstance(silent, bool):
            silent = False

        if not self.__exist:
            msg = 'Must generate a map first'
            self.__error(msg)
            return

        if not self.__fullmapx and not self.__fullmapy:
            msg = 'Map must be full globe'
            self.__error(msg)
            return

        try:

            msg = 'Blowing some wind'
            if not silent:
                self.__print(msg)

            # Allocate vector map
            __wind = np.zeros((self.__nth,self.__nch,3))

            # If not inputting nodes
            if self.__wind_nodes is None or \
               self.__wind_nnodes < 0:

                self.generate_windnodes()

            __wind_nnodes = self.__wind_nnodes
            __wind_nodes = copy.deepcopy(self.__wind_nodes)

            # For each latitude
            for ii,lat in zip(range(self.__nth),self.__lat):

                # Convert to radian
                la = lat*self.__dera

                # Cosine and sine
                ct = np.cos(la + np.pi*.5)
                st = np.sin(la + np.pi*.5)

                # For each longitude
                for jj,lon in zip(range(self.__nch),self.__lon):

                    # Convert to radian
                    lo = lon*self.__dera

                    # Get cyllindrical coordenates of point
                    x = lo
                    y = np.tan(lat*self.__dera)

                    # Initialize
                    wind0 = 0.
                    wind1 = 0.

                    # For each node
                    for kk in range(__wind_nnodes):

                        # Latitude and longitude
                        lav = __wind_nodes[kk][0]
                        lov = __wind_nodes[kk][1]

                        # Cosine and sin
                        ctv = np.cos(lav + np.pi*.5)
                        stv = np.sin(lav + np.pi*.5)

                        # Distance
                        dist = np.arccos(st*stv*np.cos(lo - lov) + \
                               ct*ctv)

                        # Cyllindrical coordinates
                        xv = lov
                        yv = np.tan(lav)

                        # Diference vector
                        dx = xv - x
                        dy = yv - y
                        if np.absolute(dx) > np.pi:
                            dx *= -1

                        # Sign
                        dx *= __wind_nodes[kk][2]
                        dy *= __wind_nodes[kk][2]

                        # Wind direction and normalization
                        xf = -dy
                        yf =  dx
                        rf = np.sqrt(xf*xf + yf*yf)

                        # Scale
                        ww = __wind_nodes[kk][3]/ \
                             (1. + dist)

                        if rf > 1e-7:
                            wind0 += yf*ww/rf
                            wind1 += xf*ww/rf

                    # For each node
                    for kk in range(__wind_nnodes):

                        xf = wind1
                        yf = wind0
                        rf = np.sqrt(xf*xf + yf*yf)

                        if rf > 1e-7:
                            __wind[ii,jj,0] += yf/rf
                            __wind[ii,jj,1] += xf/rf

                        # Module
                        __wind[ii,jj,2] = rf

            self.__wind = __wind*self.__maxwindspeed/ \
                          np.amax(__wind[:,:,2])
            self.__wind_exist = True

            # Update extremes
            self.__minwind = np.min(self.__wind[:,:,2])
            self.__maxwind = np.max(self.__wind[:,:,2])

        except:

            self.__wind_exist = False
            msg = 'Could not generate wind'
            error = sys.exc_info()[:2]
            self.__error(msg, error)

######################################################################
######################################################################
######################################################################
######################################################################

    def __redo_wind(self):
        ''' Re-generates a wind map for the current parameters.
        '''

        if not self.__exist:
            msg = 'Must generate a map first'
            self.__error(msg)
            return False

        try:

            # If not inputting nodes
            if self.__wind_nodes is None or \
               self.__wind_nnodes < 0:

                return False

            __wind_nnodes = self.__wind_nnodes
            __wind_nodes = copy.deepcopy(self.__wind_nodes)

            # Allocate vector map
            __wind = np.zeros((self.__nth,self.__nch,3))

            # For each latitude
            for ii,lat in zip(range(self.__nth),self.__lat):

                # Convert to radian
                la = lat*self.__dera

                # Cosine and sine
                ct = np.cos(la + np.pi*.5)
                st = np.sin(la + np.pi*.5)

                # For each longitude
                for jj,lon in zip(range(self.__nch),self.__lon):

                    # Convert to radian
                    lo = lon*self.__dera

                    # Get cyllindrical coordenates of point
                    x = lo
                    y = np.tan(lat*self.__dera)

                    # Initialize
                    wind0 = 0.
                    wind1 = 0.

                    # For each node
                    for kk in range(__wind_nnodes):

                        # Latitude and longitude
                        lav = __wind_nodes[kk][0]
                        lov = __wind_nodes[kk][1]

                        # Cosine and sin
                        ctv = np.cos(lav + np.pi*.5)
                        stv = np.sin(lav + np.pi*.5)

                        # Distance
                        dist = np.arccos(st*stv*np.cos(lo - lov) + \
                               ct*ctv)

                        # Cyllindrical coordinates
                        xv = lov
                        yv = np.tan(lav)

                        # Diference vector
                        dx = xv - x
                        dy = yv - y
                        if np.absolute(dx) > np.pi:
                            dx *= -1

                        # Sign
                        dx *= __wind_nodes[kk][2]
                        dy *= __wind_nodes[kk][2]

                        # Wind direction and normalization
                        xf = -dy
                        yf =  dx
                        rf = np.sqrt(xf*xf + yf*yf)

                        # Scale
                        ww = __wind_nodes[kk][3]/ \
                             (1. + dist)

                        if rf > 1e-7:
                            wind0 += yf*ww/rf
                            wind1 += xf*ww/rf

                    # For each node
                    for kk in range(__wind_nnodes):

                        xf = wind1
                        yf = wind0
                        rf = np.sqrt(xf*xf + yf*yf)

                        if rf > 1e-7:
                            __wind[ii,jj,0] += yf/rf
                            __wind[ii,jj,1] += xf/rf

                        # Module
                        __wind[ii,jj,2] = rf

            self.__wind = __wind*self.__maxwindspeed/ \
                          np.amax(__wind[:,:,2])
            return True

        except:

            return False

######################################################################
######################################################################
######################################################################
######################################################################

    def generate_temperature(self, silent=False):
        ''' Generates a temperature map for the current parameters
        '''

        # Check silent
        if not isinstance(silent, bool):
            silent = False

        if not self.__exist:
            msg = 'Must generate a map first'
            self.__error(msg)
            return

        if not self.__fullmapx and not self.__fullmapy:
            msg = 'Map must be full globe'
            self.__error(msg)
            return

        if not self.__wind_exist:
            msg = 'Must generate winds first'
            self.__error(msg)
            return

        try:

            msg = 'Heating up'
            if not silent:
                self.__print(msg)

            # Allocate vector map
            __air_temperature = np.zeros((self.__nth,self.__nch))
            __temperature = np.zeros((self.__nth,self.__nch))
            absorb = np.zeros((self.__nth,self.__nch))

            # For each latitude
            for ii,lat in zip(range(self.__nth),self.__lat):

                # Convert to radian
                la = lat*self.__dera

                # For each longitude
                for jj,lon in zip(range(self.__nch),self.__lon):

                    # Asign initial temperature
                    __air_temperature[ii,jj] = 2.*np.cos(la)/ \
                                        (1. + np.absolute(np.sin(la)))
                    __temperature[ii,jj] = np.cos(la)/ \
                                        (1. + np.absolute(np.sin(la)))

                    absorb[ii,jj] = .01*np.amax([.1, \
                                       np.amin([self.__maxwindspeed/ \
                                               self.__wind[ii,jj,2], \
                                               2.])])

                    # Land heats more
                    if self.__height[ii,jj] > 0:

                        # Except if is very height
                        if self.__height[ii,jj] > 5:

                           absorb[ii,jj] *= 0.1

                        # Intermediate
                        if self.__height[ii,jj] >= 3.5 and \
                           self.__height[ii,jj] <= 5.:

                           absorb[ii,jj] *= 0.5

                       #else:
                        elif self.__height[ii,jj] <= 2.5:

                            absorb[ii,jj] *= 2.

            # Move the air and absorb heat
            kk = 0
            while np.amax(__air_temperature) > 0.:

                kk += 1

                # Absorb

                # For each latitude
                for ii,lat in zip(range(self.__nth),self.__lat):

                    # Convert to radian
                    la = lat*self.__dera

                    # For each longitude
                    for jj,lon in zip(range(self.__nch),self.__lon):

                        # Convert to radian
                        lo = lon*self.__dera

                        # Absorbed
                        labsorb = np.amin([absorb[ii,jj], \
                                   __air_temperature[ii,jj]])

                        # Process
                        __temperature[ii,jj] += labsorb
                        __air_temperature[ii,jj] -= labsorb

                        __temperature[ii,jj] += .0001*np.cos(la)/ \
                                        (1. + np.absolute(np.sin(la)))


                # Up to 50 iterations
                if kk > 50:
                    break

                # Move heat

                __air_temperature0 = copy.deepcopy(__air_temperature)
                __air_temperature = np.zeros(__air_temperature0.shape)

                # For each latitude
                for ii,lat in zip(range(self.__nth),self.__lat):

                    # Convert to radian
                    la = lat*self.__dera

                    # For each longitude
                    for jj,lon in zip(range(self.__nch),self.__lon):

                        # Convert to radian
                        lo = lon*self.__dera

                        x = lo
                        y = np.tan(lat*self.__dera)

                        # Direction of wind
                        if self.__maxwindspeed > 0.:
                            vy = self.__wind[ii,jj,0]/ \
                                 self.__maxwindspeed
                            vx = self.__wind[ii,jj,1]/ \
                                 self.__maxwindspeed
                        else:
                            vx = 0.
                            vy = 0.
                        sx, sy = np.sign([vx,vy])
                        mx, my = np.absolute([vx,vy])

                        # No transfer
                        if mx < 1e-5 and my < 1e-5:

                            continue

                        # Horizontal transfer
                        elif my < 1e-5:

                            jja = int(jj + sx)

                            if jja > self.__nch - 1:
                                jja = 0


                            __air_temperature[ii,jja] += \
                                             __air_temperature0[ii,jj]


                        # Vertical transfer
                        elif mx < 1e-5:

                            iia = int(ii + sy)

                            if iia > self.__nth - 1:
                                iia = 0

                            __air_temperature[iia,jj] += \
                                             __air_temperature0[ii,jj]


                        # Border latitudes
                        elif (ii == 0 and sy < 0) or \
                             (ii == self.__nth-1 and sy > 0):

                            iia = ii
                            dy = 2.*(1000. - np.tan(self.__lat[ii]* \
                                                    self.__dera))

                            # Even number of longitudes
                            if self.__nch % 2 == 0:

                                jjw0 = self.__nch/2 + jj
                                jja = int(jjw0 + sx)
                                if jjw0 > self.__nch-1:
                                    jjw0 -= self.__nch
                                jjw1 = jjw0
                                if jja > self.__nch-1:
                                    jja -= self.__nch

                                x = self.__lon[jjw0]*self.__dera

                                if (jjw0 == 0 and sx < 0) or \
                                   (jjw0 == self.__nch-1 and sx > 0):

                                    dx = 2.*np.pi - \
                                        self.__lon[-1]*self.__dera + \
                                         self.__lon[0]*self.__dera

                                else:

                                    dx = \
                                       self.__lon[jja]*self.__dera - x
                                    dxb = 0.

                            # Odd number of longitudes
                            else:

                                jjw0 = self.__nch/2 + jj
                                jjw1 = jjw0 + 1
                                if jjw0 > self.__nch-1:
                                    jjw0 -= self.__nch
                                if jjw1 > self.__nch-1:
                                    jjw1 -= self.__nch
                                if sx > 0:
                                    jja = jjw1
                                else:
                                    jja = jjw0

                                if (sx > 0 and jjw1 == 0) or \
                                   (sx < 0 and jjw0 == self.__nch-1):

                                    dx = np.pi - .5* \
                                                 self.__lon[jjw0]* \
                                                 self.__dera

                                else:

                                    dx = .5*(self.__lon[jjw1] - \
                                             self.__lon[jjw0])* \
                                         self.__dera
                                    dxb = dx

                            dyv = dx*vy/vx

                            # Cut X
                            if np.absolute(dyv) > np.absolute(dy):

                                dxv = dy*vx/vy

                                dr1 = np.absolute((dxb + dxv)/ \
                                                  (dxb + dx))
                                dr0 = np.absolute((dx - dxv)/dx)

                                __air_temperature[iia,jjw0] += \
                                         dr0*__air_temperature0[ii,jj]
                                __air_temperature[iia,jjw1] += \
                                         dr1*__air_temperature0[ii,jj]

                            # Cut Y
                            else:

                                dr1 = np.absolute(dyv/dy)
                                dr0 = np.absolute((dy - dyv)/dy)

                                __air_temperature[ii,jja] += \
                                         dr0*__air_temperature0[ii,jj]
                                __air_temperature[iia,jja] += \
                                         dr1*__air_temperature0[ii,jj]

                        # General transfer
                        else:

                            jja = int(jj + sx)
                            iia = int(ii + sy)
                            if jja > self.__nch - 1:
                                jja = 0
                            if iia > self.__nth - 1:
                                iia = 0

                            if (jj == 0 and sx < 0) or \
                               (jj == self.__nch-1 and sx > 0):

                                dx = 2.*np.pi - \
                                     self.__lon[-1]*self.__dera + \
                                     self.__lon[0]*self.__dera

                                dx *= sx

                            else:

                                dx = self.__lon[jja]*self.__dera - x

                            dy = \
                               np.tan(self.__lat[iia]*self.__dera) - y

                            dyv = dx*vy/vx

                            # Cut X
                            if np.absolute(dyv) > np.absolute(dy):

                                dxv = dy*vx/vy

                                dr1 = np.absolute(dxv/dx)
                                dr0 = np.absolute((dx - dxv)/dx)

                                __air_temperature[iia,jj] += \
                                         dr0*__air_temperature0[ii,jj]
                                __air_temperature[iia,jja] += \
                                         dr1*__air_temperature0[ii,jj]

                            # Cut Y
                            else:

                                dr1 = np.absolute(dyv/dy)
                                dr0 = np.absolute((dy - dyv)/dy)

                                __air_temperature[ii,jja] += \
                                         dr0*__air_temperature0[ii,jj]
                                __air_temperature[iia,jja] += \
                                         dr1*__air_temperature0[ii,jj]


            # For each latitude
            for ii,lat in zip(range(self.__nth),self.__lat):

                # Convert to radian
                la = lat*self.__dera

                # For each longitude
                for jj,lon in zip(range(self.__nch),self.__lon):

                    # Except if is very height
                    if self.__height[ii,jj] > 5:

                        __temperature[ii,jj] -= 2.


                    if self.__height[ii,jj] >= 3.5 and \
                       self.__height[ii,jj] <= 5.:

                        __temperature[ii,jj] -= 1.

                    else:

                        __temperature[ii,jj] -= 1.*np.sin(la)* \
                                                   np.sin(la)

               #__temperature[ii,:] -= 1.*np.sin(la)*np.sin(la)

            self.__temperature = __temperature

            pointsp = np.where(self.__temperature >= 0.)
            if np.amax(self.__temperature) > 0.:
                self.__temperature[pointsp] *= \
                                              self.__maxtemperature/ \
                                np.amax(self.__temperature[pointsp])

            pointsm = np.where(self.__temperature < 0.)
            if np.amin(self.__temperature) < 0.:
                self.__temperature[pointsm] *= \
                                             -self.__mintemperature/ \
                     np.amax(np.absolute(self.__temperature[pointsm]))

            self.__temperature_exist = True

            # Update extremes
            self.__mintemperature = np.min(self.__temperature)
            self.__maxtemperature = np.max(self.__temperature)
 
        except:
            self.__temperature_exist = False
            msg = 'Could not generate temperature'
            error = sys.exc_info()[:2]
            self.__error(msg, error)

######################################################################
######################################################################
######################################################################
######################################################################

    def generate_moist(self, silent=False):
        ''' Generates a moisture map for the current parameters
            Experimental.
        '''

        # Check silent
        if not isinstance(silent, bool):
            silent = False

        if not self.__exist:
            msg = 'Must generate a map first'
            self.__error(msg)
            return

        if not self.__fullmapx and not self.__fullmapy:
            msg = 'Map must be full globe'
            self.__error(msg)
            return

        if not self.__wind_exist:
            msg = 'Must generate winds first'
            self.__error(msg)
            return

        if not self.__temperature_exist:
            msg = 'Must generate temperature first'
            self.__error(msg)
            return

        try:

            msg = 'Watering the land a bit'
            if not silent:
                self.__print(msg)

            # Resolution factor
            rfac = 1.

            # Allocate vector map
            __air_moist = np.zeros((self.__nth,self.__nch))
            __moist = np.zeros((self.__nth,self.__nch))
            absorb = np.zeros((self.__nth,self.__nch))

            # For each latitude
            for ii,lat in zip(range(self.__nth),self.__lat):

                # Convert to radian
                la = lat*self.__dera

                # For each longitude
                for jj,lon in zip(range(self.__nch),self.__lon):

                    # Asign initial moist

                    height = self.__height[ii,jj]
                    temp = self.__temperature[ii,jj]

                    if height > 0:
                        land = True
                    else:
                        land = False

                    if land:

                        if height < 1.5:
                           #absorb[ii,jj] = .005
                            absorb[ii,jj] = .01
                        elif height < 2.:
                            absorb[ii,jj] = .025
                        else:
                            absorb[ii,jj] = .075

                    else:

                        __moist[ii,jj] = 1.

                        if temp < -.5:
                            __air_moist[ii,jj] = .05
                        elif temp < 0.:
                            __air_moist[ii,jj] = .1
                        elif temp < 5.:
                            __air_moist[ii,jj] = .3
                        elif temp < 10.:
                            __air_moist[ii,jj] = .5
                        elif temp < 15.:
                            __air_moist[ii,jj] = .8
                        elif temp < 20.:
                            __air_moist[ii,jj] = .9
                        elif temp < 30.:
                            __air_moist[ii,jj] = 1.
                        else:
                            __air_moist[ii,jj] = 1.2

            # Apply resolution factor to absorption
            absorb *= rfac

            kk = 0
            while np.amax(__air_moist) > 0.:

                kk += 1

                # Absorb

                # For each latitude
                for ii,lat in zip(range(self.__nth),self.__lat):

                    # Convert to radian
                    la = lat*self.__dera

                    # For each longitude
                    for jj,lon in zip(range(self.__nch),self.__lon):

                        # Convert to radian
                        lo = lon*self.__dera

                        height = self.__height[ii,jj]

                        if height > 0:
                            land = True
                        else:
                            land = False

                        if land:

                            # Absorbed
                            labsorb = np.amin([absorb[ii,jj], \
                                       __air_moist[ii,jj], \
                                  np.amax([0.,1. - __moist[ii,jj]])])

                            # Process
                            __moist[ii,jj] += labsorb
                            __air_moist[ii,jj] -= labsorb

                if kk > 50:
                    break

                # Move moist

                __air_moist0 = copy.deepcopy(__air_moist)
                __air_moist = np.zeros(__air_moist0.shape)


                # For each latitude
                for ii,lat in zip(range(self.__nth),self.__lat):

                    # Convert to radian
                    la = lat*self.__dera

                    # For each longitude
                    for jj,lon in zip(range(self.__nch),self.__lon):

                        # Convert to radian
                        lo = lon*self.__dera

                        x = lo
                        y = np.tan(lat*self.__dera)

                        # Direction of wind
                        if self.__maxwindspeed > 0.:
                            vy = self.__wind[ii,jj,0]/ \
                                 self.__maxwindspeed
                            vx = self.__wind[ii,jj,1]/ \
                                 self.__maxwindspeed
                        else:
                            vx = 0.
                            vy = 0.
                        sx, sy = np.sign([vx,vy])
                        mx, my = np.absolute([vx,vy])

                        # No transfer
                        if mx < 1e-5 and my < 1e-5:

                            continue

                        # Horizontal transfer
                        elif my < 1e-5:

                            jja = int(jj + sx)

                            if jja > self.__nch - 1:
                                jja = 0
                                
                            __air_moist[ii,jja] += __air_moist0[ii,jj]

                        # Vertical transfer
                        elif mx < 1e-5:

                            iia = int(ii + sy)

                            if iia > self.__nth - 1:
                                iia = 0
                                
                            __air_moist[iia,jj] += __air_moist0[ii,jj]

                        # Border latitudes
                        elif (ii == 0 and sy < 0) or \
                             (ii == self.__nth-1 and sy > 0):

                            iia = ii
                            dy = 2.*(1000. - np.tan(self.__lat[ii]* \
                                                    self.__dera))

                            # Even number of longitudes
                            if self.__nch % 2 == 0:

                                jjw0 = self.__nch/2 + jj
                                jja = int(jjw0 + sx)
                                if jjw0 > self.__nch-1:
                                    jjw0 -= self.__nch
                                jjw1 = jjw0
                                if jja > self.__nch-1:
                                    jja -= self.__nch

                                x = self.__lon[jjw0]*self.__dera

                                if (jjw0 == 0 and sx < 0) or \
                                   (jjw0 == self.__nch-1 and sx > 0):

                                    dx = 2.*np.pi - \
                                         self.__lon[-1]* \
                                         self.__dera + \
                                         self.__lon[0]*self.__dera

                                else:

                                    dx = \
                                       self.__lon[jja]*self.__dera - x
                                    dxb = 0.

                            # Odd number of longitudes
                            else:

                                jjw0 = self.__nch/2 + jj
                                jjw1 = jjw0 + 1
                                if jjw0 > self.__nch-1:
                                    jjw0 -= self.__nch
                                if jjw1 > self.__nch-1:
                                    jjw1 -= self.__nch
                                if sx > 0:
                                    jja = jjw1
                                else:
                                    jja = jjw0

                                if (sx > 0 and jjw1 == 0) or \
                                   (sx < 0 and jjw0 == self.__nch-1):

                                    dx = np.pi - .5* \
                                                 self.__lon[jjw0]* \
                                                 self.__dera

                                else:

                                    dx = .5*(self.__lon[jjw1] - \
                                             self.__lon[jjw0])* \
                                         self.__dera
                                    dxb = dx

                            dyv = dx*vy/vx

                            # Cut X
                            if np.absolute(dyv) > np.absolute(dy):

                                dxv = dy*vx/vy

                                dr1 = np.absolute((dxb + dxv)/ \
                                                  (dxb + dx))
                                dr0 = np.absolute((dx - dxv)/dx)

                                __air_moist[iia,jjw0] += dr0* \
                                                   __air_moist0[ii,jj]
                                __air_moist[iia,jjw1] += dr1* \
                                                   __air_moist0[ii,jj]

                            # Cut Y
                            else:

                                dr1 = np.absolute(dyv/dy)
                                dr0 = np.absolute((dy - dyv)/dy)

                                __air_moist[ii,jja] += dr0* \
                                                   __air_moist0[ii,jj]
                                __air_moist[iia,jja] += dr1* \
                                                   __air_moist0[ii,jj]

                        # General transfer
                        else:

                            jja = int(jj + sx)
                            iia = int(ii + sy)
                            if jja > self.__nch - 1:
                                jja = 0
                            if iia > self.__nth - 1:
                                iia = 0

                            if (jj == 0 and sx < 0) or \
                               (jj == self.__nch-1 and sx > 0):

                                dx = 2.*np.pi - \
                                     self.__lon[-1]*self.__dera + \
                                     self.__lon[0]*self.__dera

                                dx *= sx

                            else:

                                dx = self.__lon[jja]*self.__dera - x

                            dy = \
                               np.tan(self.__lat[iia]*self.__dera) - y

                            dyv = dx*vy/vx

                            # Cut X
                            if np.absolute(dyv) > np.absolute(dy):

                                dxv = dy*vx/vy

                                dr1 = np.absolute(dxv/dx)
                                dr0 = np.absolute((dx - dxv)/dx)

                                __air_moist[iia,jj] += dr0* \
                                                   __air_moist0[ii,jj]
                                __air_moist[iia,jja] += dr1* \
                                                   __air_moist0[ii,jj]

                            # Cut Y
                            else:

                                dr1 = np.absolute(dyv/dy)
                                dr0 = np.absolute((dy - dyv)/dy)

                                __air_moist[ii,jja] += dr0* \
                                                   __air_moist0[ii,jj]
                                __air_moist[iia,jja] += dr1* \
                                                   __air_moist0[ii,jj]

            self.__moist = __moist*1e2

            self.__moist_exist = True

            # Update extremes
            self.__minmoist = np.min(self.__moist)
            self.__maxmoist = np.max(self.__moist)
 
        except:
            self.__moist_exist = False
            msg = 'Could not generate moist'
            error = sys.exc_info()[:2]
            self.__error(msg, error)

######################################################################
######################################################################
######################################################################
######################################################################

    def generate_biome(self, silent=False):
        ''' Generates a biome map for the current parameters
            Experimental.
        '''

        # Check silent
        if not isinstance(silent, bool):
            silent = False

        if not self.__exist:
            msg = 'Must generate a map first'
            self.__error(msg)
            return

        if not self.__wind_exist:
            msg = 'Must generate winds first'
            self.__error(msg)
            return

        if not self.__temperature_exist:
            msg = 'Must generate temperature first'
            self.__error(msg)
            return

        if not self.__moist_exist:
            msg = 'Must generate moisture first'
            self.__error(msg)
            return

        try:

            msg = 'Clasifying land'
            if not silent:
                self.__print(msg)

            # Allocate vector map
            __biome = np.zeros((self.__nth,self.__nch),dtype=np.int16)

            # For each latitude
            for ii in range(self.__nth):

                # For each longitude
                for jj in range(self.__nch):

                    # Asign initial biome

                    height = self.__height[ii,jj]
                    temp = self.__temperature[ii,jj]
                    moist = self.__moist[ii,jj]

                    if height < 0:

                        __biome[ii,jj] = -3

                    elif height > 3:

                        if temp <= 0.:
                            __biome[ii,jj] = -2
                        else:
                            __biome[ii,jj] = -1

                    else:

                        if temp < 0:
                            __biome[ii,jj] = 0

                        elif temp > 20:

                            if moist > 65:
                                __biome[ii,jj] = 3

                            elif moist > (1.8*temp-24.):
                                __biome[ii,jj] = 5

                            else:
                                __biome[ii,jj] = 8

                        else:

                            if temp < 5:

                                if moist < (0.35*temp+5):
                                    __biome[ii,jj] = 7

                                elif moist > (4.*temp+5):
                                    __biome[ii,jj] = 1

                                else:
                                    __biome[ii,jj] = 6

                            else:

                                if moist > 50:
                                    __biome[ii,jj] = 2

                                elif moist > 25:
                                    __biome[ii,jj] = 4

                                elif moist > (0.35*temp+5):
                                    __biome[ii,jj] = 6

                                else:
                                    __biome[ii,jj] = 7

            self.__biome = __biome
            self.__biome_exist = True
 
        except:
            self.__biome_exist = False
            msg = 'Could not generate biomes'
            error = sys.exc_info()[:2]
            self.__error(msg, error)


######################################################################
######################################################################
######################################################################
######################################################################

    def generate_rivers(self, detailed=None, silent=False):
        ''' Generates rivers/lakes
            Extremely experimental.
        '''

        # Check silent
        if not isinstance(silent, bool):
            silent = False

        if not self.__exist:
            msg = 'Must generate a map first'
            self.__error(msg)
            return

        if not self.__fullmapx and not self.__fullmapy:
            msg = 'Map must be full globe'
            self.__error(msg)
            return

        if not self.__wind_exist:
            msg = 'Must generate winds first'
            self.__error(msg)
            return

        if not self.__temperature_exist:
            msg = 'Must generate temperature first'
            self.__error(msg)
            return

        if not self.__moist_exist:
            msg = 'Must generate moisture first'
            self.__error(msg)
            return

        if not self.__biome_exist:
            msg = 'Must generate biome first'
            self.__error(msg)
            return

        if not self.__interp:
            msg = 'An interpolator is required'
            self.__error(msg)
            return

        if not self.__random:
            msg = 'The random generator is required'
            self.__error(msg)
            return

        if detailed is None:
            idetailed = False
        else:
            if isinstance(detailed, bool):
                idetailed = detailed
            else:
                msg = 'detailed must be bool, taking default'
                self.__warning(msg)
                idetailed = False

        try:

            msg = 'Making that water flow'
            if not silent:
                self.__print(msg)

            # Initialize river list
            __river = []
            __lake = []
            __check = [[],[]]

            # Initialize seed to the seed of the map
            random.seed(self.__seed)
            np.random.seed(self.__seed)

            # Number of sources to try
            msource = 500
            Msource = 5000
            if self.__random:
                Nsource = random.choice(range(msource,Msource+1))
            else:
                Nsource = int(round(.5*(msource + Msource)))

            # Pointer to height
            Npoints = self.__nth*self.__nch

            # Biomes that cannot be source
            BadSource = [-3,0,1,8]

            # Create an interpolation function
            self.__fint = \
                     interpolate.RegularGridInterpolator( \
                              (self.__lat, self.__lon), self.__height)

            # For each source, get a random point in the map
            for isource in range(Nsource):

                # Choose random point
                ipoint = random.choice(range(Npoints))

                # Decouple indexes
                ilat = ipoint/self.__nch
                ilon = ipoint % self.__nch

                # Get info from this point
                lh = self.__height[ilat,ilon]
                lb = self.__biome[ilat,ilon]

                # If ocean, skip
                if lh <= 0.:
                    continue

                # If desert, skip
                if lb in BadSource:
                    continue

                # Create a river and its lakes
                river, lake = self.__create_river(ilat, ilon, \
                                                  __river, idetailed)

                if len(river) > 0:
                    __river.append(river)
                if len(lake) > 0:
                    lakex = lake[0]
                    lakey = lake[1]
                    for lakx,laky in zip(lake[0],lake[1]):
                        __lake.append([lakx,laky])

            # Save rivers
            self.__river = __river
            if len(self.__river) > 0:
                self.__river_exist = True

            # Correct and save lakes, if any
            self.__lake = __lake
            if len(self.__lake) > 0:
                self.__lake_noise()
                self.__lake_exist = True

        except:
            self.__river_exist = False
            self.__lake_exist = False
            msg =  '##Error## Could not generate rivers'
            error = sys.exc_info()[:2]
            self.__error(msg, error)

######################################################################
######################################################################
######################################################################

    def __create_river(self, ilat, ilon, gcheck, detailed):
        ''' Creates a river stream
        '''

        # Distance in degrees to consider a river touching an existing
        # one
        driver = .001

        # Distance in degrees between origins to compare with a
        # previous river
        doririver = 30.

        # Number of steps ahead to check after a lake
        Nahead = 10

        # Gaussian sigma to smooth the 2D zooms after interpolation
        sigma = 10.

        # Initialize river origin
        lorigin = [ilat,ilon]

        # Go up to find a higher origin
        while True:

            # Look for the upstream
            up = self.__upstream(lorigin[0], lorigin[1])

            # If going up, update origin
            if len(up) > 0:
                lorigin = up

            # If we are done, exit the loop
            else:
                break

        # Get resolution distance squared
        driver2 = driver*driver
        doririver2 = doririver*doririver

        # Initialize origin
        lat0 = self.__lat[lorigin[0]]
        lon0 = self.__lon[lorigin[1]]

        # Initialize stream
        outx = [lat0]
        outy = [lon0]
        check = [[lat0,lon0]]

        # Create array with bools for comparison
        compare = []
        if len(gcheck) > 0:
            for riv in gcheck:
                rivx = riv[0][0]
                rivy = riv[1][0]
                dx = rivx - lat0
                dy = rivy - lon0
                d2 = dx*dx + dy*dy
                if d2 <= doririver2:
                    compare.append(True)
                else:
                    compare.append(False)

        # Initialize lake
        loutx = []
        louty = []


        #
        # Create thinner grid around origin
        #


        # Size of box
        dlat = 5.
        nlat = 128
        dlon = 5.
        nlon = 128

        # Create thinner grid around origin
        lamin = lat0 - dlat
        if lamin < self.__lat[0]:
            lamin = self.__lat[0]
        lamax = lat0 + dlat
        if lamax > self.__lat[-1]:
            lamax = self.__lat[-1]
        x = np.linspace(lamin, lamax, num=nlat, endpoint=True, \
                        dtype=np.float)

        lomin = lon0 - dlon
        lomax = lon0 + dlon
        y = np.linspace(lomin, lomax, num=nlon, endpoint=True, \
                        dtype=np.float)

        # Generate submap
        if detailed:
            zz = self.__generate_map(x,y,self.__seed, \
                                     self.__persistence, \
                                     self.__octaves, \
                                     self.__frequency, \
                                     self.__shift, \
                                     self.__maxdepth, \
                                     self.__maxheight)

        if lomin < self.__lon[0]:
            for iy in range(nlon):
                if y[iy] < self.__lon[0]:
                    y[iy] += self.__lon[-1] - self.__lon[0]
                else:
                    break
        if lomax > self.__lon[-1]:
            for iy in range(nlon-1,-1,-1):
                if y[iy] > self.__lon[-1]:
                    y[iy] -= self.__lon[-1] - self.__lon[0]
                else:
                    break

        # Interpolate heights to new grid
        if not detailed:
            xx, yy = np.meshgrid(x, y, indexing='ij')
            zz = self.__fint((xx,yy))

        # Smooth the map section
        if self.__smooth:

            if np.min(zz) <= 0.:
                sea = True
                inds = np.where(zz <= 0.)
            else:
                sea = False

            zz = filters.gaussian_filter(zz, sigma, mode='reflect', \
                                         truncate=4.0)
            if sea:
                if len(inds[0]) == 1:
                    zz[inds[0],inds[1]] = -1
                else:
                    for ind in inds:
                        zz[ind[0],ind[1]] = -1.


        # Look for indexes
        ilat0 = np.argmin(np.absolute(x-lat0))
        ilon0 = np.argmin(np.absolute(y-lon0))

        # Current point
        current = [ilat0,ilon0]

        # Create a lake in origin point
        pond = self.__create_lake(ilat0, ilon0, zz)

        # If we got a lake
        if len(pond) > 0:
            lpondx = []
            lpondy = []
            for po in pond:
                lpondx.append(x[po[0]])
                lpondy.append(y[po[1]])
                check.append([x[po[0]],y[po[1]]])
            loutx.append(lpondx)
            louty.append(lpondy)

        # Now go downstream until you find a minimum or the ocean
        while True:

            # Check if new interpolation is needed
            if current[0] <= 20 or current[1] <= 20 or \
               current[0] >= nlat-21 or current[1] >= nlon-21:

                # Coordinates
                lat1 = x[current[0]]
                lon1 = y[current[1]]

                # Create thinner grid around origin
                lamin = lat1 - dlat
                if lamin < self.__lat[0]:
                    lamin = self.__lat[0]
                lamax = lat1 + dlat
                if lamax > self.__lat[-1]:
                    lamax = self.__lat[-1]
                x = np.linspace(lamin, lamax, num=nlat, \
                                endpoint=True, dtype=np.float)

                lomin = lon1 - dlon
                lomax = lon1 + dlon
                y = np.linspace(lomin, lomax, num=nlon, \
                                endpoint=True, dtype=np.float)

                # Generate submap
                if detailed:
                    zz = self.__generate_map(x,y,self.__seed, \
                                             self.__persistence, \
                                             self.__octaves, \
                                             self.__frequency, \
                                             self.__shift, \
                                             self.__maxdepth, \
                                             self.__maxheight)

                if lomin < self.__lon[0]:
                    for iy in range(nlon):
                        if y[iy] < self.__lon[0]:
                            y[iy] += self.__lon[-1] - self.__lon[0]
                        else:
                            break
                if lomax > self.__lon[-1]:
                    for iy in range(nlon-1,-1,-1):
                        if y[iy] > self.__lon[-1]:
                            y[iy] -= self.__lon[-1] - self.__lon[0]
                        else:
                            break

                # Interpolate heights to new grid
                if not detailed:
                    xx, yy = np.meshgrid(x, y, indexing='ij')
                    zz = self.__fint((xx,yy))

                # Smooth the map section
                if np.min(zz) <= 0.:
                    sea = True
                    inds = np.where(zz <= 0.)
                else:
                    sea = False

                zz = filters.gaussian_filter(zz, sigma, \
                                             mode='reflect', \
                                             truncate=4.0)
                if sea:
                    if len(inds[0]) == 1:
                        zz[inds[0],inds[1]] = -1
                    else:
                        for ind in inds:
                            zz[ind[0],ind[1]] = -1.

                # Look for indexes
                ilat1 = np.argmin(np.absolute(x-lat1))
                ilon1 = np.argmin(np.absolute(y-lon1))

                # Current point
                current = [ilat1,ilon1]

                # Check if too much latitude abort river
                if current[0] == 0 or current[0] == nlat-1:
                    return [], []

            # Look for the next down point
            down, pond, cont = self.__downstream(current[0], \
                                                 current[1], zz)

            # Initialize uniqueness
            unique = True

            # If we got a lake
            if len(pond) > 0:
                lake = True
                lpondx = []
                lpondy = []
                for po in pond:
                    lpondx.append(x[po[0]])
                    lpondy.append(y[po[1]])
                    check.append([x[po[0]],y[po[1]]])
                loutx.append(lpondx)
                louty.append(lpondy)
            # If we did not get a lake
            else:
                lake = False

            # If going down, update current
            if len(down) > 0:
                current = down[-1]
                for dw in down:
                    # For some reason, something it can be empty
                    if len(dw) < 2:
                        continue
                    # If this is repeated, signal the end
                    if [x[dw[0]],y[dw[1]]] in check:
                        unique = False
                        continue
                    # If there are previous rivers
                    if len(gcheck) > 0:
                        # Check if we are very close to a previous
                        # river
                        for griv,comp in zip(gcheck,compare):
                            if not comp:
                                continue
                            grivx = griv[0]
                            grivy = griv[1]
                            for grivxl,grivyl in zip(grivx,grivy):
                                dx = x[dw[0]] - grivxl
                                dy = y[dw[1]] - grivyl
                                d2 = dx*dx + dy*dy
                                # If closer than resolution, join
                                if d2 <= driver2:
                                    outx.append(grivxl)
                                    outy.append(grivyl)
                                    unique = False
                                    cont = False
                                    break
                    outx.append(x[dw[0]])
                    outy.append(y[dw[1]])
                    check.append([x[dw[0]],y[dw[1]]])

            # If we are done, exit the loop
            else:
                break

            # If there was a lake, check if you are falling back in
            # the pool in the Nahead next iterations
            if lake:
                good = self.__checkahead(current[0], current[1], \
                                         x, y, zz, check, Nahead)
                # If we are falling back
                if not good:

                    # Remove stream of lake
                    for ii in range(3):
                        outx.pop(-1)
                        outy.pop(-1)
                    unique = False

            # If we feel back into a river/lake
            if not unique:
                break

            # If stream finished
            if not cont:
                break

        # Format the output rivers and lakes
        out = [outx,outy]
        lout = [loutx,louty]

        return out, lout

######################################################################
######################################################################
######################################################################

    def __create_lake(self, ilat, ilon, zz):
        ''' Create a lake around coordinate
        '''

        # Initialize output
        lout = []

        # Get size limits
        limx, limy = zz.shape

        # Add list of points not valid to fill
        invalid = [[ilat,ilon]]

        # Initialize a pool (lake)
        pool = [[[ilat,ilon]]]

        # Probability to fill new ring
        pring = 1.0

        # Probability to fill  neighbour cell
        pfill = .99

        # Probability reduction factor
        pfactor = .7

        # Fill until it is full
        while True:

            # Initialize if the point is fillable
            # and the local pool
            fillable = False
            lpool = []

            # See if we are filling a new ring
            doing = np.random.choice([True,False], p=[pring,1.-pring])

            # Not doing next ring
            if not doing:
                break

            # For each point in the last ring of the
            # pool
            for pol in pool[-1]:

                # Get indexes of the point
                ilat0 = pol[0]
                ilon0 = pol[1]

                # Look around
                for ila in range(-1,2):
                    for ilo in range(-1,2):

                        # Indexes point to try
                        ilat1 = ilat0 + ila
                        ilon1 = ilon0 + ilo

                        if ilat1 <= 0 or ilon1 <= 0 or \
                           ilat1 >= limx-1 or ilon1 >= limy-1:
                            break

                        # Check if invalid point
                        if [ilat1,ilon1] in invalid:
                            continue

                        # And this point is invalid now
                        invalid.append([ilat1,ilon1])

                        # Throw to fill
                        fill = np.random.choice([True,False], \
                                                p=[pfill,1.-pfill])

                        # If can be filled, add to pool
                        if fill:
                            pfill *= pfactor
                            lpool.append([ilat1,ilon1])
                            fillable = True

            # If it was possible to fill them, add to pool
            # and check if a stream can be recovered
            if fillable:

                # Add last ring to the pool
                pool.append(lpool)

            # If not fillable and there was pool, process
            else:

                # Find boundaries
                lout = self.__get_boundary_rough(pool)

                # And finished
                break

        return lout

######################################################################
######################################################################
######################################################################

    def __downstream(self, ilat, ilon, zz):
        ''' Find next point downstream
        '''

        # Initialize output
        out = []
        lout = []

        # Save entry point
        entry = [ilat,ilon]

        # Get current height
        lh = zz[ilat,ilon]

        # Initialize variables
        minV = lh
        minI = []
        WW = []
        pool = []
        found = False

        # Look around the current point
        for ila in range(-1,2):
            for ilo in range(-1,2):

                # Get coordinates and height
                ilat1 = ilat + ila
                ilon1 = ilon + ilo
                lh1 = zz[ilat1,ilon1]

                # To introduce randomness
                if lh1 < lh:
                    WW.append(lh-lh1)
                    pool.append([ilat1,ilon1])

                # Deterministic
                if lh1 < minV:
                    minV = lh1
                    minI = [ilat1,ilon1]
                    found = True

        # If reached water, finished
        if minV <= 0.:
            out = [minI]
            bout = False

        else:

            # If we found a valid minimum, continue the stream
            if found:

                # Randomness
                total = sum(WW)
                for ii in range(len(WW)):
                    WW[ii] /= total

                # Choice the upstream point
                ind = np.random.choice(range(len(WW)), 1, p=WW)[0]

                out = [pool[ind]]
                bout = True

            # If in a local minimum, fill up to 10m
            else:

                # Add list of points not valid for a stream
                invalid = [[ilat,ilon]]

                # Initialize a pool (lake)
                pool = [[[ilat,ilon]]]

                # Height parameters
                dz = 0.
                dzT = 0.1

                # Initialize flag to know if the stream
                # can be recovered
                recovered = False

                # Fill until you find a stream or it is full
                while True:

                    # Initialize if the point is fillable
                    # and the local pool
                    fillable = False
                    lpool = []

                    # For each point in the last ring of the
                    # pool
                    for pol in pool[-1]:

                        # Get indexes of the point
                        ilat0 = pol[0]
                        ilon0 = pol[1]

                        # Look around
                        for ila in range(-1,2):
                            for ilo in range(-1,2):

                                # Indexes point to try
                                ilat1 = ilat0 + ila
                                ilon1 = ilon0 + ilo

                                # Check if invalid point
                                if [ilat1,ilon1] in invalid:
                                    continue

                                # And this point is invalid now
                                invalid.append([ilat1,ilon1])

                                # Check height difference
                                dz1 = zz[ilat1,ilon1] - lh

                                # If can be filled, add to pool
                                if dz1 < dzT:
                                    if dz1 > dz:
                                        dz = dz1
                                    lpool.append([ilat1,ilon1])
                                    fillable = True

                    # If it was possible to fill them, add to pool
                    # and check if a stream can be recovered
                    if fillable:

                        # Add last ring to the pool
                        pool.append(lpool)

                        # For each point in the last ring
                        for pol in lpool:

                            # Get indexes of the point
                            ilat0 = pol[0]
                            ilon0 = pol[1]

                            lh0 = zz[ilat0,ilon0]
                            minV = lh0

                            # Look around
                            for ila in range(-1,2):
                                for ilo in range(-1,2):

                                    # Indexes point to try
                                    ilat1 = ilat0 + ila
                                    ilon1 = ilon0 + ilo

                                    # Check if invalid point
                                    if [ilat1,ilon1] in invalid:
                                        continue

                                    # Check height difference
                                    lh1 = zz[ilat1,ilon1]

                                    # If can be continued
                                    if lh1 < minV:
                                        minV = lh1
                                        stream = [ilat1,ilon1]
                                        destiny = [ilat0,ilon0]
                                        recovered = True

                    # If not fillable and there was pool, process
                    else:

                        # Find boundary
                        lout = self.__get_boundary_rough(pool)

                        # But the stream is finished
                        bout = False
                        break

                    # If it was possible to fill them, process the
                    # pool
                    if recovered:

                        # Check the next two steps
                        ilat0 = stream[0]
                        ilon0 = stream[1]

                        # Store current heights
                        lh0 = zz[ilat0,ilon0]
                        minV = lh0

                        # Initialize validity flag 1
                        check1 = False

                        # Look around
                        for ila in range(-1,2):
                            for ilo in range(-1,2):

                                # Indexes point to try
                                ilat1 = ilat0 + ila
                                ilon1 = ilon0 + ilo

                                # Check if valid stream
                                if [ilat1,ilon1] in invalid:
                                    continue

                                # Check height difference
                                lh1 = zz[ilat1,ilon1]

                                # If can be continued
                                if lh1 < minV:
                                    minV = lh1
                                    step1 = [ilat1,ilon1]
                                    check1 = True

                        # If we can make one step
                        if check1:

                            # Check the second step
                            ilat0 = step1[0]
                            ilon0 = step1[1]

                            # Store current heights
                            lh0 = zz[ilat0,ilon0]
                            minV = lh0

                            # Initialize validity flag 2
                            check2 = False

                            # Look around
                            for ila in range(-1,2):
                                for ilo in range(-1,2):

                                    # Indexes point to try
                                    ilat1 = ilat0 + ila
                                    ilon1 = ilon0 + ilo

                                    # Check if valid stream
                                    if [ilat1,ilon1] in invalid:
                                        continue

                                    # Check height difference
                                    lh1 = zz[ilat1,ilon1]

                                    # If can be continued
                                    if lh1 < minV:
                                        minV = lh1
                                        step1 = [ilat1,ilon1]
                                        check2 = True

                        # If we can make two steps
                        if check1 and check2:

                            # Find boundary
                            lout = self.__get_boundary_rough(pool)

                            # We have a pool, join the entry with
                            # the new output

                            # Initialize auxiliar lists
                            pool1 = []

                            # Reorder the pool in a more convenient
                            # way
                            for pol in pool:
                                for p in pol:
                                    pool1.append([p[0],p[1]])

                            # Distance between stream and origin
                            dx = destiny[0] - entry[0]
                            dy = destiny[1] - entry[1]
                            d2 = dx*dx + dy*dy
                            current = copy.deepcopy(entry)

                            # Until we reach the destination in
                            # the boundary
                            while d2 > 0:

                                # Initialize dictionary of distance
                                d2_d = {}

                                # Look around
                                for ila in range(-1,2):
                                    for ilo in range(-1,2):

                                        # Indexes point to try
                                        ilat1 = current[0] + ila
                                        ilon1 = current[1] + ilo

                                        # If in pool, add to d2
                                        if [ilat1, ilon1] in \
                                           pool1 and \
                                           [ilat1, ilon1] not in \
                                           [current] and \
                                           [ilat1, ilon1] not in \
                                           out:
                                            dx = ilat1 - destiny[0]
                                            dy = ilon1 - destiny[1]
                                            d2 = dx*dx + dy*dy
                                            d2_d[d2] = [ilat1,ilon1]

                                # Look for most efficient path
                                d2_v = d2_d.keys()
                                d2 = np.min(d2_v)
                                indm = np.argmin(d2_v)
                                current = d2_d[d2_v[indm]]
                                out.append(current)

                            # And freely advance towards stream

                            # Distance between stream and origin
                            dx = stream[0] - current[0]
                            dy = stream[1] - current[1]
                            d2 = dx*dx + dy*dy

                            # Until we reach the stream point
                            while d2 > 0:

                                # Initialize dictionary of distance
                                d2_d = {}

                                # Look around
                                for ila in range(-1,2):
                                    for ilo in range(-1,2):

                                        # Indexes point to try
                                        ilat1 = current[0] + ila
                                        ilon1 = current[1] + ilo

                                        # If in pool, add to d2
                                        if [ilat1, ilon1] not in \
                                           [current] and \
                                           [ilat1, ilon1] not in \
                                           out:
                                            dx = ilat1 - stream[0]
                                            dy = ilon1 - stream[1]
                                            d2 = dx*dx + dy*dy
                                            d2_d[d2] = [ilat1,ilon1]

                                # Look for most efficient path
                                d2_v = d2_d.keys()
                                d2 = np.min(d2_v)
                                indm = np.argmin(d2_v)
                                current = d2_d[d2_v[indm]]
                                out.append(current)

                            # And continue the stream
                            bout= True
                            break

                        else:

                            recovered = False

        return out, lout, bout

######################################################################
######################################################################
######################################################################

    def __get_boundary_detail(self, area):
        ''' Get the boundary of an area
        '''

        # Initialize auxiliar lists
        listx = []
        listy = []
        pool1 = []
        poolf = []

        # Reorder the pool in a more convenient way
        for pol in area:
            for p in pol:
                listx.append(p[0])
                listy.append(p[1])
                pool1.append([p[0],p[1]])

        # Get minimum x and indexes where x has
        # such value
        minx = np.min(listx)
        indx = np.where(listx == minx)[0]

        # Look for the minimum y if more than one
        # x candidate
        if len(indx) > 1:
            llisty = []
            for ind in indx:
                llisty.append(listy[ind])
            miny = np.min(llisty)
            indy = np.where(llisty == miny)[0][0]
            indx = indx[indy]
        else:
            indx = indx[0]

        # Set current point and add to border
        poolf.append([listx[indx],listy[indx]])
        current = poolf[-1]

        # Order of priority to look for next border
        disp = [[-1,1],[0,1],[1,1],[1,0],[1,-1], \
                [0,-1],[-1,-1],[-1,0]]

        # Current displacement
        idisp = 0
        iback = 0

        counter = 0
        maxcounter = len(pool1) + 1
        stepped = True

        # Look for the boundary
        while True:

            if counter >= maxcounter:
                break

            # The displacement index at entry of loop
            if stepped:
                eidisp = idisp
                stepped = False
            
            # Get potential next point
            pnext = [current[0] + \
                     disp[idisp][0], \
                     current[1] + \
                     disp[idisp][1]]

            # If the point is in the final pool,
            # we are done
            if pnext in [poolf[0]]:
                break

            # If the point is in the pool, add to
            # the final one
            if pnext in pool1 and \
               pnext not in poolf:

                counter =0
                stepped = True
                poolf.append(pnext)
                current = pnext
                idisp -= 2
                iback = 0
                if idisp < 0:
                    idisp = len(disp) - 1

            # If not, try other displacement
            else:

                counter += 1
                idisp += 1
                if idisp >= len(disp):
                    idisp = 0
                if idisp == eidisp:
                    iback += 1
                    current = \
                          poolf[len(poolf) -1 - iback]

        # Append in the pool
        lout = poolf

        return lout

######################################################################
######################################################################
######################################################################

    def __get_boundary_rough(self, area):
        ''' Get the boundary of an area
        '''

        # Initialize auxiliar lists
        listx = []
        listy = []
        pool1 = []
        poolf = []

        # Reorder the pool in a more convenient way
        for pol in area:
            for p in pol:
                listx.append(p[0])
                listy.append(p[1])
                pool1.append([p[0],p[1]])

        #
        # Left
        #

        # Get minimum x and indexes where x has
        # such value
        minx = np.min(listx)
        indx = np.where(listx == minx)[0]

        # Look for the y
        if len(indx) > 1:
            llistx = []
            llisty = []
            for ind in indx:
                llistx.append(listx[ind])
                llisty.append(listy[ind])
            inds = np.argsort(llisty)
            for ind in inds:
                point = [llistx[ind], llisty[ind]]
                if point not in poolf:
                    poolf.append(point)
        else:
            indx = indx[0]
            if [listx[indx], listy[indx]] not in poolf:
                poolf.append([listx[indx], listy[indx]])


        #
        # Top
        #

        # Get minimum x and indexes where x has
        # such value
        maxy = np.max(listy)
        indy = np.where(listy == maxy)[0]

        # Look for the x
        if len(indy) > 1:
            llistx = []
            llisty = []
            for ind in indy:
                llistx.append(listx[ind])
                llisty.append(listy[ind])
            inds = np.argsort(llistx)
            for ind in inds:
                point = [llistx[ind], llisty[ind]]
                if point not in poolf:
                    poolf.append(point)
        else:
            indy = indy[0]
            if [listx[indy], listy[indy]] not in poolf:
                poolf.append([listx[indy], listy[indy]])

        #
        # Right
        #

        # Get maximum x and indexes where x has
        # such value
        maxx = np.max(listx)
        indx = np.where(listx == maxx)[0]

        # Look for the y
        if len(indx) > 1:
            llistx = []
            llisty = []
            for ind in indx:
                llistx.append(listx[ind])
                llisty.append(listy[ind])
            inds = np.argsort(llisty)
            inds = inds[::-1]
            for ind in inds:
                point = [llistx[ind], llisty[ind]]
                if point not in poolf:
                    poolf.append(point)
        else:
            indx = indx[0]
            if [listx[indx], listy[indx]] not in poolf:
                poolf.append([listx[indx], listy[indx]])


        #
        # Bottom
        #

        # Get minimum x and indexes where x has
        # such value
        miny = np.min(listy)
        indy = np.where(listy == miny)[0]

        # Look for the x
        if len(indy) > 1:
            llistx = []
            llisty = []
            for ind in indy:
                llistx.append(listx[ind])
                llisty.append(listy[ind])
            inds = np.argsort(llistx)
            inds = inds[::-1]
            for ind in inds:
                point = [llistx[ind], llisty[ind]]
                if point not in poolf:
                    poolf.append(point)
        else:
            indy = indy[0]
            if [listx[indy], listy[indy]] not in poolf:
                poolf.append([listx[indy], listy[indy]])

        return poolf

######################################################################
######################################################################
######################################################################

    def __checkahead(self, ilat, ilon, x, y, z, check, Nahead):
        ''' Check Nahead iterations in from to see if we are
            falling back
        '''

        # Make a copy of inputs
        iilat = ilat
        iilon = ilon
        icheck = copy.deepcopy(check)
        current = [iilat,iilon]

        # Initialize flags
        good = True
        unique = True

        # For each step ahead specified
        for ii in range(Nahead):

            # Try to go downstream
            down, pond, cont = self.__downstream(current[0], \
                                                 current[1], z)

            # If falling in another pond so close, invalid
            if len(pond) > 0:
                good = False
                break

            # If fell in water, it is good
            if not cont:
                break

            # If going down, update current
            if len(down) > 0:
                current = down[-1]
                for dw in down:
                    if [x[dw[0]],y[dw[1]]] in icheck:
                        unique = False
                        continue
                    icheck.append([x[dw[0]],y[dw[1]]])

            # If we feel back into a river/lake
            if not unique:
                good = False
                break

        return good

######################################################################
######################################################################
######################################################################

    def __upstream(self, ilat, ilon):
        ''' From a point, goes up in height, randomly with more
            probability the higher. Returns empty if no higher points
        '''

        # Initialize
        out = []
        found = False
        pool = {}
        maxH = 0.

        # Get local data
        lh = self.__height[ilat,ilon]
        lb = self.__biome[ilat,ilon]

        # If we are in a mountain, there is a 10 per cent probability
        # of starting just here
        if lb < 0:
            prob = [.1,.9]

        # If it is not a mountain, 1 per cent probability of
        # starting here
        else:
            prob = [.01,.99]

        # Throw the die
        exit = np.random.choice([True,False], 1, p=prob)

        # If we got true, this is the origin
        if exit:
            return out

        # For every point around
        for ila in range(-1,2):
            for ilo in range(-1,2):

                # Coordinates
                ilat1 = ilat + ila
                ilon1 = ilon + ilo

                # Control real index
                while ilat1 < 0:
                    ilat1 += self.__nth
                while ilon1 < 0:
                    ilon1 += self.__nch
                while ilat1 >= self.__nth:
                    ilat1 -= self.__nth
                while ilon1 >= self.__nch:
                    ilon1 -= self.__nch

                # Get height
                lh1 = self.__height[ilat1,ilon1]

                # If point higher, add to pool of candidates
                if lh1 > lh:

                    # Add to pool
                    pool[lh1] = [ilat1,ilon1]
                    found = True

                    # Update maximum height
                    if lh1 > maxH:
                        maxH = lh1

        # If we found something higher
        if found:

            # Get normalized height vector
            WW = []
            for key in pool.keys():
                WW.append(key)
            total = sum(WW)
            for ii in range(len(WW)):
                WW[ii] /= total

            # Choice the upstream point
            ind = np.random.choice(range(len(WW)), 1, p=WW)[0]

            # Define output
            out = pool[pool.keys()[ind]]

        return out

######################################################################
######################################################################
######################################################################

    def __lake_noise(self):
        ''' Adds a bit of noise to the lakes to avoid ugly squares
        '''

        # Maximum factor of minimum distance allowed to move and
        # minimum resolution of boundary
        Mdrf = .5
        Nr = 32

        # For each lake
        for ilake in range(len(self.__lake)):

            # Copy the lake
            clake = copy.deepcopy(self.__lake[ilake])

            # Split in axes
            clakex = clake[0]
            clakey = clake[1]

            # If the lake is split, it can be problematic
            if np.min(clakey) < -170. and \
               np.max(clakey) > 170.:
                reverse = True
                for iy in range(len(clakey)):
                    if clakey[iy] < 0.:
                        clakey[iy] += 359.
            else:
                reverse = False

            # Compute a center
            x0 = .5*(np.min(clakex) + np.max(clakex))
            y0 = .5*(np.min(clakey) + np.max(clakey))

            # Compute r and theta
            r = []
            theta = []
            for fx,fy in zip(clakex,clakey):
                x = fx - x0
                y = fy - y0
                lr = np.sqrt(x*x + y*y)
                lt = np.arctan2(y,x)
                r.append(lr)
                theta.append(lt)

            # Get mean radius
            Mdr = np.mean(r)

            # If too low resolution, interpolate
            if len(r) < Nr:
                 xp = copy.deepcopy(theta)
                 theta = np.linspace(-np.pi, np.pi, num=Nr, \
                                     endpoint=False, dtype=np.float)
                 r = np.interp(theta,xp,r,period=2.*np.pi)

            # Perturbate the radius
            for ir in range(len(r)):
                ldr = Mdr*random.uniform(-Mdrf,Mdrf)
                r[ir] += ldr

            # Smooth the lake border
            dr = len(r)*.03
            r = filters.gaussian_filter1d(r, dr, mode='wrap')

            # Reconstruct axes
            clakex = []
            clakey = []
            for ii in range(len(r)):
                clakex += [x0 + r[ii]*np.cos(theta[ii])]
                clakey += [y0 + r[ii]*np.sin(theta[ii])]

            # Check latitudes do not overflow
            for ix in range(len(clakex)):
                if np.absolute(clakex[ix]) >= 89.:
                    if clakex[ix] > 0.:
                        clakex[ix] = 89.
                    else:
                        clakex[ix] = -89.

            # If the lake was split, restore the longitude axis
            if reverse:
                for iy in range(len(clakey)):
                    if clakey[iy] > 179.:
                        clakey[iy] -= 359.

            # Reconstruct lake
            self.__lake[ilake] = [clakex, clakey]
        

######################################################################
######################################################################
######################################################################
######################################################################

    def save_vtk(self, name=None, R=None, only_h=None):
        ''' Stores a generated map in a vtk file
        '''

        if not self.__exist:
            msg = 'Must generate a map first'
            self.__error(msg)
            return

        if name is None:
            name = self.__name + '.vtk'

        if R is None:
            R = 6000.0

        if only_h is None:
            only_h = False

        # Wrap the map on longitudes
        if self.__fullmapx:
            pheight, elon = addcyclic(self.__height, self.__lon)
            elon[-1] = 179.99
            nch = len(elon)
        else:
            pheight = self.__height
            elon = self.__lon
            nch = self.__nch

        # Manual vtk generation
        f = open(name,'w')
        f.write("# vtk DataFile Version 2.0\n")
        f.write("# Map "+name+" \n")
        f.write("ASCII\n")
        f.write("DATASET STRUCTURED_GRID\n");
        f.write("DIMENSIONS {0} {1} {2} \n".format( \
                 self.__nth,nch,1))
        f.write("POINTS  {0} float\n".format(self.__nth* \
                                             nch))

        for ii,lat in zip(range(self.__nth),self.__lat):
            la = (90. + lat)*self.__dera
            for jj,lon in zip(range(nch),elon):
                lo = (lon + 180.)*self.__dera
                x = (R + pheight[ii,jj])* \
                    np.sin(la)*np.cos(lo)
                y = (R + pheight[ii,jj])* \
                    np.sin(la)*np.sin(lo)
                z = (R + pheight[ii,jj])*np.cos(la)
                f.write("{0} {1} {2} \n".format(x,y,z))

        f.write("\n")
        f.write("POINT_DATA {0}\n".format(nch*self.__nth))
        f.write("SCALARS Height float\n");
        f.write("LOOKUP_TABLE default\n");
        for ii in range(self.__nth):
            for jj in range(nch):
                 f.write("{0} ".format(pheight[ii,jj]))
        f.write("\n")
        f.close()

        if only_h:
            return

        # Wind
        if self.__wind_exist:

            if self.__fullmapx:
                plonv, elon = addcyclic(self.__wind[:,:,1], \
                                        self.__lon)
                platv, elon = addcyclic(self.__wind[:,:,0], \
                                        self.__lon)
                pmodv, elon = addcyclic(self.__wind[:,:,2], \
                                        self.__lon)
                elon[-1] = 179.99
                nch = len(elon)
            else:
                plonv = self.__wind[:,:,1]
                platv = self.__wind[:,:,0]
                pmodv = self.__wind[:,:,2]
                elon = self.__lon
                nch = len(elon)

            # Manual vtk generation
            try:

                if '.vtk' in name[-4:]:
                    namev = name[0:-4] + '-v.vtk'
                else:
                    namev = name+'-v'

            except:
                namev = name+'-v'

            f = open(namev,'w')
            f.write("# vtk DataFile Version 2.0\n")
            f.write("# Map "+name+" \n")
            f.write("ASCII\n")
            f.write("DATASET STRUCTURED_GRID\n");
            f.write("DIMENSIONS {0} {1} {2} \n".format( \
                     self.__nth,nch,1))
            f.write("POINTS  {0} float\n".format(self.__nth* \
                                                 nch))

            for ii,lat in zip(range(self.__nth),self.__lat):
                la = (90. + lat)*self.__dera
                for jj,lon in zip(range(nch),elon):
                    lo = (lon + 180.)*self.__dera
                    x = (R + pheight[ii,jj])* \
                        np.sin(la)*np.cos(lo)
                    y = (R + pheight[ii,jj])* \
                        np.sin(la)*np.sin(lo)
                    z = (R + pheight[ii,jj])*np.cos(la)
                    f.write("{0} {1} {2} \n".format(x,y,z))

            f.write("POINT_DATA {0}\n".format(nch*self.__nth))
            f.write("VECTORS v float\n");

            for ii,lat in zip(range(self.__nth),self.__lat):
                la = (90. + lat)*self.__dera
                ct = np.cos(la)
                st = np.sin(la)
                for jj,lon in zip(range(nch),elon):
                    lo = (lon + 180.)*self.__dera
                    cc = np.cos(lo)
                    sc = np.sin(lo)

                    lav = platv[ii,jj]
                    lov = plonv[ii,jj]

                    x = st*lav
                    y = lov
                    z = ct*lav

                    vx = cc*x - sc*y
                    vy = sc*x + cc*y
                    vz = z

                    f.write("{0} {1} {2}\n".format(vx,vy,vz))

            f.close()

        # Temperature
        if self.__temperature_exist:

            if self.__fullmapx:
                ptemp, elon = addcyclic(self.__temperature, \
                                        self.__lon)
                elon[-1] = 179.99
                nch = len(elon)
            else:
                ptemp = self.__temperature
                elon = self.__lon
                nch = len(elon)

            # Manual vtk generation
            try:

                if '.vtk' in name[-4:]:
                    namev = name[0:-4] + '-T.vtk'
                else:
                    namev = name+'-T'

            except:
                namev = name+'-T'

            f = open(namev,'w')
            f.write("# vtk DataFile Version 2.0\n")
            f.write("# Map "+name+" \n")
            f.write("ASCII\n")
            f.write("DATASET STRUCTURED_GRID\n");
            f.write("DIMENSIONS {0} {1} {2} \n".format( \
                     self.__nth,nch,1))
            f.write("POINTS  {0} float\n".format(self.__nth* \
                                                 nch))

            for ii,lat in zip(range(self.__nth),self.__lat):
                la = (90. + lat)*self.__dera
                for jj,lon in zip(range(nch),elon):
                    lo = (lon + 180.)*self.__dera
                    x = (R + pheight[ii,jj])* \
                        np.sin(la)*np.cos(lo)
                    y = (R + pheight[ii,jj])* \
                        np.sin(la)*np.sin(lo)
                    z = (R + pheight[ii,jj])*np.cos(la)
                    f.write("{0} {1} {2} \n".format(x,y,z))

            f.write("POINT_DATA {0}\n".format(nch*self.__nth))
            f.write("SCALARS Temperature float\n");
            f.write("LOOKUP_TABLE default\n");

            for ii in range(self.__nth):
                for jj in range(nch):
                     f.write("{0} ".format(ptemp[ii,jj]))
            f.write("\n")
            f.close()

        # Moisture
        if self.__moist_exist:

            if self.__fullmapx:
                pmoist, elon = addcyclic(self.__moist, \
                                        self.__lon)
                elon[-1] = 179.99
                nch = len(elon)
            else:
                pmoist = self.__moist
                elon = self.__lon
                nch = len(elon)

            # Manual vtk generation
            try:

                if '.vtk' in name[-4:]:
                    namev = name[0:-4] + '-h.vtk'
                else:
                    namev = name+'-h'

            except:
                namev = name+'-h'

            f = open(namev,'w')
            f.write("# vtk DataFile Version 2.0\n")
            f.write("# Map "+name+" \n")
            f.write("ASCII\n")
            f.write("DATASET STRUCTURED_GRID\n");
            f.write("DIMENSIONS {0} {1} {2} \n".format( \
                     self.__nth,nch,1))
            f.write("POINTS  {0} float\n".format(self.__nth* \
                                                 nch))

            for ii,lat in zip(range(self.__nth),self.__lat):
                la = (90. + lat)*self.__dera
                for jj,lon in zip(range(nch),elon):
                    lo = (lon + 180.)*self.__dera
                    x = (R + pheight[ii,jj])* \
                        np.sin(la)*np.cos(lo)
                    y = (R + pheight[ii,jj])* \
                        np.sin(la)*np.sin(lo)
                    z = (R + pheight[ii,jj])*np.cos(la)
                    f.write("{0} {1} {2} \n".format(x,y,z))

            f.write("POINT_DATA {0}\n".format(nch*self.__nth))
            f.write("SCALARS Moisture float\n");
            f.write("LOOKUP_TABLE default\n");

            for ii in range(self.__nth):
                for jj in range(nch):
                     f.write("{0} ".format(pmoist[ii,jj]))
            f.write("\n")
            f.close()

        # Biome
        if self.__biome_exist:

            if self.__fullmapx:
                pbiome, elon = addcyclic(self.__biome, \
                                         self.__lon)
                elon[-1] = 179.99
                nch = len(elon)
            else:
                pbiome = self.__biome
                elon = self.__lon
                nch = len(elon)

            # Manual vtk generation
            try:

                if '.vtk' in name[-4:]:
                    namev = name[0:-4] + '-b.vtk'
                else:
                    namev = name+'-b'

            except:
                namev = name+'-b'

            f = open(namev,'w')
            f.write("# vtk DataFile Version 2.0\n")
            f.write("# Map "+name+" \n")
            f.write("ASCII\n")
            f.write("DATASET STRUCTURED_GRID\n");
            f.write("DIMENSIONS {0} {1} {2} \n".format( \
                     self.__nth,nch,1))
            f.write("POINTS  {0} float\n".format(self.__nth* \
                                                 nch))

            for ii,lat in zip(range(self.__nth),self.__lat):
                la = (90. + lat)*self.__dera
                for jj,lon in zip(range(nch),elon):
                    lo = (lon + 180.)*self.__dera
                    x = (R + pheight[ii,jj])* \
                        np.sin(la)*np.cos(lo)
                    y = (R + pheight[ii,jj])* \
                        np.sin(la)*np.sin(lo)
                    z = (R + pheight[ii,jj])*np.cos(la)
                    f.write("{0} {1} {2} \n".format(x,y,z))

            f.write("POINT_DATA {0}\n".format(nch*self.__nth))
            f.write("SCALARS Biome float\n");
            f.write("LOOKUP_TABLE default\n");

            for ii in range(self.__nth):
                for jj in range(nch):
                     f.write("{0} ".format(pbiome[ii,jj]))
            f.write("\n")
            f.close()


######################################################################
######################################################################
######################################################################
######################################################################

    def save_parameters(self, name=None):
        ''' Stores the parameters in a .par file
        '''

        if name is None:
            name = self.__name + '.par'

        f = open(name,'wb')
        self.__save_parameters(f)
        f.close()

######################################################################
######################################################################
######################################################################
######################################################################

    def save_map(self, name=None):
        ''' Stores the parameters in a .map file
        '''

        if not self.__exist:
            msg = 'Must generate a map first'
            self.__error(msg)
            return

        if name is None:
            iname = self.__name + '.map'
        else:
            if isinstance(name, str):
                iname = name
            else:
                iname = self.__name + '.map'

        f = open(iname,'wb')
        self.__save_parameters(f)
        self.__save_map(f)
        f.close()

######################################################################
######################################################################
######################################################################
######################################################################

    def __save_parameters(self, f):
        ''' Really stores parameters
        '''

        f.write(struct.pack('<i', self.__nth))
        f.write(struct.pack('<i', self.__nch))
        f.write(struct.pack('<dd', *self.__thrange))
        f.write(struct.pack('<dd', *self.__chrange))
        f.write(struct.pack('<dd', *self.__pthrange))
        f.write(struct.pack('<dd', *self.__pchrange))
        f.write(struct.pack('<i', self.__seed))
        f.write(struct.pack('<i', self.__octaves))
        f.write(struct.pack('<d', self.__frequency))
        f.write(struct.pack('<d', self.__persistence))
        f.write(struct.pack('<d', self.__water))
        f.write(struct.pack('<d', self.__maxdepth))
        f.write(struct.pack('<d', self.__maxheight))
        f.write(struct.pack('<i', self.__wind_nnodes))
        if self.__wind_nodes is None:
            f.write(struct.pack('<i', 0))
        else:
            if len(self.__wind_nodes) < 1:
                f.write(struct.pack('<i', 0))
            else:
                f.write(struct.pack('<i', len(self.__wind_nodes)))
                for node in self.__wind_nodes:
                    f.write(struct.pack('<d', node[0]))
                    f.write(struct.pack('<d', node[1]))
                    f.write(struct.pack('<d', node[2]))
                    f.write(struct.pack('<d', node[3]))
        f.write(struct.pack('<d', self.__maxwindspeed))
        f.write(struct.pack('<d', self.__mintemperature))
        f.write(struct.pack('<d', self.__maxtemperature))
        mode = self.__projection + \
               ' '*np.amax([0,6 - len(list(self.__projection))])
        f.write(bytearray(mode))
        name = self.__name + \
               ' '*np.amax([0,2048 - len(list(self.__name))])
        f.write(bytearray(name))

######################################################################
######################################################################
######################################################################
######################################################################

    def __save_map(self, f, only_h=None):
        ''' Really stores map
        '''

        if only_h is None:
            ionly_h = False
        else:
            if isinstance(only_h, bool):
                ionly_h = only_h
            else:
                ionly_h = False

        form = '<'+'d'*self.__nth
        f.write(struct.pack(form, *self.__lat))
        form = '<'+'d'*self.__nch
        f.write(struct.pack(form, *self.__lon))
        form = '<'+'d'*self.__nth*self.__nch
        f.write(struct.pack(form, *(self.__height.flatten())))
        form = '<d'
        f.write(struct.pack(form, self.__shift))
        f.write(struct.pack(form, self.__minz))
        f.write(struct.pack(form, self.__maxz))

        if ionly_h:
            return

        if self.__wind_exist:
            f.write(struct.pack('<i', 1))
            form = '<'+'d'*self.__nth*self.__nch*3
            f.write(struct.pack(form, *(self.__wind.flatten())))
            form = '<d'
            f.write(struct.pack(form, self.__minwind))
            f.write(struct.pack(form, self.__maxwind))
        else:
            f.write(struct.pack('<i', -1))

        if self.__temperature_exist:
            f.write(struct.pack('<i', 1))
            form = '<'+'d'*self.__nth*self.__nch
            f.write(struct.pack(form, \
                                *(self.__temperature.flatten())))
            form = '<d'
            f.write(struct.pack(form, self.__mintemperature))
            f.write(struct.pack(form, self.__maxtemperature))
        else:
            f.write(struct.pack('<i', -1))

        if self.__moist_exist:
            f.write(struct.pack('<i', 1))
            form = '<'+'d'*self.__nth*self.__nch
            f.write(struct.pack(form, *(self.__moist.flatten())))
            form = '<d'
            f.write(struct.pack(form, self.__minmoist))
            f.write(struct.pack(form, self.__maxmoist))
        else:
            f.write(struct.pack('<i', -1))

        if self.__biome_exist:
            f.write(struct.pack('<i', 1))
            form = '<'+'d'*self.__nth*self.__nch
            f.write(struct.pack(form, *(self.__biome.flatten())))
        else:
            f.write(struct.pack('<i', -1))

        if self.__river_exist:
            nn = len(self.__river)
            f.write(struct.pack('<i', nn))
            if nn > 0:
                for river in self.__river:
                    length = len(river[0])
                    f.write(struct.pack('<i', length))
                    for il in range(length):
                        f.write(struct.pack('<d', river[0][il]))
                        f.write(struct.pack('<d', river[1][il]))
        else:
            f.write(struct.pack('<i', -1))

        if self.__lake_exist:
            nn = len(self.__lake)
            f.write(struct.pack('<i', nn))
            if nn > 0:
                for lake in self.__lake:
                    length = len(lake[0])
                    f.write(struct.pack('<i', length))
                    for il in range(length):
                        f.write(struct.pack('<d', lake[0][il]))
                        f.write(struct.pack('<d', lake[1][il]))
        else:
            f.write(struct.pack('<i', -1))

######################################################################
######################################################################
######################################################################
######################################################################

    def load_parameters(self, name=None):
        ''' Load the parameters from a .par file
        '''

        if name is None:
            name = self.__name + '.par'

        try:
            f = open(name,'rb')
            check = self.__load_parameters(f)
            f.close()
        except ValueError:
            msg = 'Invalid inputs in load_parameters'
            self.__error(msg)
        except:
            msg = 'Unexpected error in load_parameters'
            error = sys.exc_info()[:2]
            self.__error(msg, error)

######################################################################
######################################################################
######################################################################
######################################################################

    def load_map(self, name=None):
        ''' Load a .map file
        '''

        if name is None:
            name = self.__name + '.map'

        try:
            f = open(name,'rb')
            check = self.__load_parameters(f)

            if check == 0:

                checkb = self.__load_map(f)

                if checkb == 0:
                    self.__exist = True
                else:
                    if self.__exist:
                        msg = 'There is no map anymore'
                        self.__warning(msg)
                        self.__exist = False

                # Check full map
                thrange = [self.__thrange[0]*self.__dera, \
                           self.__thrange[1]*self.__dera]
                chrange = [self.__chrange[0]*self.__dera, \
                           self.__chrange[1]*self.__dera]
                if np.absolute(thrange[1] - \
                               thrange[0] - np.pi) < 1e-3:
                    self.__fullmapy = True
                else:
                    self.__fullmapy = False
                if np.absolute(chrange[1] - \
                               chrange[0] - np.pi*2.) < 1e-3:
                    step = 2.*np.pi/self.__nch
                    Tch = chrange[1] - step
                    Lon = np.linspace(chrange[0],Tch,self.__nch)
                    self.__fullmapx = True
                else:
                    Lon = np.linspace(chrange[0], \
                                      chrange[1],self.__nch)
                    self.__fullmapx = False


                cw, ct, cm, cb = self.__load_map_weather(f)

                if cw > 0:
                    self.__wind_exist = True
                else:
                    if self.__wind_exist:
                        msg = 'There is no wind map ' + \
                              'anymore'
                        self.__warning(msg)
                        self.__wind_exist = False

                if ct > 0:
                    self.__temperature_exist = True
                else:
                    if self.__temperature_exist:
                        msg = 'There is no ' + \
                              'temperature map anymore'
                        self.__warning(msg)
                        self.__temperature_exist = False

                if cm > 0:
                    self.__moist_exist = True
                else:
                    if self.__moist_exist:
                        msg = 'There is no ' + \
                              'moisture map anymore'
                        self.__warning(msg)
                        self.__moist_exist = False

                if cb > 0:
                    self.__biome_exist = True
                else:
                    if self.__biome_exist:
                        msg = 'There is no ' + \
                              'biome map anymore'
                        self.__warning(msg)
                        self.__biome_exist = False


                cr, cl = self.__load_map_water(f)

                if cr == 0:
                    self.__river_exist = True
                else:
                    if self.__river_exist:
                        msg = 'There are no rivers ' + \
                              'anymore'
                        self.__warning(msg)
                        self.__river_exist = False

                if cl == 0:
                    self.__lake_exist = True
                else:
                    if self.__lake_exist:
                        msg = 'There are no ' + \
                              'lakes anymore'
                        self.__warning(msg)
                        self.__lake_exist = False
            f.close()
        except ValueError:
            msg = 'Invalid inputs in load_map'
            self.__error(msg)
        except:
            msg = 'Unexpected error in load_map'
            error = sys.exc_info()[:2]
            self.__error(msg, error)

######################################################################
######################################################################
######################################################################
######################################################################

    def __load_parameters(self, f):
        ''' Really loads the parameters
        '''

        try:
            bit = f.read(4)
            self.__nth = struct.unpack('<i', bit)[0]
            if self.__nth < 3:
                raise ValueError()
            bit = f.read(4)
            self.__nch = struct.unpack('<i', bit)[0]
            if self.__nch < 3:
                raise ValueError()
            bit = f.read(16)
            self.__thrange = list(struct.unpack('<dd', bit)[0:2])
            l0, l1 = self.__thrange
            if l0 >= l1 or \
               l0 < -90. or \
               l0 > 90. or \
               l1 < -90. or \
               l1 > 90.:
                raise ValueError()
            bit = f.read(16)
            self.__chrange = list(struct.unpack('<dd', bit)[0:2])
            l0, l1 = self.__chrange
            if l0 >= l1 or \
               l0 < -180. or \
               l0 > 180. or \
               l1 < -180. or \
               l1 > 180.:
                raise ValueError()
            bit = f.read(16)
            self.__pthrange = list(struct.unpack('<dd', bit)[0:2])
            l0, l1 = self.__pthrange
            if l0 >= l1 or \
               l0 < -90. or \
               l0 > 90. or \
               l1 < -90. or \
               l1 > 90.:
                raise ValueError()
            bit = f.read(16)
            self.__pchrange = list(struct.unpack('<dd', bit)[0:2])
            l0, l1 = self.__pchrange
            if l0 >= l1 or \
               l0 < -180. or \
               l0 > 180. or \
               l1 < -180. or \
               l1 > 180.:
                raise ValueError()
            bit = f.read(4)
            self.__seed = struct.unpack('<i', bit)[0]
            bit = f.read(4)
            self.__octaves = struct.unpack('<i', bit)[0]
            if self.__octaves < 1:
                raise ValueError()
            bit = f.read(8)
            self.__frequency = struct.unpack('<d', bit)[0]
            if self.__frequency < 0. or self.__frequency > 1.:
                raise ValueError()
            bit = f.read(8)
            self.__persistence = struct.unpack('<d', bit)[0]
            if self.__persistence < 0.:
                raise ValueError()
            bit = f.read(8)
            self.__water = struct.unpack('<d', bit)[0]
            if self.__water >= 0.:
                if self.__water > 100.:
                    raise ValueError()
            else:
                self.__water = -1
            bit = f.read(8)
            self.__maxdepth = struct.unpack('<d', bit)[0]
            if self.__maxdepth <= 0.:
                raise ValueError()
            bit = f.read(8)
            self.__maxheight = struct.unpack('<d', bit)[0]
            if self.__maxheight <= 0.:
                raise ValueError()
            bit = f.read(4)
            self.__wind_nnodes = struct.unpack('<i', bit)[0]
            bit = f.read(4)
            rwind_nnodes = struct.unpack('<i', bit)[0]
            if rwind_nnodes > 0:
                self.__wind_nodes = []
                for inod in range(rwind_nnodes):
                    bit = f.read(8)
                    lat = struct.unpack('<d', bit)[0]
                    bit = f.read(8)
                    lon = struct.unpack('<d', bit)[0]
                    bit = f.read(8)
                    sig = struct.unpack('<d', bit)[0]
                    bit = f.read(8)
                    wei = struct.unpack('<d', bit)[0]
                    if lat < -90. or lat > 90.:
                        raise ValueError()
                    if lon < -180. or lon > 180.:
                        raise ValueError()
                    if int(np.absolute(sig)) != 1:
                        raise ValueError()
                    if wei < 0.:
                        raise ValueError()
                    self.__wind_nodes.append([lat,lon,sig,wei])
            bit = f.read(8)
            self.__maxwindspeed = struct.unpack('<d', bit)[0]
            if self.__maxwindspeed <= 0.:
                raise ValueError()
            bit = f.read(8)
            self.__mintemperature = struct.unpack('<d', bit)[0]
            if self.__mintemperature > 0.:
                raise ValueError()
            bit = f.read(8)
            self.__maxtemperature = struct.unpack('<d', bit)[0]
            if self.__maxtemperature < 0.:
                raise ValueError()
            bit = f.read(6)
            self.__projection = str(bit).strip()
            if self.__projection not in self.__modes.keys():
                raise ValueError()
            bit = f.read(2048)
            self.__name = str(bit).strip()
        except ValueError:
            msg = 'Invalid inputs loading parameters'
            self.__error(msg)
            return -1
        except:
            msg = 'Unexpected error loading parameters'
            error = sys.exc_info()[:2]
            self.__error(msg, error)
            return -1

        return 0

######################################################################
######################################################################
######################################################################
######################################################################

    def __load_map(self, f):
        ''' Really loads the map
        '''

        try:
            bit = f.read(8*self.__nth)
            form = '<'+'d'*self.__nth
            self.__lat = np.array(struct.unpack(form, bit))
            bit = f.read(8*self.__nch)
            form = '<'+'d'*self.__nch
            self.__lon = np.array(struct.unpack(form, bit))
            bit = f.read(8*self.__nth*self.__nch)
            form = '<'+'d'*self.__nth*self.__nch
            self.__height = np.array(struct.unpack(form, bit)). \
                               reshape(self.__nth,self.__nch)
            bit = f.read(8)
            form = '<d'
            self.__shift = float(struct.unpack(form, bit)[0])
            bit = f.read(8)
            form = '<d'
            self.__minheight = float(struct.unpack(form, bit)[0])
            bit = f.read(8)
            form = '<d'
            self.__maxheight = float(struct.unpack(form, bit)[0])
            return 0
        except:
            msg = 'Unexpected error loading map'
            error = sys.exc_info()[:2]
            self.__error(msg, error)
            return -1

######################################################################
######################################################################
######################################################################
######################################################################

    def __load_map_weather(self, f):
        ''' Really loads the map weather
        '''

        try:
            bit = f.read(4)
            cw = struct.unpack('<i', bit)
            if cw > 0:
                bit = f.read(8*self.__nth*self.__nch*3)
                form = '<'+'d'*self.__nth*self.__nch*3
                self.__wind = np.array(struct.unpack(form, bit)). \
                                     reshape(self.__nth,self.__nch,3)
                bit = f.read(8)
                form = '<d'
                self.__minwind = float(struct.unpack(form, bit)[0])
                bit = f.read(8)
                form = '<d'
                self.__maxwind = float(struct.unpack(form, bit)[0])
        except:
            cw = -1

        try:
            bit = f.read(4)
            ct = struct.unpack('<i', bit)
            if ct > 0:
                bit = f.read(8*self.__nth*self.__nch)
                form = '<'+'d'*self.__nth*self.__nch
                self.__temperature = \
                               np.array(struct.unpack(form, bit)). \
                                        reshape(self.__nth,self.__nch)
                bit = f.read(8)
                form = '<d'
                self.__mintemperature = \
                                    float(struct.unpack(form, bit)[0])
                bit = f.read(8)
                form = '<d'
                self.__maxtemperature = \
                                    float(struct.unpack(form, bit)[0])
        except:
            ct = -1

        try:
            bit = f.read(4)
            cm = struct.unpack('<i', bit)
            if cm > 0:
                bit = f.read(8*self.__nth*self.__nch)
                form = '<'+'d'*self.__nth*self.__nch
                self.__moist = np.array(struct.unpack(form, bit)). \
                                  reshape(self.__nth,self.__nch)
                bit = f.read(8)
                form = '<d'
                self.__minmoist = float(struct.unpack(form, bit)[0])
                bit = f.read(8)
                form = '<d'
                self.__maxmoist = float(struct.unpack(form, bit)[0])
        except:
            cm = -1

        try:
            bit = f.read(4)
            cb = struct.unpack('<i', bit)
            if cb > 0:
                bit = f.read(8*self.__nth*self.__nch)
                form = '<'+'d'*self.__nth*self.__nch
                self.__biome = np.array(struct.unpack(form, bit)). \
                                  reshape(self.__nth,self.__nch)
        except:
            cb = -1

        return cw, ct, cm, cb

######################################################################
######################################################################
######################################################################
######################################################################

    def __load_map_water(self, f):
        ''' Really loads the map rivers and lakes
        '''

        # Rivers
        try:
            bit = f.read(4)
            form = '<i'
            nrivers = struct.unpack(form,bit)[0]
            if nrivers < 1:
                riv = -1
            else:
                self.__river = []
                riv = 0
                for ir in range(nrivers):
                    bit = f.read(4)
                    form = '<i'
                    length = struct.unpack(form,bit)[0]
                    riverx = []
                    rivery = []
                    for il in range(length):
                        form = '<d'
                        bit = f.read(8)
                        riverx.append(struct.unpack(form,bit)[0])
                        bit = f.read(8)
                        rivery.append(struct.unpack(form,bit)[0])
                    self.__river.append([riverx,rivery])
        except:
            riv = -1

        # Lakes
        try:
            bit = f.read(4)
            form = '<i'
            nlakes = struct.unpack(form,bit)[0]
            if nlakes < 1:
                lak = -1
            else:
                self.__lake = []
                lak = 0
                for ir in range(nlakes):
                    bit = f.read(4)
                    form = '<i'
                    length = struct.unpack(form,bit)[0]
                    lakex = []
                    lakey = []
                    for il in range(length):
                        form = '<d'
                        bit = f.read(8)
                        lakex.append(struct.unpack(form,bit)[0])
                        bit = f.read(8)
                        lakey.append(struct.unpack(form,bit)[0])
                    self.__lake.append([lakex,lakey])
        except:
            lak = -1

        return riv, lak

######################################################################
######################################################################
######################################################################
######################################################################

    def refine(self, nth=None, nch=None, thrange=None, chrange=None, \
               force=None):
        ''' Rotates the longitude axis the quantity specified
        '''

        # Check there is a map
        if not self.__exist:
            msg = 'Must generate a map first'
            self.__error(msg)
            return

        # Check there is interpolation
        if not self.__interp:
            if force is None:
                iforce = False
            else:
                if isinstance(force, bool):
                    iforce = force
                else:
                    iforce = False
            if not iforce:
                msg = 'Without interpolation, the weather would ' + \
                      'be destroyed. Use force=True to do it'
                self.__error(msg)
                return

        # Polar nodes
        if nth is None:
            msg = 'nth is a required parameter'
            self.__error(msg)
            return
        else:
            try:
                __nth = int(nth)
                if __nth < 3:
                    msg = 'nth must be at least 3.'
                    self.__error(msg)
                    return
            except ValueError:
                msg = 'ValueError in nth.'
                self.__error(msg)
                return
            except:
                msg = 'Unexpected error with nth'
                error = sys.exc_info()[:2]
                self.__error(msg, error)
                return

        # Azimuthal nodes
        if nth is None:
            msg = 'nch is a required parameter'
            self.__error(msg)
            return
        else:
            try:
                __nch = int(nch)
                if __nch < 3:
                    msg = 'nch must be at least 3.'
                    self.__error(msg)
                    return
            except ValueError:
                msg = 'ValueError in nch.'
                self.__error(msg)
                return
            except:
                msg = 'Unexpected error with nch'
                error = sys.exc_info()[:2]
                self.__error(msg, error)
                return

        # Polar range
        if thrange is None:
            msg = 'thrange is a required parameter'
            self.__error(msg)
            return
        else:
            if len(thrange) == 2:
                try:
                    __thrange = []
                    for th in thrange:
                        __thrange.append(float(th))
                    if __thrange[0] >= __thrange[1] or \
                       __thrange[0] < -90. or \
                       __thrange[0] > 90. or \
                       __thrange[1] < -90. or \
                       __thrange[1] > 90.:
                        msg = 'Bad range in thrange.'
                        self.__error(msg)
                        return
                except ValueError:
                    msg = 'ValueError in thrange.'
                    self.__error(msg)
                    return
                except:
                    msg = 'Unexpected error with thrange'
                    error = sys.exc_info()[:2]
                    self.__error(msg, error)
                    return
            else:
                msg = 'Wrong format thrange. Set default'
                self.__warning(msg)
                return


        # Azimuthal range
        if chrange is None:
            msg = 'chrange is a required parameter'
            self.__error(msg)
            return
        else:
            if len(chrange) == 2:
                try:
                    __chrange = []
                    for ch in chrange:
                        __chrange.append(float(ch))
                    if __chrange[0] >= __chrange[1] or \
                       __chrange[0] < -180. or \
                       __chrange[0] > 180. or \
                       __chrange[1] < -180. or \
                       __chrange[1] > 180.:
                        msg = 'Bad range in chrange.'
                        self.__error(msg)
                        return
                except ValueError:
                    msg = 'ValueError in chrange.'
                    self.__error(msg)
                    return
                except:
                    msg = 'Unexpected error with chrange'
                    error = sys.exc_info()[:2]
                    self.__error(msg, error)
                    return
            else:
                msg = 'Wrong format chrange. Set default'
                self.__error(msg)
                return

        # Check that you changed anything
        if __chrange[0] == self.__chrange[0] and \
           __chrange[1] == self.__chrange[1] and \
           __thrange[0] == self.__thrange[0] and \
           __thrange[1] == self.__thrange[1] and \
           __nth == self.__nth and __nch == self.__nch:
            msg = 'You have not changed the ranges at all'
            self.__error(msg)
            return

        # Store old quantities
        __lat_old = copy.deepcopy(self.__lat)
        __lon_old = copy.deepcopy(self.__lon)
        __nth_old = copy.deepcopy(self.__nth)
        __nch_old = copy.deepcopy(self.__nch)
        __thrange_old = copy.deepcopy(self.__thrange)
        __chrange_old = copy.deepcopy(self.__chrange)
        __height_old = copy.deepcopy(self.__height)
        __temperature_old = copy.deepcopy(self.__temperature)

        # Update to new quantities
        self.__nth = __nth
        self.__nch = __nch
        self.__thrange = __thrange
        self.__chrange = __chrange

        # Check plot ranges
        if self.__pthrange[0] < self.__thrange[0]:
            self.__pthrange[0] = self.__thrange[0]
            msg = 'Adjusted lower limit of pthrange'
           #self.__warning(msg)
        if self.__pthrange[1] > self.__thrange[1]:
            self.__pthrange[1] = self.__thrange[1]
            msg = 'Adjusted upper limit of pthrange'
           #self.__warning(msg)
        if self.__pchrange[0] < self.__chrange[0]:
            self.__pchrange[0] = self.__chrange[0]
            msg = 'Adjusted lower limit of pchrange'
           #self.__warning(msg)
        if self.__pchrange[1] > self.__chrange[1]:
            self.__pchrange[1] = self.__chrange[1]
            msg = 'Adjusted upper limit of pchrange'
           #self.__warning(msg)

        # Print message
        msg = 'Refining the map...'
        self.__print(msg)

        # Do map new
        self.generate_map(silent=True, refine=True)

        # Check if failed
        if not self.__exist:
            msg = 'The map failed, destroying everything'
            self.__error(msg)
            self.__wind_exist = False
            self.__temperature_exist = False
            self.__moist_exist = False
            self.__biome_exist = False
            self.__river_exist = False
            self.__lake_exist = False

        # Do the weather

        # Wind
        if self.__wind_exist:

            check = self.__redo_wind()

            if not check:

                w0 = copy.deepcopy(self.__wind[:,:,0])
                w1 = copy.deepcopy(self.__wind[:,:,1])
                w2 = copy.deepcopy(self.__wind[:,:,2])

                self.__interN(w0, __lat_old, __lon_old, __height_old)
                self.__interN(w1, __lat_old, __lon_old, __height_old)
                self.__interN(w2, __lat_old, __lon_old, __height_old)

                self.__wind = np.zeros((self.__nth,self.__nch,3))
                self.__wind[:,:,0] = w0
                self.__wind[:,:,1] = w1
                self.__wind[:,:,2] = w2

        # Temperature
        if self.__temperature_exist:

            self.__temperature = self.__interN(self.__temperature, \
                                               __lat_old, __lon_old, \
                                               __height_old)

        # Moisture
        if self.__moist_exist:

            self.__moist = self.__interN(self.__moist, \
                                         __lat_old, __lon_old, \
                                         __height_old, \
                                         __temperature_old)

        # Biome
        self.generate_biome(silent=True)

######################################################################
######################################################################

    def __interN(self, data, xold, yold, zold, wold=None):
        ''' Pseudo-interpolation to a different grid. It looks for
            the closest neighbour
        '''

        # Weights for distances
        wh = 1.
        wz = 2.
        ww = 1.

        # Limit the old maps to a 50% extra of the refined

        # Lat
        dlat = self.__thrange[1] - self.__thrange[0]
        dlat *= 0.25
        lat0 = self.__thrange[0] - dlat
        lat1 = self.__thrange[1] + dlat
        lat0 = np.max([-90., lat0])
        lat1 = np.min([ 90., lat1])

        # Lon
        dlon = self.__chrange[1] - self.__chrange[0]
        dlon *= 0.25
        lon0 = self.__chrange[0] - dlon
        lon1 = self.__chrange[1] + dlon
        lon0 = np.max([-180., lon0])
        lon1 = np.min([ 179., lon1])

        # Limit the old arrays
        ix0 = np.argmin(np.absolute(xold - lat0))
        ix1 = np.argmin(np.absolute(xold - lat1))
        __xold = xold[ix0:ix1]
        iy0 = np.argmin(np.absolute(yold - lon0))
        iy1 = np.argmin(np.absolute(yold - lon1))
        __yold = yold[iy0:iy1]
        __zold = zold[ix0:ix1,iy0:iy1]
        if wold is not None:
            __wold = wold[ix0:ix1,iy0:iy1]
        data = data[ix0:ix1,iy0:iy1]

        # Maximums
        maxx = lat1 - lat0
        maxx = 1./maxx
        maxy = lon1 - lon0
        maxy = 1./maxy
        maxz = np.max(__zold)
        wz /= maxz
        if wold is not None:
            maxw = np.max(self.__temperature) - \
                   np.min(self.__temperature)
            ww /= maxw

        # Allocate new data
        datanew = np.zeros((self.__nth, self.__nch))

        # For each new point
        for ilat, xnew in zip(range(self.__nth), self.__lat):
            for ilon, ynew in zip(range(self.__nch), self.__lon):

                znew = self.__height[ilat, ilon]

                dx = (__xold - xnew)*maxx
                dx *= dx
                dy = (__yold - ynew)*maxy
                dy *= dy
                ddx, ddy = np.meshgrid(dx, dy, indexing='ij')
                dz = __zold - znew
                dz *= dz

                if wold is None:
                    dw = np.zeros(dz.shape)
                else:
                    wnew = self.__temperature[ilat,ilon]
                    dw = __wold - wnew
                    dw *= dw

                Delta = wh*(ddx + ddy) + wz*dz + ww*dw

                inds = np.argmin(Delta)

                ix, iy = np.unravel_index(inds, dz.shape)

                datanew[ilat,ilon] = data[ix,iy]

        # Return
        return datanew

######################################################################
######################################################################

    def __interpol3D(self, data, xold, yold):
        ''' Interpolates 3D array with third coordinate length 3
        '''

        f0 = interpolate.RegularGridInterpolator((xold, yold), \
                                                 data[:,:,0])
        f1 = interpolate.RegularGridInterpolator((xold, yold), \
                                                 data[:,:,1])
        f2 = interpolate.RegularGridInterpolator((xold, yold), \
                                                 data[:,:,2])

        data = np.zeros((self.__nth, self.__nch, 3))

        x, y = np.meshgrid(self.__lat, self.__lon, indexing='ij')
        data[:,:,0] = f0((x,y))
        data[:,:,1] = f1((x,y))
        data[:,:,2] = f2((x,y))

        return data

######################################################################
######################################################################

    def __interpol2D(self, data, xold, yold):
        ''' Interpolates 2D array
        '''

        f = interpolate.RegularGridInterpolator((xold, yold), \
                                                data)

        data = np.zeros((self.__nth, self.__nch))

        x, y = np.meshgrid(self.__lat, self.__lon, indexing='ij')
        data = f((x,y))

        return data

######################################################################
######################################################################
######################################################################
######################################################################

    def rotate(self, ang=None):
        ''' Rotates the longitude axis the quantity specified
        '''

        if not self.__exist:
            msg = 'Must generate a map first'
            self.__error(msg)
            return

        if ang is None or ang < -180 or ang > 180:
            msg = 'Must specify an angle between -180 ' + \
                  'and 180 in degrees'
            self.__error(msg)
            return

        for ii in range(len(self.__lon)):

            self.__lon[ii] += ang

            if self.__lon[ii] < -180.:
                self.__lon[ii] += 359.

            if self.__lon[ii] > 179.:
                self.__lon[ii] -= 359.

        ind = self.__lon.argsort()
        self.__lon.sort()

        # Rotate height
        arraycopy = np.copy(self.__height)

        for jj in range(len(self.__lat)):
            for ii in range(len(self.__lon)):
                self.__height[jj,ii] = arraycopy[jj,ind[ii]]

        # Rotate wind
        if self.__wind_exist:
            arraycopy = np.copy(self.__wind)

            for kk in range(2):
                for jj in range(len(self.__lat)):
                    for ii in range(len(self.__lon)):
                        self.__wind[jj,ii,kk] = \
                                              arraycopy[jj,ind[ii],kk]

        # Rotate temperature
        if self.__temperature_exist:
            arraycopy = np.copy(self.__temperature)

            for jj in range(len(self.__lat)):
                for ii in range(len(self.__lon)):
                    self.__temperature[jj,ii] = arraycopy[jj,ind[ii]]

        # Rotate moist
        if self.__moist_exist:
            arraycopy = np.copy(self.__moist)

            for jj in range(len(self.__lat)):
                for ii in range(len(self.__lon)):
                    self.__moist[jj,ii] = arraycopy[jj,ind[ii]]

        # Rotate biome
        if self.__biome_exist:
            arraycopy = np.copy(self.__biome)

            for jj in range(len(self.__lat)):
                for ii in range(len(self.__lon)):
                    self.__biome[jj,ii] = arraycopy[jj,ind[ii]]

        # Rotate lakes
        if self.__lake_exist:
            for il in range(len(self.__lake)):
                for ix in range(len(self.__lake[il][1])):
                    newx = self.__lake[il][1][ix] + ang
                    if newx < -180.:
                        newx += 359.
                    if newx > 179.:
                        newx -= 359.
                    self.__lake[il][1][ix] = newx

        # Rotate rivers
        if self.__river_exist:
            for ir in range(len(self.__river)):
                for ix in range(len(self.__river[ir][1])):
                    newx = self.__river[ir][1][ix] + ang
                    if newx < -180.:
                        newx += 359.
                    if newx > 179.:
                        newx -= 359.
                    self.__river[ir][1][ix] = newx

######################################################################
######################################################################
######################################################################
######################################################################

    def __color_def(self, minv=None, maxv=None):
        ''' Defines the color bar given the limits of the map and the
            predefined height color scale
        '''

        if minv is None:
            minv = -10.0
        if maxv is None:
            maxv = 20.0

        bit_rgb = np.linspace(0,1,256)

        # Predefined color scale
        precolo = [(  0, 51,102),(  0,102,204),( 51,153,255), \
                   (102,178,255),(204,255,255),(220,255,255), \
                   (220,255,255),(  0,153,  0),(  0,204,  0), \
                   (  0,204,  0),(  0,204,102),(  0,255,128), \
                   (225,225,123),(215,215,111),(255,178,102), \
                   (204,102,  0),(153, 76,  0),(102, 51,  0), \
                   (102,  0,204),(153, 51,255),(178,102,255), \
                   (204,153,255),(229,204,255)]
        preleve = [-10.0,-8.0,-6.0, \
                   -4.0,-2.0,-0.2, \
                   -0.001,0.0,0.001, \
                   0.25,0.5,0.75, \
                   1.0,1.5,1.75, \
                   3.0,3.5,6.0, \
                   6.001,9.0,10.0, \
                   15.0,20.0]

        # Introduce the extremes of the data
        level = [minv]
        if minv >= preleve[0]:
            if minv == preleve[0]:
                color = [precolo[0]]
                ini=1
            else:
                ini = -1
                for i in range(len(preleve)-1):
                    if minv >= preleve[i] and minv < preleve[i+1]:
                        R = [precolo[i][0],precolo[i+1][0]]
                        G = [precolo[i][1],precolo[i+1][1]]
                        B = [precolo[i][2],precolo[i+1][2]]
                        color = [(np.interp([minv],[preleve[i], \
                                            preleve[i+1]],R)[0],
                                  np.interp([minv],[preleve[i], \
                                            preleve[i+1]],G)[0],
                                  np.interp([minv],[preleve[i], \
                                            preleve[i+1]],B)[0])]
                        ini = i+1
                        break
        else:
            color = [(  0,  0,  0)]
            ini = 0
        if maxv < preleve[len(preleve)-1]:
            for j in range(len(preleve[ini:])-1):
                if maxv > preleve[ini+j] and maxv <= preleve[ini+j+1]:
                    for k in range(len(preleve[ini:ini+j+1])):
                        level.append(preleve[k+ini])
                        color.append(precolo[k+ini])
                    level.append(maxv)
                    fin = ini+j
                    R = [precolo[fin][0],precolo[fin+1][0]]
                    G = [precolo[fin][1],precolo[fin+1][1]]
                    B = [precolo[fin][2],precolo[fin+1][2]]
                    color.append((np.interp([maxv],[preleve[fin], \
                                            preleve[fin+1]],R)[0], \
                                  np.interp([maxv],[preleve[fin], \
                                            preleve[fin+1]],G)[0],
                                  np.interp([maxv],[preleve[fin], \
                                            preleve[fin+1]],B)[0]))
        else:
            for j in range(len(preleve[ini:])):
                level.append(preleve[ini+j])
                color.append(precolo[ini+j])
                fin = ini+j
            if maxv != preleve[len(preleve)-1]:
                level.append(maxv)
                color.append((255,255,255))

        # Normalize the height scale
        a = 1./(maxv - minv)
        b = minv/(minv - maxv)
        for i in range(len(level)):
            level[i] = a*level[i] + b
        level[0] = 0
        level[len(level)-1]=1

        # Transform from bits to 0-1
        for i in range(len(color)):
            color[i] = (bit_rgb[int(round(color[i][0]))],
                        bit_rgb[int(round(color[i][1]))],
                        bit_rgb[int(round(color[i][2]))])

        # Generate the dictionary
        cdict = {'red':[], 'green':[], 'blue':[]}
        for lv, col in zip(level, color):
            cdict[ 'red' ].append((lv, col[0], col[0]))
            cdict['green'].append((lv, col[1], col[1]))
            cdict[ 'blue'].append((lv, col[2], col[2]))

        # Define map and return
        cmap = LinearSegmentedColormap('physicheigh',cdict,256)

        return cmap

######################################################################
######################################################################
######################################################################
######################################################################

    def __color_def_T(self, minv=None, maxv=None):
        ''' Defines the color bar given the limits of the map and the
            predefined temperature color scale
        '''

        if minv is None:
            minv = -10.0
        if maxv is None:
            maxv = 20.0

        bit_rgb = np.linspace(0,1,256)

        # Predefined color scale
        precolo = [(  0,  0,255), \
                   (100,100,255), \
                   (255,255,255), \
                   (255,200,  0), \
                   (255,100,100), \
                   (255,  0,  0), \
                   (255,  0,  0)]
        preleve = [-10, \
                    -5, \
                     0, \
                     5, \
                    20, \
                    30, \
                    100]

        # Introduce the extremes of the data
        level = [minv]
        if minv >= preleve[0]:
            if minv == preleve[0]:
                color = [precolo[0]]
                ini=1
            else:
                ini = -1
                for i in range(len(preleve)-1):
                    if minv >= preleve[i] and minv < preleve[i+1]:
                        R = [precolo[i][0],precolo[i+1][0]]
                        G = [precolo[i][1],precolo[i+1][1]]
                        B = [precolo[i][2],precolo[i+1][2]]
                        color = [(np.interp([minv],[preleve[i], \
                                            preleve[i+1]],R)[0],
                                  np.interp([minv],[preleve[i], \
                                            preleve[i+1]],G)[0],
                                  np.interp([minv],[preleve[i], \
                                            preleve[i+1]],B)[0])]
                        ini = i+1
                        break
        else:
            color = [(  0,  0,  0)]
            ini = 0
        if maxv < preleve[len(preleve)-1]:
            for j in range(len(preleve[ini:])-1):
                if maxv > preleve[ini+j] and maxv <= preleve[ini+j+1]:
                    for k in range(len(preleve[ini:ini+j+1])):
                        level.append(preleve[k+ini])
                        color.append(precolo[k+ini])
                    level.append(maxv)
                    fin = ini+j
                    R = [precolo[fin][0],precolo[fin+1][0]]
                    G = [precolo[fin][1],precolo[fin+1][1]]
                    B = [precolo[fin][2],precolo[fin+1][2]]
                    color.append((np.interp([maxv],[preleve[fin], \
                                            preleve[fin+1]],R)[0], \
                                  np.interp([maxv],[preleve[fin], \
                                            preleve[fin+1]],G)[0],
                                  np.interp([maxv],[preleve[fin], \
                                            preleve[fin+1]],B)[0]))
        else:
            for j in range(len(preleve[ini:])):
                level.append(preleve[ini+j])
                color.append(precolo[ini+j])
                fin = ini+j
            if maxv != preleve[len(preleve)-1]:
                level.append(maxv)
                color.append((255,255,255))

        # Normalize the height scale
        a = 1./(maxv - minv)
        b = minv/(minv - maxv)
        for i in range(len(level)):
            level[i] = a*level[i] + b
        level[0] = 0
        level[len(level)-1]=1

        # Transform from bits to 0-1
        for i in range(len(color)):
            color[i] = (bit_rgb[int(round(color[i][0]))],
                        bit_rgb[int(round(color[i][1]))],
                        bit_rgb[int(round(color[i][2]))])

        # Generate the dictionary    
        cdict = {'red':[], 'green':[], 'blue':[]}
        for lv, col in zip(level, color):
            cdict[ 'red' ].append((lv, col[0], col[0]))
            cdict['green'].append((lv, col[1], col[1]))
            cdict[ 'blue'].append((lv, col[2], col[2]))

        # Define map and return    
        cmap = LinearSegmentedColormap('physicheigh',cdict,256)    

        return cmap

######################################################################
######################################################################
######################################################################
######################################################################

    def __color_adjust(self, cmap, minv, maxv, ominv, omaxv):
        ''' Redefines a new color bar from a previous one respecting
            the limits
        '''

        bit_rgb = np.linspace(0,1,256)

        # Color map
        precolo = []
        preleve = []
        for i in range(256):
            precolo.append((cmap(i)[0]*255., \
                            cmap(i)[1]*255., \
                            cmap(i)[2]*255.))
            preleve.append((omaxv-ominv)*float(i)/255. + ominv)

        # Introduce the extremes of the data
        level = [minv]
        if minv >= preleve[0]:
            if minv == preleve[0]:
                color = [precolo[0]]
                ini=1
            else:
                ini = -1
                for i in range(len(preleve)-1):
                    if minv >= preleve[i] and minv < preleve[i+1]:
                        R = [precolo[i][0],precolo[i+1][0]]
                        G = [precolo[i][1],precolo[i+1][1]]
                        B = [precolo[i][2],precolo[i+1][2]]
                        color = [(np.interp([minv],[preleve[i], \
                                            preleve[i+1]],R)[0],
                                  np.interp([minv],[preleve[i], \
                                            preleve[i+1]],G)[0],
                                  np.interp([minv],[preleve[i], \
                                            preleve[i+1]],B)[0])]
                        ini = i+1
                        break
        else:
            color = [(  0,  0,  0)]
            ini = 0
        if maxv < preleve[-1]:
            for j in range(len(preleve[ini:])-1):
                if maxv > preleve[ini+j] and maxv <= preleve[ini+j+1]:
                    for k in range(len(preleve[ini:ini+j+1])):
                        level.append(preleve[k+ini])
                        color.append(precolo[k+ini])
                    level.append(maxv)
                    fin = ini+j
                    R = [precolo[fin][0],precolo[fin+1][0]]
                    G = [precolo[fin][1],precolo[fin+1][1]]
                    B = [precolo[fin][2],precolo[fin+1][2]]
                    color.append((np.interp([maxv],[preleve[fin], \
                                            preleve[fin+1]],R)[0], \
                                  np.interp([maxv],[preleve[fin], \
                                            preleve[fin+1]],G)[0],
                                  np.interp([maxv],[preleve[fin], \
                                            preleve[fin+1]],B)[0]))
        else:
            for j in range(len(preleve[ini:])):
                level.append(preleve[ini+j])
                color.append(precolo[ini+j])
                fin = ini+j
            if maxv != preleve[len(preleve)-1]:
                level.append(maxv)
                color.append((255,255,255))

        # Normalize the height scale
        a = 1./(maxv - minv)
        b = minv/(minv - maxv)
        for i in range(len(level)):
            level[i] = a*level[i] + b
        level[0] = 0
        level[len(level)-1]=1

        # Transform from bits to 0-1
        for i in range(len(color)):
            color[i] = (bit_rgb[int(round(color[i][0]))],
                        bit_rgb[int(round(color[i][1]))],
                        bit_rgb[int(round(color[i][2]))])

        # Generate the dictionary    
        cdict = {'red':[], 'green':[], 'blue':[]}
        for lv, col in zip(level, color):
            cdict[ 'red' ].append((lv, col[0], col[0]))
            cdict['green'].append((lv, col[1], col[1]))
            cdict[ 'blue'].append((lv, col[2], col[2]))

        # Define map and return    
        cmap = LinearSegmentedColormap('physicheigh',cdict,256)    

        return cmap

######################################################################
######################################################################
######################################################################
######################################################################

    def __color_def_B(self):
        ''' Defines the color bar for biomes
        '''

        minv = -3.
        maxv = 8.

        bit_rgb = np.linspace(0,1,256)

        # Predefined color scale
        color = [(  0,  0,204), \
                 (255,255,255), \
                 (102, 51,  0), \
                 (204,255,255), \
                 ( 51,204,255), \
                 (  0,153, 51), \
                 (  0,102,  0), \
                 (  0,204,202), \
                 (204,255, 51), \
                 (  0,204,102), \
                 (  0,153, 51), \
                 (255,204,  0)]
        level = [-3.0000, \
                 -2.0000, \
                 -1.0000, \
                  0.0000, \
                  1.0000, \
                  2.0000, \
                  3.0000, \
                  4.0000, \
                  5.0000, \
                  6.0000, \
                  7.0000, \
                  8.0000]

        # Normalize the height scale
        a = 1./(maxv - minv)
        b = minv/(minv - maxv)
        for i in range(len(level)):
            level[i] = a*level[i] + b
        level[0] = 0
        level[len(level)-1]=1

        # Transform from bits to 0-1
        for i in range(len(color)):
            color[i] = (bit_rgb[int(round(color[i][0]))],
                        bit_rgb[int(round(color[i][1]))],
                        bit_rgb[int(round(color[i][2]))])

        # Generate the dictionary    
        cdict = {'red':[], 'green':[], 'blue':[]}
        for lv, col in zip(level, color):
            cdict[ 'red' ].append((lv, col[0], col[0]))
            cdict['green'].append((lv, col[1], col[1]))
            cdict[ 'blue'].append((lv, col[2], col[2]))

        # Define map and return    
        cmap = LinearSegmentedColormap('biomes',cdict,12)

        return cmap

######################################################################
######################################################################
######################################################################
######################################################################

    def draw(self,resolution=None,lat_0=None, \
             lon_0=None,width=None,height=None,lat_ts=None, \
             lat_1=None,lat_2=None,lon_1=None,lon_2=None, \
             k_0=None,no_rot=None,boundinglat=None,round=None, \
             satellite_height=None,latdel=None,londel=None, \
             lat0=None,lat1=None,lon0=None,lon1=None,h0=None, \
             h1=None, fileout=None, dpi=None):
        ''' Visualize a map in 2D
        '''

        if not self.__can_plot_2D:
            msg = 'pyplot/basemap not available'
            self.__error(msg)
            return

        if not self.__exist:
            msg = 'Must generate a map first'
            self.__error(msg)
            return

        if lat0 is None:
            lat0 = self.__pthrange[0]
        if lat1 is None:
            lat1 = self.__pthrange[1]
        pthrange = [lat0,lat1]
        if lon0 is None:
            lon0 = self.__pchrange[0]
        if lon1 is None:
            lon1 = self.__pchrange[1]
        pchrange = [lon0,lon1]
        if h0 is None:
            h0 = np.min(self.__height)
        if h1 is None:
            h1 = np.max(self.__height)
        phrange = [h0,h1]
        if lat_0 is None:
            lat_0 = .5*(lat0+lat1)
        if lon_0 is None:
            lon_0 = .5*(lon0+lon1)
        if latdel is None:
            latdel = 30
        if londel is None:
            londel = 30

        # If Full map, cycle
        if self.__fullmapx:
            pheight, lon = addcyclic(self.__height,self.__lon)
            lat = self.__lat
        else:
            pheight, lon, lat = self.__height, self.__lon, self.__lat

        mapp = Basemap(projection=self.__projection, \
                       lat_0=lat_0, \
                       lon_0=lon_0, \
                       resolution=resolution, \
                       llcrnrlon=pchrange[0], \
                       urcrnrlon=pchrange[1], \
                       llcrnrlat=pthrange[0], \
                       urcrnrlat=pthrange[1], \
                       width=width, \
                       height=height, \
                       lat_ts=lat_ts, \
                       lat_1=lat_1, \
                       lat_2=lat_2, \
                       lon_1=lon_1, \
                       lon_2=lon_2, \
                       k_0=k_0, \
                       no_rot=no_rot, \
                       boundinglat=boundinglat, \
                       round=round, \
                       satellite_height=satellite_height)

        # Compute native map projection coordinates of lat/lon grid.
        x, y = mapp(*np.meshgrid(lon, lat))

        # Load color for terrain
        cco =  self.__color_def(minv=phrange[0], maxv=phrange[1])

        # Trues
        trues = [1]
        if self.__wind_exist:
            trues.append(1)
        else:
            trues.append(0)
        if self.__temperature_exist:
            trues.append(1)
        else:
            trues.append(0)
        if self.__moist_exist:
            trues.append(1)
        else:
            trues.append(0)
        if self.__biome_exist:
            trues.append(1)
        else:
            trues.append(0)

        # Subplots
        nplots = sum(trues)
        rows = [1, 1, 2, 2, 3, 3, 4, 4, 5]
        cols = [1, 2, 2, 2, 2, 2, 2, 2, 2]
        irow = [0, 0, 1, 1, 2, 2, 3, 3, 4]
        icol = [0, 1, 0, 1, 0, 1, 0, 1, 0]
        nrow = rows[nplots-1]
        ncol = cols[nplots-1]
        fsize = [(17,4),(17,4),(20,7),(20,7),(20,9),(20,9), \
                 (20,12),(20,12),(20,15)]

        fig = plt.figure(figsize=fsize[nplots-1])
        gs = gridspec.GridSpec(nrow, ncol, \
                               wspace=0.01, hspace=0.01)
        fig.clf()

        kk = -1
        # For each variable
        for ii in range(5):

            if trues[ii] < .5:
                continue

            kk += 1

            ax = plt.subplot(gs[irow[ii], icol[ii]])

            # Draw lat/lon grid lines
            mapp.drawmeridians(np.arange(pchrange[0], \
                                         pchrange[1], \
                                         londel))
            mapp.drawparallels(np.arange(pthrange[0], \
                                         pthrange[1], \
                                         latdel))

            if ii == 0:

                # Contour data over the map.
                im = ax.imshow(pheight)
                cs = mapp.contourf(x,y,pheight,255,cmap=cco)

                # Color bar
                ticks = []
                minc = np.min(pheight)
                maxc = np.max(pheight)
                mini = int(minc)-2
                maxi = int(maxc)+2
                for ff in range(mini,maxi):
                    if ff > minc:
                        ticks.append(ff)
                    elif ff > maxc:
                        break
                divider = make_axes_locatable(ax)
                cax = divider.append_axes("left", size="4%", pad=0.05)
                cbar = plt.colorbar(cs,orientation='vertical', \
                                    cax=cax,ticks=ticks)
                cax.yaxis.set_ticks_position('left')
                cbar.set_label('km')
                cax.yaxis.set_label_position('left')

                # If there are lakes
                if self.__lake_exist:
                    jj = -1
                    # For each lake
                    for lake in self.__lake:
                        jj += 1
                        lakex = lake[1]
                        lakey = lake[0]

                        # Borders
                        xl = np.min(lakex)
                        xr = np.max(lakex)
                        yl = np.min(lakey)
                        yr = np.max(lakey)

                        # Check if completely out of bounds
                        if xr < pchrange[0] or \
                           xl > pchrange[1] or \
                           yr < pthrange[0] or \
                           yl > pthrange[1]:
                            continue

                        # If the lake crosses the globe
                        if np.min(lakex) < -170. and \
                           np.max(lakex) > 170.:
                            lakexl = []
                            lakeyl = []
                            lakexr = []
                            lakeyr = []
                            for ii in range(len(lakex)):
                                if lakex[ii] > 0.:
                                    lakexl += [lakex[ii]]
                                    lakeyl += [lakey[ii]]
                                else:
                                    lakexr += [lakex[ii]]
                                    lakeyr += [lakey[ii]]
                            xla,yla = mapp(lakexr,lakeyr)
                            ax.fill(xla, yla, \
                                    color=(0.86,1.0,1.0,1.0))
                            xla,yla = mapp(lakexl,lakeyl)
                            ax.fill(xla, yla, \
                                    color=(0.86,1.0,1.0,1.0))
                        # If the lake is fine in the periodic boundary
                        else:

                            # Check that the lake is not split
                            if xl < pchrange[0] or \
                               xr > pchrange[1] or \
                               yl < pthrange[0] or \
                               yr > pthrange[1]:

                                lakexl = []
                                lakeyl = []
                                for ii in range(len(lakex)):
                                    if lakex[ii] > xl and \
                                       lakex[ii] < xr and \
                                       lakey[ii] > yl and \
                                       lakey[ii] < yr:
                                        lakexl += [lakex[ii]]
                                        lakeyl += [lakey[ii]]
                                xla,yla = mapp(lakexl,lakeyl)
                                ax.fill(xla, yla, \
                                        color=(0.86,1.0,1.0,1.0))

                            # If no split
                            else:
                                xla,yla = mapp(lakex,lakey)
                                ax.fill(xla, yla, \
                                        color=(0.86,1.0,1.0,1.0))

                # If there are rivers
                if self.__river_exist:
                    for river in self.__river:
                        riverx = river[1]
                        rivery = river[0]
                        # If the river crosses the globe
                        if np.min(riverx) < -170. and \
                           np.max(riverx) > 170.:

                            lriversx = []
                            lriversy = []
                            llriversx = [riverx[0]]
                            llriversy = [rivery[0]]

                            if riverx[0] < 0.:
                                pos = -1.
                            else:
                                pos = 1

                            for ii in range(1,len(riverx)):

                                if pos < 0. and riverx[ii] > 0.:
                                    pos = 1.
                                    lriversx.append(llriversx)
                                    lriversy.append(llriversy)
                                    llriversx = [riverx[ii]]
                                    llriversy = [rivery[ii]]

                                elif pos > 0. and riverx[ii] < 0.:
                                    pos = -1.
                                    lriversx.append(llriversx)
                                    lriversy.append(llriversy)
                                    llriversx = [riverx[ii]]
                                    llriversy = [rivery[ii]]

                                else:
                                    llriversx.append(riverx[ii])
                                    llriversy.append(rivery[ii])
                                
                                if len(llriversx) > 0:
                                    lriversx.append(llriversx)
                                    lriversy.append(llriversy)

                            for rivx,rivy in zip(lriversx,lriversy):
                                xrv,yrv = mapp(rivx,rivy)
                                plt.plot(xrv,yrv, \
                                         color=(0.86,1.0,1.0,1.0))

                        # If the river is fine
                        else:
                            xrv,yrv = mapp(riverx,rivery)
                            ax.plot(xrv,yrv,color=(0.86,1.0,1.0,1.0))


            elif ii == 1:

                # Limits arrows
               #skipth = 20
                dTh =  (self.__pthrange[1] - self.__pthrange[0])/8.
                dth = self.__lat[1] - self.__lat[0]
                skipth = 0
                while skipth*dth < dTh:
                    skipth += 1
                    if skipth >= self.__nth-1:
                        break
               #skipch = 50
               #skipth = np.amax([self.__nth/skipth,1])
                dCh =  (self.__pchrange[1] - self.__pchrange[0])/4.
                dch = self.__lon[1] - self.__lon[0]
                skipch = 0
                while skipch*dch < dCh:
                    skipch += 1
                    if skipch >= self.__nch-1:
                        break

                # If Full map, cycle
                if self.__fullmapx:
                    plonv, lon = addcyclic(self.__wind[:,:,1], \
                                           self.__lon)
                    platv, lon = addcyclic(self.__wind[:,:,0], \
                                           self.__lon)
                    pmodv, lon = addcyclic(self.__wind[:,:,2], \
                                           self.__lon)
                    skipch = np.amax([(self.__nch+1)/skipch,1])
                else:
                    plonv = self.__wind[:,:,1]
                    platv = self.__wind[:,:,0]
                    pmodv = self.__wind[:,:,2]
                    skipch = np.amax([self.__nch/skipch,1])

                Lon, Lat = np.meshgrid(lon,lat)
                Lonv, Latv = plonv, platv
                xv, yv = mapp.rotate_vector(Lonv, Latv, Lon, Lat)

                xp = x[::skipth,::skipch]
                yp = y[::skipth,::skipch]
                xv = xv[::skipth,::skipch]
                yv = yv[::skipth,::skipch]
                scale = pmodv[::skipth,::skipch]

                # Adjust color
                cmap = self.__color_adjust(cm.rainbow, \
                                           np.min(scale), \
                                           np.max(scale), \
                                           self.__minwind, \
                                           self.__maxwind)

                # Contour data over the map.
                mapp.contour(x,y,pheight,[0.])

                cs = mapp.quiver(xp, yp, xv, yv, scale, cmap=cmap)

                # Color bar
                divider = make_axes_locatable(ax)
                cax = divider.append_axes("right", size="4%", \
                                          pad=0.05)
                cbar = plt.colorbar(cs,orientation='vertical',cax=cax)
                cax.yaxis.set_ticks_position('right')
                cbar.set_label('km/h')
                cax.yaxis.set_label_position('right')

            elif ii == 2:

                # If Full map, cycle
                if self.__fullmapx:
                    ptemp, lon = addcyclic(self.__temperature, \
                                           self.__lon)
                else:
                    ptemp = self.__temperature

                # Contour data over the map.
                mapp.contour(x,y,pheight,[0.])

                # Load color for terrain
                T0, T1 = self.__mintemperature, self.__maxtemperature
                cct =  self.__color_def_T(minv=T0, maxv=T1)

                # Contour data over the map.
                if T0 < 0:
                    cs = mapp.contourf(x,y,ptemp,255,cmap=cct)
                else:
                    cs = mapp.contourf(x,y,ptemp,255,cmap='Reds_r')

                # Color bar
                ticks = []
                minc = np.min(ptemp)
                maxc = np.max(ptemp)
                if minc < 0:
                    dec = int(minc)/10
                    mini = (dec+1)*10
                else:
                    mini = int(minc)-2
                maxi = int(maxc)+2
                for ff in range(mini,maxi,10):
                    if ff > minc:
                        ticks.append(ff)
                    elif ff > maxc:
                        break
                divider = make_axes_locatable(ax)
                cax = divider.append_axes("left", size="4%", pad=0.05)
                cbar = plt.colorbar(cs,orientation='vertical', \
                                    cax=cax,ticks=ticks)
                cax.yaxis.set_ticks_position('left')
                cbar.set_label('T')
                cax.yaxis.set_label_position('left')


            elif ii == 3:

                # If Full map, cycle
                if self.__fullmapx:
                    pmoist, lon = addcyclic(self.__moist, \
                                            self.__lon)
                else:
                    pmoist = self.__moist

                # Contour data over the map.
                mapp.contour(x,y,pheight,[0.])

                # Load color for terrain
                M0, M1 = np.amin(pmoist), np.amax(pmoist)

                # Adjust color
                cmap = self.__color_adjust(cm.plasma, M0, M1, \
                                           self.__minmoist, \
                                           self.__maxmoist)

                # Contour data over the map.
                cs = mapp.contourf(x,y,pmoist,255, \
                                   cmap=cmap)

                # Color bar
                ticks = []
                minc = np.min(pmoist)
                maxc = np.max(pmoist)
                mini = 0
                maxi = 100
                for ff in range(mini,maxi,10):
                    if ff > minc:
                        ticks.append(ff)
                    elif ff > maxc:
                        break
                divider = make_axes_locatable(ax)
                cax = divider.append_axes("right", size="4%", \
                                          pad=0.05)
                cbar = plt.colorbar(cs,orientation='vertical', \
                                    cax=cax,ticks=ticks)
                cax.yaxis.set_ticks_position('right')
                cbar.set_label('Moisture')
                cax.yaxis.set_label_position('right')

            elif ii == 4:

                # If Full map, cycle
                if self.__fullmapx:
                    pbiome, lon = addcyclic(self.__biome, \
                                            self.__lon)
                else:
                    pbiome = self.__biome

                # Contour data over the map.
                mapp.contour(x,y,pheight,[0.])

                # Load color and texturefor terrain
                cct =  self.__color_def_B()
                hatches = [None, \
                           None, \
                           None, \
                           None, \
                           None, \
                           '//' \
                           '//' \
                           '//', \
                           '//', \
                           '//', \
                           '//', \
                           None, \
                           None, \
                           None, \
                           None]

                # Contour data over the map.
                levels = [-3.5]
                for ii in range(12):
                    levels.append(levels[-1] + 1.)
                cs = mapp.contourf(x,y,pbiome, \
                                   levels=levels,vmin=-3.5,vmax=8.5, \
                                   cmap=cct,hatches=hatches)

                # Color bar
                ticks = [-3.0]
                for ii in range(12):
                    ticks += [ticks[-1] + 1.]
                divider = make_axes_locatable(ax)
                cax = divider.append_axes("left", size="4%", pad=0.05)
                cbar = plt.colorbar(cs,orientation='vertical', \
                                    cax=cax,ticks=ticks)
                cax.yaxis.set_ticks_position('left')
                cax.yaxis.set_label_position('left')
                labels = [' Ocean', \
                          'Co. Mt.', \
                          '    Mt.', \
                          ' Tundra', \
                          'Bor. F.', \
                          'Tp. Rf.', \
                          'Tr. Rf.', \
                          'Tp. Sf.', \
                          'Savanna', \
                          'Woodla.', \
                          'Grassl.', \
                          'Dessert' \
                         ]
                cbar.ax.set_yticklabels(labels)


        # Plot
        if fileout is not None:
            try:
                plt.savefig(fileout, format='png', dpi=dpi)
            except:
                plt.show()        
        else:
            plt.show()        
        plt.close('all')

