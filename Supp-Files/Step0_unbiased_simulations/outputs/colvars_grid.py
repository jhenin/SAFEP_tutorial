import numpy as np

class colvars_grid:
    '''
    Manipulate gridded data from the Colvars Module

    The data lists may have one or several elements depending on the number of data series.
    PMF files contain one data series.
    Gradient files contain as many data series as there are dimensions

    # 3d surface plot for a PMF
    fig = plt.figure()
    from mpl_toolkits.mplot3d import Axes3D
    ax = fig.add_subplot(111, projection='3d')
    ax.plot_surface(*pmf.meshgrid(), np.clip(pmf.data[0], None, 40))

    Args:
        filename (str): an optional filename to parse
    
    Attributes:
        filename (str): name of source file
        data (list of array): data for last frame
        histdata (list of array): data for all frames (array shape: nframes, nx[0], nx[1]...)
        dim (int): array dimension
        nx (list of int): data sizes
        nframes (int): n
    '''

    def __init__(self, filename=None):
        self.reset()
        if filename is not None:
            self.read(filename)


    def reset(self, ):
        '''Reset all data'''
        self.filename = None
        self.data = []
        self.dim = 0
        self.nx = []
        self.xmin = []
        self.dx = []
        self.pbc = []
        self.histdata = []
        self.nframes = 0


    def summary(self):
        '''Print a short summary of the contents of the object'''
        print('Source file:', self.filename)
        print('Grid dimension:', self.dim, ['PBC' if p else 'non-PBC' for p in self.pbc])
        print('Number of data series:', len(self.data))
        print('Number of grid points:', np.prod(self.nx), self.nx)
        print('Grid spacings:', self.dx)
        print('Number of time frames:', self.nframes)


    def axes(self):
        '''Returns the axes as a list of 1d meshes'''
        return [np.array([self.xmin[i] + (k+0.5)*self.dx[i] for k in range(self.nx[i])]) for i in range(self.dim)]


    def meshgrid(self, ):
        '''Returns a mesh grid suitable for plotting the data'''
        return np.meshgrid(*self.axes(), copy=True, indexing='ij')


    def read(self, filename):
        '''Read data from a Colvars multicolumn file'''
        self.reset()
        self.filename = filename
        with open(filename) as f:
            l = f.readline().split()
            assert len(l) == 2 
            self.dim = int(l[1])
            for _ in range(self.dim):
                l = f.readline().split()
                assert(len(l) == 5)
                self.xmin.append(float(l[1]))
                self.dx.append(float(l[2]))
                self.nx.append(int(l[3]))
                self.pbc.append(l[4] == "1")
            f.close()
        rawdata = np.loadtxt(filename).transpose()
        nf = len(rawdata[0]) / np.prod(self.nx)
        assert nf == int(nf), 'Data size {len(rawdata[0])} is not a multiple of grid size {np.prod(self.nx)}' # I deleted f before string
        self.nframes = int(nf)
        # keep only values of gridded quantities
        for d in rawdata[self.dim:]:
            self.histdata.append(np.copy(d.reshape([self.nframes] + self.nx)))
            # Current data is last frame of last inserted data series
            self.data.append(self.histdata[-1][-1])
