import os
os.chdir(os.getcwd())
import xarray as xr
import numpy as np

class PlaRegion:
    """
    Geographical region definitions for PlaSim-LSG.

    This class requires the files
        - wet.nc (land-sea mask)
        - "basins_scalar.nc" (ocean basin definitions on T point grid)
        - "basins_vector.nc" (ocean basin definitions on U point grid)

    Initializing the class provides predefined areas and masks, for example:

    mask = PlaRegion()
    atl = mask.Atlantic3D() # 3D Atlantic basin mask
    """
    def __init__(self, path="../../grid/"):
        self.path = path
        self.lsm3d = xr.load_dataset(path+"wet.nc").wet.isel(time=0)
        self.lsm2d = self.lsm3d.isel(depth=0)
        self.lat = self.lsm2d.lat
        self.lon = self.lsm2d.lon
        self.depth = self.lsm3d.depth
        self.basins_sca = xr.load_dataset(path+"basins_scalar.nc").om 
        self.basins_vec = xr.load_dataset(path+"basins_vector.nc").om
        self.area = {}
        self.area["Labrador"] = ((360-70, 360-45), (52, 70))
        self.area["Irminger"] = ((360-45, 360-13), (55, 65))
        self.area["NorwegianW"] = ((360-8, 360), (55, 75))
        self.area["NorwegianE"] = ((0, 23), (55, 75))
        self.area["GreenSeaW"] = ((335, 360), (68, 80))
        self.area["GreenSeaE"] = ((0, 10), (68, 80))
        self.area["Europe"] = ((0, 35), (35, 70))
        self.area["EuropeN"] = ((0, 35), (50, 70))
        self.area["EuropeS"] = ((0, 35), (35, 50))
        self.area["England_ipla"] = (0, 6) # indices with respect to pla.lon, pla.lat
        self.area["Iceland_ipla"] = (61, 4)
        self.area["Nuuk_ipla"] = (55, 4)
        self.area["Moscow_ipla"] = (7, 6)
        self.area["Kuroshio"] = ((135,165), (22,40))
        self.area["KuroshioS"] = ((120, 130), (10, 30))
        self.area["KuroshioN"] = ((130, 150), (30, 40))
        
        
    def boolmask(self, mask, thresh=0.5):
        """
        Turns a mask of numerical values into boolean values.
        All values larger than `thresh` become True, otherwise False.
        """
        return np.ma.masked_array(mask.to_numpy(), mask.to_numpy() > thresh).mask
    
    def basins(self, points):
        if points == 'vec':
            return self.basins_vec
        elif points =='sca':
            return self.basins_sca
    
    def Atlantic(self, type='sca', southern_border=-34, hudson_bay=False):
        mask = np.full((self.basins(type)).shape, False)
        mask[np.where(self.basins(type) == 1)] = True
        mask[np.where(self.basins(type) == 4)] = True
        mask[np.where((self.lat < southern_border))] = False
        self.exclude_area(mask, (0,40), (20,40)) # Mediterranean I
        self.exclude_area(mask, (355,360), (20,40)) # Mediterranean II
        if not hudson_bay:
            self.exclude_area(mask, (265,283), (50,64)) # Hudson Bay
        self.exclude_area(mask, (0,360), (80,90)) # north of 80N
        self.exclude_area(mask, (20,285), (65,90)) # Arctic Ocean
        return mask

    def AtlanticArctic(self, type='sca', southern_border=-34, hudson_bay=False):
        mask = np.full((self.basins(type)).shape, False)
        mask[np.where(self.basins(type) == 1)] = True
        mask[np.where(self.basins(type) == 4)] = True
        mask[np.where((self.lat < southern_border))] = False
        self.exclude_area(mask, (0,40), (20,40)) # Mediterranean I
        self.exclude_area(mask, (355,360), (20,40)) # Mediterranean II
        if not hudson_bay:
            self.exclude_area(mask, (265,283), (50,64)) # Hudson Bay
        return mask
    
    def Pacific(self, type='sca', southern_border=-90):
        mask = np.full((self.basins(type)).shape, False)
        mask[np.where(self.basins(type) == 2)] = True
        mask[np.where((self.lat < southern_border))] = False
        return mask
    
    def Indian(self, type='sca', southern_border=-90):
        mask = np.full((self.basins(type)).shape, False)
        mask[np.where(self.basins(type) == 3)] = True
        mask[np.where((self.lat < southern_border))] = False
        return mask
    
    def Arctic(self, type='sca'):
        mask = np.full((self.basins(type)).shape, False)
        mask[np.where(self.basins(type) == 4)] = True
        self.exclude_area(mask, (315,360), (40,80)) # SW of Svalbard
        self.exclude_area(mask, (0,20), (40,80)) # SW of Svalbard
        return mask
    
    def StommelPolar(self, type='sca', southern_border=40, arctic=False):
        mask = np.full((self.basins(type)).shape, False)
        if arctic:
            mask = self.AtlanticArctic(type=type, southern_border=southern_border)
        else:
            mask = self.Atlantic(type=type, southern_border=southern_border)
        return mask
    
    def StommelEquatorial(self, type='sca',
        northern_border=40, southern_border=-34, maxdepth=1000):
        mask = np.full((self.basins(type)).shape, False)
        mask = self.Atlantic(type=type, southern_border=southern_border)
        mask[np.where((self.lat>=northern_border))] = False
        mask[:,:,np.where(self.depth>maxdepth)] = False
        return mask
    
    def Southern(self, type='sca', northern_lat=-60):
        mask = np.full((self.basins(type)).shape, False)
        mask[np.where(self.basins(type) > 0)] = True
        mask[np.where((self.lat>=northern_lat))] = False
        return mask

    def Labrador(self, type='sca'):
        mask = np.full((self.basins(type)).shape, False)
        mask = self.Atlantic(type=type)
        mask[np.where((self.lat<=52)),:] = False
        mask[np.where((self.lat>=70)),:] = False
        mask[:,np.where((self.lon<=360-70))] = False
        mask[:,np.where((self.lon>=360-45))] = False
        return mask
    
    def Irminger(self, type='sca'):
        mask = np.full((self.basins(type)).shape, False)
        mask = self.Atlantic(type=type)
        mask[np.where((self.lat<=55)),:] = False
        mask[np.where((self.lat>=65)),:] = False
        mask[:,np.where((self.lon<=360-45))] = False
        mask[:,np.where((self.lon>=360-13))] = False
        return mask
    
    def Norwegian(self, type='sca'):
        mask = np.full((self.basins(type)).shape, False)
        mask = self.Atlantic(type=type)
        mask[np.where((self.lat<=55)),:] = False
        mask[np.where((self.lat>=75)),:] = False
        mask[:,np.where((self.lon<=360-8) & (self.lon>23))] = False
        return mask
    
    def GreenSea(self, type='sca'):
        mask = np.full((self.basins(type)).shape, False)
        mask = self.Atlantic(type=type)
        mask[np.where((self.lat<=68)),:] = False
        mask[np.where((self.lat>=80)),:] = False
        mask[:,np.where((self.lon<=360-25) & (self.lon>10))] = False
        return mask
    
    def GoodHope(self, type='sca'):
        mask = np.full((self.basins(type)).shape, False)
        mask = self.Oceans(type=type)
        mask[np.where((self.lat<=-40)),:] = False
        mask[np.where((self.lat>=-30)),:] = False
        mask[:,np.where((self.lon<=20))] = False
        mask[:,np.where((self.lon>=40))] = False
        return mask
    
    def Kuroshio(self, type='sca'):
        mask = np.full((self.basins(type)).shape, False)
        mask = self.Pacific(type=type)
        mask[np.where((self.lat<=22)),:] = False
        mask[np.where((self.lat>=40)),:] = False
        mask[:,np.where((self.lon<=135))] = False
        mask[:,np.where((self.lon>=165))] = False
        return mask

    def Oceans(self, type='sca'):
        mask = np.full((self.basins(type)).shape, False)
        mask[np.where(self.basins(type) > 0)] = True
        return mask
    
    def Oceans3D(self, type='sca'):
        return self.extend_vertically(self.Oceans(type=type))
    
    def Atlantic3D(self, type='sca', southern_border=-34):
        return self.extend_vertically(self.Atlantic(type=type, southern_border=southern_border))
    
    def Pacific3D(self, type='sca', southern_border=-90):
        return self.extend_vertically(self.Pacific(type=type, southern_border=southern_border))
    
    def Indian3D(self, type='sca', southern_border=-90):
        return self.extend_vertically(self.Indian(type=type, southern_border=southern_border))
    
    def extend_vertically(self, mask):
        mask3d = np.repeat(mask[np.newaxis,:,:], len(self.depth), axis=0)
        return np.multiply(self.lsm3d, mask3d)
    
    def exclude_area(self, mask, lon, lat):
        for i in range(len(self.lat)):
            if (self.lat[i] >= lat[0]) & (self.lat[i] < lat[1]):
                mask[i, np.where((self.lon>=lon[0]) & (self.lon<lon[1]))] = False
    
    def dataarray(self, array):
        return xr.DataArray(array, dims=(["lat", "lon"]),
            coords={"lat" : self.lat, "lon" : self.lon})
    
    def nanmask(self, mask):
        return np.where(mask, 1.0, np.nan)
    
    def to3d(self, mask):
        """alias for extend_vertically"""
        return self.extend_vertically(mask)
    
    def lsgslice(self, areaname):
        coords = self.area[areaname]
        return slice(coords[0][0], coords[0][1]), slice(coords[1][0], coords[1][1])

    def plaslice(self, areaname):
        coords = self.area[areaname]
        return slice(coords[0][0], coords[0][1]), slice(coords[1][1], coords[1][0])