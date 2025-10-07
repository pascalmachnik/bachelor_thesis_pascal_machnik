from numpy import array as npArray
from numpy import mean as npMean

class PlotHelpers:

    @staticmethod    
    def compareDict(d1,d2):
        '''
        Takes two dictionaries and compares their keys.
        It also returns the list of shared keys with different content and 
        the list of shared keys with the same content
        '''
        d1_keys = set(d1.keys())
        d2_keys = set(d2.keys())
        shared_keys = d1_keys.intersection(d2_keys)
        d1_extra_keys = d1_keys - d2_keys
        d2_extra_keys = d2_keys - d1_keys
        modified = {o : [d1[o], d2[o]] for o in shared_keys if d1[o] != d2[o]}
        same = set(o for o in shared_keys if d1[o] == d2[o])
        return d1_extra_keys, d2_extra_keys, modified, same


    @staticmethod    
    def GetBinCenterWithErr(boundaries):
        '''
            Gets the bin boundariens and returns the bin center and half the bin width.
            Returns an array.
        '''
        boundaries = npArray(boundaries).T  
        center = npMean(boundaries,axis=0)   
        return center.tolist(), (center-boundaries[0]).tolist()
   
    @staticmethod    
    def getBinEdges(center_list):
        '''
            Takes the boundaries from the file and creates the list of bin edges
        '''
        center_list = npArray(center_list)  
        if (len(center_list.shape)==2): #if 1D
            return [list(dict.fromkeys(center_list.flatten()))]
        else: #2D
            x_edges = list(dict.fromkeys(center_list[:,0,:].flatten()))
            y_edges = list(dict.fromkeys(center_list[:,1,:].flatten()))
            return [x_edges,y_edges]
        # center_list[:,0,:] returns the 2D array of all the x-axis boundaries. 
        # flatten puts them into a 1D list
        # dict.fromkeys takes the 1D list and makes them into a key list. THIS KEEPS THE ORDER OF THE ORIGINAL LIST. Which is important in case something went wrong somewhere
        # Finally make it into a list and voila
    
    @staticmethod 
    def protectNegatives(val):
        if (val[0] == -1): return [0.0, 1.0, -1.0]
        else: return val

    @staticmethod    
    def get2Deff(eff_list, xBins):
        '''
            Tranforms the list of efficiencies into a 2D array
        '''
        return npArray(eff_list).reshape(xBins,-1).T 

    @staticmethod    
    def get2DeffErr(eff_err_list, xBins):
        '''
            Tranforms the list of efficiency errors into a 2D array
        '''
        new_err_list = []
        new_err_list.append(PlotHelpers.get2Deff(eff_err_list[0],xBins))
        new_err_list.append(PlotHelpers.get2Deff(eff_err_list[1],xBins))
        return new_err_list