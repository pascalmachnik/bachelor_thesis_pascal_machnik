from Utilities import OutputColoring as OC
import Dictionaries
from pickle import load as pkl_load
import re

### output text functions
line_width = 64

standard_text = {}
standard_text["line"] = '#' + '-'*(line_width-2)+'#'
standard_text["empty_line"] = '#' + ' '*(line_width-2)+'#'
standard_text["bold_line"] = '#'*line_width


def center_message(text: str, linewidth: int = line_width)->str:
    return "#"+text.center(linewidth-2)+"#"


def left_message(text: str, linewidth: int = line_width)->str:
    return "#"+text.ljust(linewidth-2)+"#"


def title(text: str):
    OC.get_info_text(standard_text['bold_line'])
    OC.get_info_text(center_message(text))
    OC.get_info_text(standard_text['bold_line'])


def multiple_lines_title(*args):
    OC.get_info_text(standard_text['bold_line'])
    for text in args:
        OC.get_info_text(center_message(text))
    OC.get_info_text(standard_text['bold_line'])


def intro():
    print('''\n
    ████████╗██████╗  █████╗  ██████╗██╗  ██╗ ██████╗ █████╗ ██╗     ██╗██████╗ ██████╗ 
    ╚══██╔══╝██╔══██╗██╔══██╗██╔════╝██║ ██╔╝██╔════╝██╔══██╗██║     ██║██╔══██╗╚════██╗
       ██║   ██████╔╝███████║██║     █████╔╝ ██║     ███████║██║     ██║██████╔╝ █████╔╝
       ██║   ██╔══██╗██╔══██║██║     ██╔═██╗ ██║     ██╔══██║██║     ██║██╔══██╗██╔═══╝ 
       ██║   ██║  ██║██║  ██║╚██████╗██║  ██╗╚██████╗██║  ██║███████╗██║██████╔╝███████╗
       ╚═╝   ╚═╝  ╚═╝╚═╝  ╚═╝ ╚═════╝╚═╝  ╚═╝ ╚═════╝╚═╝  ╚═╝╚══════╝╚═╝╚═════╝ ╚══════╝
    ''')
    
    print("version: %s" % Dictionaries.version)
    
    separator = ", "
    print("authors: ", separator.join(Dictionaries.authors))

    print("contacts: ", end = '')
    print(*[str(name) + ': ' + str(email) for name,email in Dictionaries.contacts.items()], sep = separator)
    print("")

def argument_dict(parser):
    '''    
    Returns a dictionary with all the possible arguments and their corresponding default values and help. The sturcture is a nested dictionary:
    dictionary = { 'verbose': {
                        'default' = False,
                        'help' = "Increase output verbosity"
                    }
                }
    '''
    dict_tmp = {}

    for arg in parser._actions:
        dict_arg = {}
        dict_arg['name']    = arg.option_strings
        dict_arg['default'] = arg.default
        dict_arg['help']   = arg.help
        #dict_arg['type'] = arg.type
        #Boolean is not a supported type by argparse, so the type is None. Actually all types are useless in argparse, see https://bugs.python.org/issue24754
        dict_tmp[arg.dest] = dict_arg
    return dict_tmp


def insert_single_cut(in_cut, name, new_cut):
    '''
        Takes the cut string (in dict shape or not) and inserts the new cut at the correct place.
        New cut has to include only type of cut; eg add a Global cut, but never try to add Global+Probe cut!
    '''
    name = name + ":"
    
    # if dictionary split into strings for easier manipulation
    # checks for }; but only splits on ;
    sep_cuts = re.split(r'(?<=})\s*;', in_cut)
    out_cuts = []
        
    for old_cut in sep_cuts:
        # Case for only one cut
        # case for another Global cut already present 
        if (name in old_cut):
            str_idx = old_cut.find(name) + len(name)
            old_cut = old_cut[:str_idx] + new_cut + "," + old_cut[str_idx:]
        # Case for no Global cut present but encased in brackets
        elif "}" in old_cut:
            old_cut = old_cut[:-1] + ";" + name + new_cut + "}"
        # Case for no Global cut present and no brackets
        else:
            old_cut = old_cut + ";" +name + new_cut
        
        out_cuts.append(old_cut)
        
    out_cut = ";".join(out_cuts)
    OC.get_debug_text(f"old_cut {out_cut}")
    
    return out_cut

def insert_Global_cut(old_cut,new_cut):
    '''
        Adds a new Global cut to the old cut STRING
    '''
    return insert_single_cut(old_cut,"Global",new_cut)

def insert_Probe_cut(old_cut,new_cut):
    '''
        Adds a new Probe cut to the old cut STRING
    '''
    return insert_single_cut(old_cut,"Probe",new_cut)

def insert_Tag_cut(old_cut,new_cut):
    '''
        Adds a new Tag cut to the old cut STRING
    '''
    return insert_single_cut(old_cut,"Tag",new_cut)   

def insert_Mother_cut(old_cut,new_cut):
    '''
        Adds a new Mother cut to the old cut STRING
    '''
    return insert_single_cut(old_cut,"Tag",new_cut)         

   
def readPickleFile(filename):
    with open(filename, "rb") as f:
        while True:
            try:
                yield pkl_load(f)
            except EOFError:
                break            
    return


class binEdges:

    @staticmethod
    def getXYBin(bin_i, list2D_of_edges):
        '''
            Takes 2D list of bin edges and bin number i,
            returns the x_bin and y_bin coordinates
        '''
        bin1len = len(list2D_of_edges[0])-1
        bin_x = bin_i % bin1len
        bin_y = int(bin_i/bin1len)
        return bin_x, bin_y


    @staticmethod
    def getBoundaries2D(bin_i, list2D_of_edges):
        '''
            Takes 2D list of bin edges and bin number i,
            returns the corresponding boundaries
        '''        
        bin_x, bin_y = binEdges.getXYBin(bin_i,list2D_of_edges)     
        return [ 
            [list2D_of_edges[0][bin_x],list2D_of_edges[0][bin_x+1]],            
            [list2D_of_edges[1][bin_y],list2D_of_edges[1][bin_y+1]]
        ]       
    @staticmethod
    def getBoundaries1D(bin_i, list1D_of_edges):
        '''
            Takes 1D list of bin edges and bin number i,
            returns the corresponding boundaries
        '''
        return [list1D_of_edges[0][bin_i],list1D_of_edges[0][bin_i+1]]

    @staticmethod
    def getBoundaries(bin_i, list_of_edges):        
        '''
            Takes either a 1D or 2D list of bin edges and bin number i,
            returns the corresponding boundaries
        '''
        if (len(list_of_edges)==1):
            return binEdges.getBoundaries1D(bin_i,list_of_edges)
        elif (len(list_of_edges)==2):
            return binEdges.getBoundaries2D(bin_i,list_of_edges)
        else:
            raise(OC.get_error_text("Your bin edges do not have dimension of one or two! Abort."))            