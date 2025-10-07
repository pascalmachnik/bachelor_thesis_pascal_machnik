import sys
from ArgParser import parse_args
from Manual import printManual
from Utils import title, intro
from Prepare import Prepare
from Fit import Fit
from Plot import Plot
from Systematics import runSystematics


def main():
    options = parse_args(sys.argv[1:])
    if hasattr(options, 'task'):
        intro()
        #Printing the manual
        if options.task == "manual":
            printManual()
            return
        title(options.task.upper() + " action has been selected!")
        
        #Running the systematics
        if (options.systematics):            
            runSystematics(options)
            return

        if options.task == "prepare":
            for mode in options.mode:
                Prepare(mode=mode, opts = options)           
        elif options.task == "fit":
            for method in options.method:
                for mode in options.mode:
                    Fit(method=method, mode=mode, opts=options)
        elif options.task == "plot":      
            Plot(options)

    else:
        print("trackcalib manual")


if __name__ == '__main__':
    main()
