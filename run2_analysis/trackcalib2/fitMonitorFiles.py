from Utils import standard_text, center_message
from Utilities import OutputColoring as OC
from Paths import Paths  


# If the chi2 or relative error for any fitted parameter exceeds this value, warnings are printed and saved into warningFile
chi2_limit = 5.0 
relative_error_limit = 0.5
#If the error and value is smaller than the limit, warnings are printed and saved into warningFile
abs_error_limit = 0.00001
par_value_limit = 0.00001


class outputFile:

    @staticmethod
    def writeIntro(file,opts):
        file.write(standard_text["bold_line"] + '\n')
        file.write(center_message('Welcome to TrackEff Fitter') + '\n')
        file.write(standard_text["empty_line"] + '\n')
        file.write(center_message('TrackEff status') + '\n')
        if opts.sim_fit:
            file.write(center_message('Simultaneous fit active') + '\n')
        if opts.simple_sig:
            file.write(center_message('Simple signal fit active: Gauss pdf') + '\n')
        if opts.simple_bkg:
            file.write(center_message('Simple background fit active: linear pdf') + '\n')
        file.write(standard_text["empty_line"] + '\n')
        file.write(standard_text["bold_line"] + '\n\n\n')
        return
    
    def writeBinInfo(file,var):
        file.write(standard_text["line"] + '\n')
        if (var == ""):
            file.write(center_message('Integrated dataset fits') + '\n')
        else:
            file.write(center_message('Fit in bins of %s' % var.ljust(8)) + '\n')
        file.write(standard_text["line"] + '\n')
       
    def close(file):
        if file:
            file.close()
            file = None
        return

class monitorFiles:
    def __init__(self,opts,mode,method):
        self.status = open(Paths.getStatusFile(opts,mode, method), 'w')
        self.warnings = open(Paths.getWarningsFile(opts,mode, method), 'w')
    def writeIntro(self, opts):
        outputFile.writeIntro(self.status,opts)
        outputFile.writeIntro(self.warnings,opts)
    def writeBinInfo(self, var):
        outputFile.writeBinInfo(self.status,var)
        outputFile.writeBinInfo(self.warnings,var)
    def writeFitStatus(self,tag: str, string_bin: str,status: bool):
        self.status.write(f"Fit status {tag} \t{string_bin} status : {status}\n")
        if status != True:
            self.warnings.write(f"Fit status {tag} \t{string_bin} status : {status}\n")
            OC.get_warning_text(f"FAILED fit status {tag} \t{string_bin} status : {status}\n")
    def writeCHI2(self, tag, chi2):
        self.status.write(f"{tag} Chi2/NDOF \t{chi2}\n")
        if chi2 > chi2_limit:
            OC.get_warning_text(f"{tag} Chi2/NDOF too big!\tChi2 ={chi2}\n")
            self.warnings.write(f"Chi2 {tag} \t{chi2}\n")
    def writeEfficiency(self, result: list):
        """_summary_
        Writes the efficiency value with its erros into the status output file
        Args:
            result (list): [efficiency, high error, low error]
        """
        self.status.write(f"Efficiency \t{str(round(result[0], 4))}\t\t+{str(round(result[1], 4))}\t\t-{str(round(result[2], 4))}\n\n")

    def writeFitResult(self,minuit,
                     binx = -1,
                     biny = -1,
                     string_name = ""):
        string_bin = ""
        if binx != -1:
            if biny != -1: string_bin = f"bin {binx}\t bin \t {biny}\t"
            else:  string_bin = "bin {binx}\t"

        # write out the values
        for par in minuit.params:
            if par.merror is None:
                self.status.write(par.name.ljust(20) +
                                str(round(par.value, 3)).ljust(12) +
                                str(round(par.error, 3)).ljust(12))
            else:
                self.status.write(par.name.ljust(20) +
                                str(round(par.value, 3)).ljust(12) +
                                str(round(par.merror[0], 3)).ljust(12) +
                                str(round(par.merror[1], 3)).ljust(12))

            # check if a value or its error is too small
            if abs(par.value)< par_value_limit:
                self.warnings.write(par.name.ljust(20) + string_bin + " is zero!\n")
            if abs(par.error) < abs_error_limit:
                self.warnings.write(par.name.ljust(20) + string_bin + " error is zero!\n")

            # evaluate relative error and throw a warning in case of it is too large
            relative_error = 1
            if par.value != 0:
                relative_error = abs(par.error/par.value) 
            self.status.write("    RelErr".ljust(20) + str(round(relative_error, 3)).ljust(12) + "\n")
            if relative_error > relative_error_limit:
                self.warnings.write(string_name.ljust(8) + string_bin)
                self.warnings.write("RelErr " + par.name.ljust(20) + str(round(relative_error, 3)).ljust(12) + "\n")        
        return

    def close(self):
        outputFile.close(self.status)
        outputFile.close(self.warnings)
        return
