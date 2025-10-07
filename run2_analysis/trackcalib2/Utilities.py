from colorama import init

#TODO: change to False when releasing
DebugMode = True

class OutputColoring:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    DEBUGPINK = '\033[35m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'
    # Background Colors:
    GREYBG = '\033[100m'
    REDBG = '\033[101m'
    GREENBG = '\033[102m'
    YELLOWBG = '\033[103m'
    BLUEBG = '\033[104m'
    PINKBG = '\033[105m'
    CYANBG = '\033[106m'
    # Special characters
    DOT = '\u2022'
    PM = '\u00B1'
    TRKCLB = '\33[92mTrackCalib2\033[0m'



    def __init__(self):
        init()
        
    @staticmethod
    def get_info_text(text: str) -> None: #I would call it print_info, as you don't return anything, but okay...
        print(f"[INFO]: {text}")

    @staticmethod
    def get_info_text_blue(text: str) -> None:
        print(f"{OutputColoring.OKBLUE}[INFO]: {text}{OutputColoring.ENDC}")

    @staticmethod
    def get_warning_text(text: str) -> None:
        print(f"{OutputColoring.WARNING}[WARNING]: {text}{OutputColoring.ENDC}")

    @staticmethod
    def get_error_text(text: str) -> None:
        print(f"{OutputColoring.FAIL}[ERROR]: {text}{OutputColoring.ENDC}")

    @staticmethod
    def get_ok_text(text: str) -> None:
        print(f"{OutputColoring.OKGREEN}[OK]: {text}{OutputColoring.ENDC}")
    
    @staticmethod
    def get_debug_text(text: str) -> None:
        if (DebugMode):
            print(f"{OutputColoring.DEBUGPINK}[DEBUG]: {text}{OutputColoring.ENDC}")


    @staticmethod
    def bold_text(text: str) -> str:
        return f"{OutputColoring.BOLD}{text}{OutputColoring.ENDC}"


class TextManipulation:
    @staticmethod
    def concat_list_to_string(string_list: list) -> str:
        string = ""
        try:
            for i in range(len(string_list)):
                if i == len(string_list)-1:
                    string += string_list[i]
                    break
                string += string_list[i] + ", "
        except:
            OutputColoring.get_error_text("Something is wrong with the list. The string should only contain strings!")
        
        return string

    @staticmethod
    def concat_dict_to_string(string_dict: dict) -> str:
        string = ""
        try:
            for key, value in string_dict.items():
                string += key + ": " + value + "\n"
        except:
            OutputColoring().get_error_text("Something is wrong with the dict. The dict should only contain strings as key/value pairs!")
        
        return string
