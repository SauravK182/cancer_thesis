def parse_bbduk():
    '''Function will parse the BBDuk stderr file and return a .txt file containing only the 4 lines recognized by MultiQC
    
    Parameters:
        bbduk:      A text file (captured by stderr stream) to isolate the lines deemed important by MultiQC
    
    Return value:   A text file with only the lines recognized by MultiQC for the BBDuk package
    '''