from __future__ import print_function



# # # # # # # # # # # # # # # # # # # 
#                                   #
#   Main part of the plot maker!    #
#                                   #
# # # # # # # # # # # # # # # # # # #
def make_figure_slides(folder, mainName, listOfEntries, defaultTitle, figsPerSlide, figsPerRow, figType, outputFile):
    
    #Get width of the pic based on the number of plots per row
    picWidth = 0.95/figsPerRow    
    
    #Pass in the mainName of the figure with __LOOP__ where the loop should happen
    lineList = []
    for num,x in enumerate(listOfEntries,start =1):
        name = mainName.replace("__LOOP__",x,1)
        newLine = "\\\\" if (num%figsPerRow==0) else "" #Put // at the end of each row
        lineList.append("\t\\includegraphics[width=" + "%.2f" %picWidth  +"\\textwidth]{" + folder + name + figType + "}"+newLine)

    #Open the file
    f = open(outputFile,'w')
    for num,line in enumerate(lineList, start = 1):
        if ((num-1)%figsPerSlide==0): #Print the begining of a slide
            print("\\begin{frame}{"+defaultTitle+"}",file = f) #TODO: Replace the list of titles by something sensible
            print("",file = f)
            print("\\centering",file = f)
        print(line, file =f) #Print the includegraphics
        if (num%figsPerSlide==0):  #Print the end of a slide
            print("\\end{frame}",file = f)
            print("%-----------------------------------------",file = f)
            print("",file = f)           
    if (len(lineList)%figsPerSlide !=0):
        print("\\end{frame}",file = f)
    f.close()
    
    #Drop everything after "." in the file name and print
    print("\\include{"+outputFile.split(".",1)[0]+"}")