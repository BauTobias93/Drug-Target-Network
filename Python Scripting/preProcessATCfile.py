
import config

def preProcessFile(filePath, newFilePath):
    file = open(filePath, "r")
    newFile = open(newFilePath, "w")
    for line in file:
        newLine = ""
        lineSplit = line.split("$")
        lineSplit[0] = lineSplit[0].replace(" ", "_")
        lineSplit[0] = lineSplit[0].replace("-", "_")
        for entry in lineSplit:
            newLine += str(entry) + str("$")

        newLine = newLine[:-1]  #delete last character
        newFile.write(newLine)

    file.close()
    newFile.close()

preProcessFile(config.ATC_CODE_FILE_PATH, config.ATC_CODE_FILE_PATH_NEW)