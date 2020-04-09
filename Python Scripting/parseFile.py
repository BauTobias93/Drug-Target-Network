from Drug import Drug
from Protein import Protein
import config
import fileLoader as FL
import sys

GLOABALCOUNTER = 0


# global GLOABALCOUNTER
# change GLOABALCOUNTER

def getProtInfo(matches):
    splitMatches = matches.split(";")
    proteins = {}
    # preprocess target groups
    for i in range(0, len(splitMatches)):
        if "," in splitMatches[i]:
            tempTargetsGroup = []
            split = splitMatches[i].split(":")
            for entry in split:
                if "," in entry:
                    partEntry = entry.split(",")
                    if len(tempTargetsGroup) != 0:
                        for j in range(0, len(partEntry)):
                            tempTargetsGroup[j] += ":" + str(partEntry[j])
                    else:
                        for j in range(0, len(partEntry)):
                            tempTargetsGroup.append(str(partEntry[j]))
                else:
                    if len(tempTargetsGroup) != 0:
                        for j in range(0, len(tempTargetsGroup)):
                            tempTargetsGroup[j] += ":" + str(entry)
                    else:
                        tempTargetsGroup.append(str(entry))
            # add all splited up targets to the matches, the first one overwrites the old one
            for j in range(0, len(tempTargetsGroup)):
                if j == 0:
                    splitMatches[i] = str(tempTargetsGroup[j])
                else:
                    splitMatches.append(str(tempTargetsGroup[j]))

    # start putting the targets in
    for match in splitMatches:
        hasValue = False
        match = match.split(":")
        protein = None
        if match[0] != '':
            protein = Protein(match[0])
            if len(match) > 1:
                if match[1] != '':
                    protein.mutation = match[1]
                if match[2] != '':
                    protein.species = match[2]
                if match[3] != '':
                    protein.type = match[3]
                if match[4] != '':
                    if '<' in match[4]:
                        match[4] = match[4][1:]
                        protein.value = 4
                    else:
                        if '>' in match[4]:
                            match[4] = match[4][1:]
                        protein.value = float(match[4])
                    hasValue = True
                if match[5] != '':
                    protein.dbSource = match[5]

        if protein != None and not hasValue:
            protein.value = 0
            protein.experimentalValue = True
            hasValue = True

        if protein != None and hasValue:  # proteins with no activity value will not be saved
            proteinName = protein.name.upper() + "_" + protein.species.upper()
            proteins[proteinName] = protein
    return proteins


def findATCforDrugInFile(path, drug):
    file = open(path, "r")
    for line in file:
        lineSplit = line.split("$")
        if len(lineSplit) == 5 and drug.name.lower() == str(lineSplit[0]).lower() and lineSplit[3] != "":
            drug.atc = lineSplit[3].split(",")
            file.close()
            return 0

    file.close()
    return 1


def findATCforDrug(atcPath, atcAdditionalPath, drug, drugsWithoutATC):
    count = 0
    count += findATCforDrugInFile(atcPath, drug)
    count += findATCforDrugInFile(atcAdditionalPath, drug)
    if count == 2:
        drugsWithoutATC.append(drug.name)
        return 1
    else:
        return 0


def appendDrugName(allDrugs, drugName):
    if drugName.upper() not in allDrugs:
        allDrugs.append(drugName.upper())


def appendTargetNames(allTargets, targetDict):
    for target in targetDict:
        if target not in allTargets:
            allTargets.append(target)


def parseDrugs(drugsPath, atcPath, atcAdditionalPath):
    drugs = {}
    allDrugs = []
    allTargets = []
    try:
        drugFile = open(drugsPath, "r")
        # counter for statistics
        count = 0
        noProteinInformationCount = 0
        errorFormatCount = 0
        parsingErrorCount = 0
        withTargetCount = 0
        withoutATCcodeCount = 0
        drugsWithoutATC = []
        drugFile.readline()  # read the first line, it is unnecessary header
        for line in drugFile:
            drug = None
            line = line.split("|")
            if line[0] == 'Idelalisib':
                print("What the fuck, dude?!")
            if len(line) != 6:  # check for correct format
                errorFormatCount += 1
                break
            if line[0] != '' and line[1] != '':  # check for name and status
                drug = Drug(line[0].upper(), line[1].upper())
                try:
                    if line[2] != '':  # check for target
                        tar = getProtInfo(line[2])
                        drug.tar = tar
                    if line[3] != '':  # check for enzyme
                        enz = getProtInfo(line[3])
                        drug.enz = enz
                    if line[4] != '':  # check for transporter
                        tra = getProtInfo(line[4])
                        drug.tra = tra
                    if line[5] != '\n':  # check for carrier
                        car = getProtInfo(line[5].rstrip())
                        drug.car = car
                except:
                    parsingErrorCount += 1
                    print(line)
                    print("Unexpected error:", sys.exc_info()[0])

                if not drug.tar and not drug.enz and not drug.tra and not drug.car:  # check if dicts are empty
                    noProteinInformationCount += 1
            else:
                errorFormatCount += 1

            if drug != None and drug.hasTargetAndActivityValue():
                withoutATCcodeCount += findATCforDrug(atcPath, atcAdditionalPath, drug, drugsWithoutATC)
                drugs[drug.name.upper()] = drug
                appendDrugName(allDrugs, drug.name)
                appendTargetNames(allTargets, drug.tar)
                withTargetCount += 1

            count += 1

        printStatistics(count, errorFormatCount, noProteinInformationCount, parsingErrorCount, withTargetCount,
                        withoutATCcodeCount, drugsWithoutATC)

        drugFile.close()
    except:
        print("There was a problem processing the file: " + str(drugsPath) + " and/or " + str(
            atcPath) + " and/or " + str(atcAdditionalPath))
        print("Unexpected error:", sys.exc_info()[0])

    return drugs, allDrugs, allTargets


def printStatistics(count, errorFormatCount, noProteinInformationCount, parsingErrorCount, withTargetCount,
                    withoutATCcode, drugsWithoutATC):
    print("Drugs processed: " + str(count))
    print("Drugs processed with target: " + str(withTargetCount))
    print("Drugs with parsing errors: " + str(parsingErrorCount))
    print("Drugs without further information: " + str(noProteinInformationCount))
    print("Lines with format errors: " + str(errorFormatCount))
    print("Drugs without atc code: " + str(withoutATCcode))
    print("Drugs: " + str(drugsWithoutATC))


drugs, allDrugs, allTargets = parseDrugs(config.DRUGBANK_FILE_PATH, config.ATC_CODE_FILE_PATH,
                                         config.ATC_CODE_ADDITION_FILE_PATH)
print("I found " + str(len(allDrugs)) + " drugs")
print("I found " + str(len(allTargets)) + " targets")
FL.saveDrugsToFiles(config.DRUG_FILE_PATH, drugs)
FL.saveToFile(config.DRUGS_LIST, allDrugs)
FL.saveToFile(config.TARGETS_LIST, allTargets)

#print("My gloabal counter found: " + str(GLOABALCOUNTER) + " appearances")
