
import numpy as np
import config
import fileLoader as FL

def getAllGOTermNames(allTargets, dictGO):
    allGOTerms = []
    for i in range(0, len(allTargets)):
        if allTargets[i] in dictGO:
            for value in dictGO.get(allTargets[i]):
                if value not in allGOTerms:
                    allGOTerms.append(value)

    return allGOTerms

def buildMatrixGO(allTargets, allGOTerms, dictGO):
    countTargetNotFound = 0
    listTargetNotFound = []
    matrixGO = np.zeros((len(allTargets),len(allGOTerms)))
    for i in range(0, len(allTargets)):
        if i % 500 == 0:
            print("Reached Nr: " + str(i) + "/" + str(len(allTargets)))
        if allTargets[i] in dictGO:
            for value in dictGO.get(allTargets[i]):
                x = allGOTerms.index(value)
                matrixGO[i, x] = 1
        else:
            countTargetNotFound += 1
            listTargetNotFound.append(allTargets[i])

    print("Targets not found in GO terms list: " + str(countTargetNotFound))
    print("Targets not found:\n" + str(listTargetNotFound))
    return matrixGO


allDrugs, allTargets = FL.loadFromFiles(config.DRUGS_LIST, config.TARGETS_LIST)
dictGO = FL.loadGODictionary(config.GO_TERM_FILE_PATH)
allGOTerms = getAllGOTermNames(allTargets, dictGO)
FL.saveToFile(config.GO_TERMS_LIST, allGOTerms)
matrixGO = buildMatrixGO(allTargets, allGOTerms, dictGO)
np.save(config.MATRIX_GO_TERMS, matrixGO)

