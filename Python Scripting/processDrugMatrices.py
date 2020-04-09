import config
import numpy as np
import json
from Drug import Drug
from Protein import Protein
import fileLoader as FL


def getATCValue(atcCodes1, atcCodes2):
    bestValue = 0
    for atcCode1 in atcCodes1:
        for atcCode2 in atcCodes2:
            value = 0
            for i in range(0, 7):
                if i < len(atcCode1) and i < len(atcCode2):
                    if i != 1 and i != 5:
                        if i == 2 or i == 6:
                            nbr1 = int(atcCode1[i - 1]) * 10 + int(atcCode1[i])
                            nbr2 = int(atcCode2[i - 1]) * 10 + int(atcCode2[i])
                            if nbr1 == nbr2:
                                value += 0.2
                        elif atcCode1[i] == atcCode2[i]:
                            value += 0.2
                        else:
                            break
                else:
                    break
            bestValue = value if value > bestValue else bestValue

    return bestValue


def getATCDistance(drugs, allDrugs, i, j):
    if allDrugs[i] in allDrugs and allDrugs[j] in allDrugs:
        drug1 = drugs[allDrugs[i]]
        drug2 = drugs[allDrugs[j]]
        if drug1 != None and drug1.atc and drug2 != None and drug2.atc:
            return getATCValue(drug1.atc, drug2.atc)
    return None


def getDrugDistance(matrix, i, j):
    magnitudes = np.linalg.norm(matrix[j]) * np.linalg.norm(matrix[i])
    return 1 - np.dot(matrix[j], matrix[i]) / (magnitudes) if magnitudes > 0 else 1


def getWeightedDrugDistance(matrix, i, j):
    experimentalWeight = 0.25
    sumOfWeights = 0
    sumOfDiffences = 0
    for col in range(0, np.size(matrix, 1)):
        if matrix[i, col] == 0 or matrix[j, col] == 0:  # one of them or both are unknown values
            firstValue = matrix[i, col] if matrix[i, col] >= 4 else 4
            secondValue = matrix[j, col] if matrix[j, col] >= 4 else 4
            diff = abs(firstValue - secondValue)
            diff *= experimentalWeight
            sumOfWeights += 0.25
        else:
            diff = abs(matrix[i, col] - matrix[j, col])
            sumOfWeights += 1

        sumOfDiffences += diff

    return sumOfDiffences / sumOfWeights


def buildMatrixDrugDistance(drugs, matrixPAct, allDrugs, withATCcode):
    matrixDrugs = np.zeros((len(allDrugs), len(allDrugs)))
    for i in range(0, np.size(matrixDrugs, 0)):  # iterate over row
        for j in range(i, np.size(matrixDrugs, 0)):  # iterate over row
            drugDistance = getDrugDistance(matrixPAct, j, i)
            matrixDrugs[i, j] = drugDistance
            matrixDrugs[j, i] = drugDistance
            if withATCcode:
                ATCdistance = getATCDistance(drugs, allDrugs, j, i)
                if ATCdistance != None:  # check if both drugs have an atc code
                    matrixDrugs[i, j] = 0.5 * matrixDrugs[i, j] + 0.5 * ATCdistance
                    matrixDrugs[j, i] = 0.5 * matrixDrugs[j, i] + 0.5 * ATCdistance
    return matrixDrugs


drugs = FL.getAllDrugs(config.DRUG_FILE_PATH)
allDrugs, allTargets = FL.loadFromFiles(config.DRUGS_LIST, config.TARGETS_LIST)
matrixPAct = FL.loadMatrix(config.MATRIX_PACT)

matrixDrugDistance = buildMatrixDrugDistance(drugs, matrixPAct, allDrugs, False)
matrixDrugATCDistance = buildMatrixDrugDistance(drugs, matrixPAct, allDrugs, True)

np.savetxt("./csv-files/drugDistance.csv", matrixDrugDistance, delimiter=";", fmt='%.3f')
np.savetxt("./csv-files/drugDistance_atc.csv", matrixDrugATCDistance, delimiter=";", fmt='%.3f')

np.save(config.MATRIX_DRUG_DISTANCE, matrixDrugDistance)
np.save(config.MATRIX_DRUG_ATC_DISTANCE, matrixDrugATCDistance)
