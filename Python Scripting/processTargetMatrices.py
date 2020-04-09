import numpy as np
import config
import fileLoader as FL


def getTargetDistance(matrix, i, j):
    magnitudes = np.linalg.norm(matrix[:, j]) * np.linalg.norm(matrix[:, i])
    return 1 - np.dot(matrix[:, j], matrix[:, i]) / (magnitudes) if magnitudes > 0 else 1


def getGODistance(matrix, i, j):
    magnitudes = np.linalg.norm(matrix[j]) * np.linalg.norm(matrix[i])
    return 1 - np.dot(matrix[j], matrix[i]) / (magnitudes) if magnitudes > 0 else 1


def buildMatrixTargetDistance(matrixPAct, allTargets, matrixGO):
    matrixTargets = np.zeros((len(allTargets), len(allTargets)))
    print("Matrix Target Size: Row:" + str(np.size(matrixTargets, 0)) + "Col:" + str(np.size(matrixTargets, 1)))
    for i in range(0, np.size(matrixTargets, 1)):  # iterate over column
        for j in range(i, np.size(matrixTargets, 1)):  # iterate over column
            targetDistance = getTargetDistance(matrixPAct, i, j) * 0.5
            goDistance = getGODistance(matrixGO, i, j) * 0.5
            matrixTargets[i, j] = targetDistance + goDistance
            matrixTargets[j, i] = targetDistance + goDistance
    return matrixTargets


drugs = FL.getAllDrugs(config.DRUG_FILE_PATH)
allDrugs, allTargets = FL.loadFromFiles(config.DRUGS_LIST, config.TARGETS_LIST)
matrixPAct = FL.loadMatrix(config.MATRIX_PACT)
print("Matrix PACT Size: Row:" + str(np.size(matrixPAct, 0)) + "Col:" + str(np.size(matrixPAct, 1)))
matrixGO = FL.loadMatrix(config.MATRIX_GO_TERMS)
matrixTargetGODistance = buildMatrixTargetDistance(matrixPAct, allTargets, matrixGO)
np.save(config.MATRIX_TARGET_DISTANCE, matrixTargetGODistance)
np.savetxt("./csv-files/foo.csv", matrixTargetGODistance, delimiter=";")
