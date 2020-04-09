import config
import numpy as np
import json
from Drug import Drug
from Protein import Protein
import fileLoader as FL


def getValue(drugs, drugName, proteinName):
    if drugName in drugs:
        drug = drugs[drugName]
        if proteinName in drug.tar:
            value = drug.tar[proteinName].value - 4 if drug.tar[
                                                       proteinName].value >= 4 else 0  # all values underneath 4 are raised to 4

            return value if not drug.tar[
                proteinName].experimentalValue else 0  # experimental values will be handled with different weight

    return 0    #lowest value is 4


def buildMatrixPAct(drugs, allDrugs, allTargets):
    matrixPAct = np.zeros((len(allDrugs), len(allTargets)))
    for y in range(0, np.size(matrixPAct, 0)):  # iterate over column
        for x in range(0, np.size(matrixPAct, 1)):  # iterate over row
            matrixPAct[y, x] = getValue(drugs, allDrugs[y], allTargets[x])
    return matrixPAct


drugs = FL.getAllDrugs(config.DRUG_FILE_PATH)
allDrugs, allTargets = FL.loadFromFiles(config.DRUGS_LIST, config.TARGETS_LIST)
matrixPAct = buildMatrixPAct(drugs, allDrugs, allTargets)
np.save(config.MATRIX_PACT, matrixPAct)

#np.savetxt("./csv-files/pActMatrix.csv", matrixPAct, delimiter=";", fmt='%.3f')
print("Finished processing PAct Matrix")
