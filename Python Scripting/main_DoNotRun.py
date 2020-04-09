
import numpy as np
import config
import fileLoader as FL
import parseFile
import processPActMatrix
import processGOMatrix
import processTargetMatrices
import processDrugMatrices
import preProcessATCfile

#DO NOT RUN, NOT FINISHED YET

#TODO check if folder where files are going to be written is actually there

'''
#preprocess ATC file
preProcessATCfile.preProcessFile(config.ATC_CODE_FILE_PATH_OLD, config.ATC_CODE_FILE_PATH)

#parse drugs
drugs, allDrugs, allTargets = parseFile.parseDrugs(config.DRUGBANK_FILE_PATH, config.ATC_CODE_FILE_PATH, config.ATC_CODE_ADDITION_FILE_PATH)
FL.saveDrugsToFiles(config.DRUG_FILE_PATH, drugs)
FL.saveToFile(config.DRUGS_LIST, allDrugs)
FL.saveToFile(config.TARGETS_LIST, allTargets)

#build PAct value matrix
matrixPAct = processPActMatrix.buildMatrixPAct(drugs, allDrugs, allTargets)
np.save(config.MATRIX_PACT, matrixPAct)

np.savetxt("pActMatrix.csv", matrixPAct, delimiter=";")

#build Drug Distance matrix
matrixDrugDistance = processDrugMatrices.buildMatrixDrugDistance(drugs, matrixPAct, allDrugs, False)
matrixDrugATCDistance = processDrugMatrices.buildMatrixDrugDistance(drugs, matrixPAct, allDrugs, True)
np.save(config.MATRIX_DRUG_DISTANCE, matrixDrugDistance)
np.save(config.MATRIX_DRUG_ATC_DISTANCE, matrixDrugATCDistance)

np.savetxt("drugDistance.csv", matrixDrugDistance, delimiter=";")
np.savetxt("drugDistance_atc.csv", matrixDrugATCDistance, delimiter=";")


#build GO term matrix
dictGO = FL.loadGODictionary(config.GO_TERM_FILE_PATH)
allGOTerms = processGOMatrix.getAllGOTermNames(allTargets, dictGO)
FL.saveToFile(config.TARGETS_LIST, allTargets)
matrixGO = processGOMatrix.buildMatrixGO(allTargets, allGOTerms, dictGO)
np.save(config.MATRIX_GO_TERMS, matrixGO)

#build target matrix
matrixTargetGODistance = processTargetMatrices.buildMatrixTargetDistance(drugs, matrixPAct, allTargets, True)
np.save(config.MATRIX_TARGET_DISTANCE, matrixTargetGODistance)
'''