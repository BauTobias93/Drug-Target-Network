
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
#TODO delete all computing parts of files so they can be imported without running it
#TODO get rid of all the not used imports

'''
#preprocess ATC file
preProcessATCfile.preProcessFile(config.ATC_CODE_FILE_PATH_OLD, config.ATC_CODE_FILE_PATH)

#parse drugs
drugs, all_drugs, all_targets = parseFile.parseDrugs(config.DRUGBANK_FILE_PATH, config.ATC_CODE_FILE_PATH, config.ATC_CODE_ADDITION_FILE_PATH)
FL.saveDrugsToFiles(config.DRUG_FILE_PATH, drugs)
FL.saveToFile(config.DRUGS_LIST, all_drugs)
FL.saveToFile(config.TARGETS_LIST, all_targets)

#build PAct value matrix
matrix_pact = processPActMatrix.buildMatrixPAct(drugs, all_drugs, all_targets)
np.save(config.MATRIX_PACT, matrix_pact)

np.savetxt("pActMatrix.csv", matrix_pact, delimiter=";")

#build Drug Distance matrix
matrix_drug_distance = processDrugMatrices.buildMatrixDrugDistance(drugs, matrix_pact, all_drugs, False)
matrix_drug_atc_distance = processDrugMatrices.buildMatrixDrugDistance(drugs, matrix_pact, all_drugs, True)
np.save(config.MATRIX_DRUG_DISTANCE, matrix_drug_distance)
np.save(config.MATRIX_DRUG_ATC_DISTANCE, matrix_drug_atc_distance)

np.savetxt("drugDistance.csv", matrix_drug_distance, delimiter=";")
np.savetxt("drugDistance_atc.csv", matrix_drug_atc_distance, delimiter=";")


#build GO term matrix
dictGO = FL.loadGODictionary(config.GO_TERM_FILE_PATH)
allGOTerms = processGOMatrix.getAllGOTermNames(all_targets, dictGO)
FL.saveToFile(config.TARGETS_LIST, all_targets)
matrixGO = processGOMatrix.buildMatrixGO(all_targets, allGOTerms, dictGO)
np.save(config.MATRIX_GO_TERMS, matrixGO)

#build target matrix
matrixTargetGODistance = processTargetMatrices.buildMatrixTargetDistance(drugs, matrix_pact, all_targets, True)
np.save(config.MATRIX_TARGET_DISTANCE, matrixTargetGODistance)
'''