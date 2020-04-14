import numpy as np
import config
import fileLoader as FL
from Protein import Protein
from processDrugMatrices import get_distance


def get_go_distance(matrix, i, j):
    magnitudes = np.linalg.norm(matrix[j]) * np.linalg.norm(matrix[i])
    return 1 - np.dot(matrix[j], matrix[i]) / (magnitudes) if magnitudes > 0 else 1


def build_matrix_target_distance(drugs, all_drugs, all_targets, matrix_pact_transposed, matrix_go):
    matrix_targets = np.zeros((len(all_targets), len(all_targets)))
    print("Matrix Target Size: Row:" + str(np.size(matrix_targets, 0)) + "Col:" + str(np.size(matrix_targets, 1)))
    for i in range(0, np.size(matrix_targets, 1)):  # iterate over column
        for j in range(i, np.size(matrix_targets, 1)):  # iterate over column
            target_distance = get_distance(drugs, all_drugs, all_targets, matrix_pact_transposed, i, j, False) * 0.5
            go_distance = get_go_distance(matrix_go, i, j) * 0.5
            matrix_targets[i, j] = target_distance + go_distance
            matrix_targets[j, i] = target_distance + go_distance
    return matrix_targets


drugs = FL.get_all_drugs(config.DRUG_FILE_PATH)
all_drugs, all_targets = FL.load_from_files(config.DRUGS_LIST, config.TARGETS_LIST)
matrix_pact_transposed = FL.load_matrix(config.MATRIX_PACT_TRANSPOSED)
print("Matrix PACT Size: Row:" + str(np.size(matrix_pact_transposed, 0)) + "Col:" + str(
    np.size(matrix_pact_transposed, 1)))
matrix_go = FL.load_matrix(config.MATRIX_GO_TERMS)
matrix_target_go_distance = build_matrix_target_distance(drugs, all_drugs, all_targets, matrix_pact_transposed,
                                                         matrix_go)
np.save(config.MATRIX_TARGET_DISTANCE, matrix_target_go_distance)
np.savetxt("./csv-files/target_distance.csv", matrix_target_go_distance, delimiter=";", fmt='%.3f')
