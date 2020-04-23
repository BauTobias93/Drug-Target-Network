import numpy as np
import config
import fileLoader as FL
from Protein import Protein


def build_go_distance_matrix(matrix_go, all_targets):
    matrix_go_terms = np.zeros((len(all_targets), len(all_targets)))
    list_go_distances = []
    for i in range(0, np.size(matrix_go_terms, 0)):  # iterate row
        for j in range(i, np.size(matrix_go_terms, 0)):  # iterate row
            go_distance = get_target_distance(matrix_go, i,
                                              j)  # target distance and go distance share the same calculation
            if go_distance is not None:  # check if both drugs have an atc code
                matrix_go_terms[i, j] = go_distance
                matrix_go_terms[j, i] = go_distance
                list_go_distances.append(go_distance)

    list_go_distances.sort()
    go_median = list_go_distances[int(len(list_go_distances) / 2)]
    return matrix_go_terms, go_median


'''
list_go_distances.sort()
    medianIdx = int(len(list_go_distances) / 2)
    print("Median position is: " + str(medianIdx) + "/" + str(len(list_go_distances)))
    for i in range(0,len(list_go_distances)):
        if list_go_distances[i] == 1:
            print("Index found first 1: " +str(i) + "/" + str(len(list_go_distances)))
            break

    go_median = list_go_distances[int(len(list_go_distances) / 2)]
    print(go_median)
'''


def get_target_distance(matrix, i, j):
    magnitudes = np.linalg.norm(matrix[i]) * np.linalg.norm(matrix[j])
    return 1 - np.dot(matrix[i], matrix[j]) / magnitudes if magnitudes > 0 else 1


def get_weight_row_target(drugs, all_drugs, all_targets, matrix, i, j):
    weight_row = np.full(np.size(matrix, 1), 1)
    for col in range(0, np.size(matrix, 1)):
        if all_targets[i] in drugs[all_drugs[col]].tar or all_targets[j] in drugs[all_drugs[col]].tar:
            # find protein targets in drugs or use default constructor
            if all_targets[i] in drugs[all_drugs[col]].tar:
                first_target = drugs[all_drugs[col]].tar[all_targets[i]]
            else:
                first_target = Protein()
            if all_targets[j] in drugs[all_drugs[col]].tar:
                second_target = drugs[all_drugs[col]].tar[all_targets[j]]
            else:
                second_target = Protein()

            # check how much the range will differ with respect to the signs
            if (first_target.less_greater_than_sign != 0 and second_target.less_greater_than_sign != 0) or \
                    first_target.experimental_value or second_target.experimental_value:
                weight_row[col] = config.SMALL_WEIGHT
            elif first_target.less_greater_than_sign != 0 or second_target.less_greater_than_sign != 0:
                weight_row[col] = config.MEDIUM_WEIGHT

    return weight_row


def get_weighted_target_distance(drugs, all_drugs, all_targets, matrix, i, j):
    sum_of_weights = 0
    sum_of_differences = 0
    weight_row = get_weight_row_target(drugs, all_drugs, all_targets, matrix, i, j)
    for col in range(0, np.size(matrix, 1)):
        if matrix[i][col] != 0 or matrix[j][col] != 0:
            sum_of_weights += weight_row[col]
            sum_of_differences += abs(matrix[i][col] - matrix[j][col]) * weight_row[col]

    return sum_of_differences / (sum_of_weights * config.SCALE_DOWN_FACTOR) \
        if sum_of_weights > 0 and config.SCALE_DOWN_FACTOR > 0 else 1


def build_matrix_target_distance(drugs, all_drugs, all_targets, matrix_pact_transposed, matrix_go_terms, go_median,
                                 with_go_term):
    matrix_targets = np.zeros((len(all_targets), len(all_targets)))
    print("Matrix Target Size: Row:" + str(np.size(matrix_targets, 0)) + " Col:" + str(np.size(matrix_targets, 1)))
    for i in range(0, np.size(matrix_targets, 0)):  # iterate row
        for j in range(i, np.size(matrix_targets, 0)):  # iterate row
            #target_distance = get_target_distance(matrix_pact_transposed, i, j) * 0.5
            target_distance = get_weighted_target_distance(drugs, all_drugs, all_targets, matrix_pact_transposed, i,
                                                           j)
            matrix_targets[i, j] = target_distance
            matrix_targets[j, i] = target_distance
            if with_go_term:
                go_distance = matrix_go_terms[i, j] if matrix_go_terms[i, j] != 0 else go_median
                matrix_targets[i, j] = matrix_targets[i, j] * 0.5 + go_distance * 0.5
                matrix_targets[j, i] = matrix_targets[j, i] * 0.5 + go_distance * 0.5

    return matrix_targets


drugs = FL.get_all_drugs(config.DRUG_FILE_PATH)
all_drugs, all_targets = FL.load_from_files(config.DRUGS_LIST, config.TARGETS_LIST)
matrix_pact_transposed = FL.load_matrix(config.MATRIX_PACT_TRANSPOSED)
print("Matrix PACT Size: Row:" + str(np.size(matrix_pact_transposed, 0)) + " Col:" + str(
    np.size(matrix_pact_transposed, 1)))
matrix_go = FL.load_matrix(config.MATRIX_GO_TERMS)

matrix_go_terms, go_median = build_go_distance_matrix(matrix_go, all_targets)
matrix_target_distance = build_matrix_target_distance(drugs, all_drugs, all_targets, matrix_pact_transposed,
                                                      matrix_go_terms, go_median, False)
matrix_target_go_distance = build_matrix_target_distance(drugs, all_drugs, all_targets, matrix_pact_transposed,
                                                         matrix_go_terms, go_median, True)

np.save(config.MATRIX_TARGET_DISTANCE, matrix_target_distance)
np.save(config.MATRIX_TARGET_GO_DISTANCE, matrix_target_go_distance)

np.savetxt("./csv-files/target_distance.csv", matrix_target_distance, delimiter=";", fmt='%.3f')
np.savetxt("./csv-files/target_go_distance.csv", matrix_target_go_distance, delimiter=";", fmt='%.3f')
