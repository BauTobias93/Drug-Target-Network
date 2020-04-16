import config
import numpy as np
import json
from Drug import Drug
from Protein import Protein
import fileLoader as FL


def get_atc_value(first_atc_codes, second_atc_codes):
    best_value = 1.0
    for first_atc_code in first_atc_codes:
        for second_atc_code in second_atc_codes:
            value = 1.0
            for i in range(0, 7):
                if i < len(first_atc_code) and i < len(second_atc_code):
                    if i != 1 and i != 5:
                        if i == 2 or i == 6:
                            nbr1 = int(first_atc_code[i - 1]) * 10 + int(first_atc_code[i])
                            nbr2 = int(second_atc_code[i - 1]) * 10 + int(second_atc_code[i])
                            if nbr1 == nbr2:
                                value -= 0.2
                        elif first_atc_code[i] == second_atc_code[i]:
                            value -= 0.2
                        else:
                            break
                else:
                    break
            best_value = value if value < best_value else best_value

    return best_value


def get_atc_distance(drugs, all_drugs, i, j):
    if all_drugs[i] in all_drugs and all_drugs[j] in all_drugs:
        drug1 = drugs[all_drugs[i]]
        drug2 = drugs[all_drugs[j]]
        if drug1 is not None and drug1.atc and drug2 is not None and drug2.atc:
            return get_atc_value(drug1.atc, drug2.atc)
    return None


def get_processed_vectors(drugs, all_drugs, all_targets, matrix, i, j, process_drugs):
    first_row = np.copy(matrix[i])
    second_row = np.copy(matrix[j])
    experimental_value_row = np.full(np.size(matrix, 1), False)
    for col in range(0, np.size(matrix, 1)):
        if process_drugs:  # processing drug distances
            if all_targets[col] in drugs[all_drugs[i]].tar:
                first_target = drugs[all_drugs[i]].tar[all_targets[col]]
            else:
                first_target = Protein()
            if all_targets[col] in drugs[all_drugs[j]].tar:
                second_target = drugs[all_drugs[j]].tar[all_targets[col]]
            else:
                second_target = Protein()
        else:  # processing target distances
            if all_targets[i] in drugs[all_drugs[col]].tar:
                first_target = drugs[all_drugs[col]].tar[all_targets[i]]
            else:
                first_target = Protein()
            if all_targets[j] in drugs[all_drugs[col]].tar:
                second_target = drugs[all_drugs[col]].tar[all_targets[j]]
            else:
                second_target = Protein()
        if first_target.experimental_value or second_target.experimental_value:  # one or both values are experimentalValues
            experimental_value_row[col] = True
        if first_target.less_greater_than_sign != 0 or second_target.less_greater_than_sign != 0:
            if first_target.less_greater_than_sign == -1 and second_target.less_greater_than_sign == -1:  # < and <
                if first_target.value > second_target.value:
                    first_row[col] = second_row[col]
                else:
                    second_row[col] = first_row[col]
            elif first_target.less_greater_than_sign == -1 and first_target.value > second_target.value \
                    and (
                    second_target.less_greater_than_sign == 0 or second_target.less_greater_than_sign == 1):  # < and fixed && < and >
                first_row[col] = second_row[col]
            elif second_target.less_greater_than_sign == -1 and first_target.value < second_target.value \
                    and (
                    first_target.less_greater_than_sign == 0 or first_target.less_greater_than_sign == 1):  # fixed and < && < and >
                second_row[col] = first_row[col]

    return first_row, second_row, experimental_value_row

def get_weight_row_drug(drugs, all_drugs, all_targets, matrix, i, j):
    weight_row = np.full(np.size(matrix, 1), 1.0)
    for col in range(0, np.size(matrix, 1)):
        if all_targets[col] in drugs[all_drugs[i]].tar or all_targets[col] in drugs[all_drugs[j]].tar:
            # find protein targets in drugs or use default constructor
            if all_targets[col] in drugs[all_drugs[i]].tar:
                first_target = drugs[all_drugs[i]].tar[all_targets[col]]
            else:
                first_target = Protein()
            if all_targets[col] in drugs[all_drugs[j]].tar:
                second_target = drugs[all_drugs[j]].tar[all_targets[col]]
            else:
                second_target = Protein()

            #check how much the range will differ with respect to the signs
            if (first_target.less_greater_than_sign != 0 and second_target.less_greater_than_sign != 0) or \
                    first_target.experimental_value or second_target.experimental_value:
                weight_row[col] = config.SMALL_WEIGHT
            elif first_target.less_greater_than_sign != 0 or second_target.less_greater_than_sign != 0:
                weight_row[col] = config.MEDIUM_WEIGHT


    return weight_row



def get_drug_distance(matrix, i, j):
    magnitudes = np.linalg.norm(matrix[i]) * np.linalg.norm(matrix[j])
    return 1 - np.dot(matrix[i], matrix[j]) / magnitudes if magnitudes > 0 else 1


def get_weighted_drug_distance(drugs, all_drugs, all_targets, matrix, i, j):
    sum_of_weights = 0
    sum_of_differences = 0
    weight_row = get_weight_row_drug(drugs, all_drugs, all_targets, matrix, i, j)
    for col in range(0, np.size(matrix, 1)):
        sum_of_weights += weight_row[col]
        sum_of_differences += abs(matrix[i][col] - matrix[j][col]) * weight_row[col]

    return sum_of_differences / (sum_of_weights * config.SCALE_DOWN_FACTOR)\
        if sum_of_weights > 0 and config.SCALE_DOWN_FACTOR > 0 else 1


def build_matrix_drug_distance(drugs, matrix_pact, all_drugs, all_targets, with_atc_code):
    matrix_drugs = np.zeros((len(all_drugs), len(all_drugs)))
    for i in range(0, np.size(matrix_drugs, 0)):  # iterate row
        for j in range(i, np.size(matrix_drugs, 0)):  # iterate row
            # drug_distance = get_distance(matrix_pact, j, i)
            drug_distance = get_weighted_drug_distance(drugs, all_drugs, all_targets, matrix_pact, j, i)
            matrix_drugs[i, j] = drug_distance
            matrix_drugs[j, i] = drug_distance
            if with_atc_code:
                atc_distance = get_atc_distance(drugs, all_drugs, j, i)
                if atc_distance is not None:  # check if both drugs have an atc code
                    matrix_drugs[i, j] = 0.5 * matrix_drugs[i, j] + 0.5 * atc_distance
                    matrix_drugs[j, i] = 0.5 * matrix_drugs[j, i] + 0.5 * atc_distance
    return matrix_drugs


drugs = FL.get_all_drugs(config.DRUG_FILE_PATH)
all_drugs, all_targets = FL.load_from_files(config.DRUGS_LIST, config.TARGETS_LIST)
matrix_pact = FL.load_matrix(config.MATRIX_PACT)

matrix_drug_distance = build_matrix_drug_distance(drugs, matrix_pact, all_drugs, all_targets, False)
# matrix_drug_atc_distance = build_matrix_drug_distance(drugs, matrix_pact, all_drugs, all_targets, True)

np.savetxt("./csv-files/drugDistance.csv", matrix_drug_distance, delimiter=";", fmt='%.3f')
# np.savetxt("./csv-files/drugDistance_atc.csv", matrix_drug_atc_distance, delimiter=";", fmt='%.3f')

np.save(config.MATRIX_DRUG_DISTANCE, matrix_drug_distance)
# np.save(config.MATRIX_DRUG_ATC_DISTANCE, matrix_drug_atc_distance)

print("Drug Matrices were calculated")
