import config
import numpy as np
import json
from Drug import Drug
from Protein import Protein
import fileLoader as FL


def get_value(drugs, drug_name, protein_name, pact_median):
    if drug_name in drugs:
        drug = drugs[drug_name]
        if protein_name in drug.tar:
            if drug.tar[protein_name].experimental_value:
                return pact_median - config.PACT_THRESHOLD if pact_median > config.PACT_THRESHOLD else 0
            return drug.tar[protein_name].value - config.PACT_THRESHOLD if drug.tar[
                                                                               protein_name].value >= config.PACT_THRESHOLD else 0  # threshold can be found in config file

    return 0


def build_matrix_pact(drugs, all_drugs, all_targets, pact_median):
    matrix_pact = np.zeros((len(all_drugs), len(all_targets)))
    for y in range(0, np.size(matrix_pact, 0)):  # iterate over column
        for x in range(0, np.size(matrix_pact, 1)):  # iterate over row
            matrix_pact[y, x] = get_value(drugs, all_drugs[y], all_targets[x], pact_median)
    return matrix_pact



pact_median = FL.load_single_value_file(config.PACT_MEDIAN)
drugs = FL.get_all_drugs(config.DRUG_FILE_PATH)
all_drugs, all_targets = FL.load_from_files(config.DRUGS_LIST, config.TARGETS_LIST)
matrix_pact = build_matrix_pact(drugs, all_drugs, all_targets, pact_median)
matrix_pact_transposed = matrix_pact.transpose()
np.save(config.MATRIX_PACT, matrix_pact)
np.save(config.MATRIX_PACT_TRANSPOSED, matrix_pact_transposed)
np.savetxt("./csv-files/pact_matrix.csv", matrix_pact, delimiter=";", fmt='%.3f')
np.savetxt("./csv-files/pact_matrix_transposed.csv", matrix_pact_transposed, delimiter=";", fmt='%.3f')
print("Finished processing PAct Matrix")
