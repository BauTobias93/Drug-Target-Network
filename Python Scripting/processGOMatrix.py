import numpy as np
import config
import fileLoader as FL


def get_all_go_term_names(all_targets, dict_go):
    all_go_terms = []
    for i in range(0, len(all_targets)):
        if all_targets[i] in dict_go:
            for value in dict_go.get(all_targets[i]):
                if value not in all_go_terms:
                    all_go_terms.append(value)

    return all_go_terms


def build_matrix_go(all_targets, all_go_terms, dict_go):
    count_target_not_found = 0
    list_target_not_found = []
    matrix_go = np.zeros((len(all_targets), len(all_go_terms)))
    for i in range(0, len(all_targets)):
        if i % 500 == 0:
            print("Reached Nr: " + str(i) + "/" + str(len(all_targets)))
        if all_targets[i] in dict_go:
            for value in dict_go.get(all_targets[i]):
                x = all_go_terms.index(value)
                matrix_go[i, x] = 1
        else:
            count_target_not_found += 1
            list_target_not_found.append(all_targets[i])

    print("Targets not found in GO terms list: " + str(count_target_not_found))
    print("Targets not found:\n" + str(list_target_not_found))
    return matrix_go


all_drugs, all_targets = FL.load_from_files(config.DRUGS_LIST, config.TARGETS_LIST)
dict_go = FL.load_go_dictionary(config.GO_TERM_FILE_PATH)
all_go_terms = get_all_go_term_names(all_targets, dict_go)
FL.save_to_file(config.GO_TERMS_LIST, all_go_terms)
matrix_go = build_matrix_go(all_targets, all_go_terms, dict_go)
np.save(config.MATRIX_GO_TERMS, matrix_go)
