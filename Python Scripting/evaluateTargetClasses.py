import config
import fileLoader as FL
import operator
import matplotlib.pyplot as plt


def get_all_target_classes_sorted(all_targets):
    not_in_file_counter = 0
    all_target_classes = FL.load_and_build_dictionary(config.TARGET_CLASSES_PATH, 4, "|")
    target_classes = {}
    for target in all_targets:
        if target in all_target_classes:
            split_classes = all_target_classes[target].split(".")
            domain_class = split_classes[0]
            if domain_class in target_classes:
                target_classes[domain_class] += 1
            else:
                if domain_class == "Protein Kinase Superfamily":
                    print(target)
                target_classes[domain_class] = 1
        else:
            not_in_file_counter += 1

    print("Targets not found in file: " + str(not_in_file_counter))
    return sorted(target_classes.items(), key=operator.itemgetter(1), reverse=True)


def write_classes_to_file(path, target_classes):
    file = open(path, "w")
    file.write("target class|generalized target class\n")
    for tClass, amount in target_classes:
        if amount >= 10 and tClass != '':
            file.write(tClass + "|" + "UNKNOWN\n")
    file.close()


def print_target_classes(target_classes):
    for entry in target_classes:
        print(entry)

    counter = 0
    other_values = 0
    limited_target_classes = {}
    for key, value in target_classes:
        if key != '' and 3 <= value:
            limited_target_classes[key] = value
            counter += 1
        else:
            other_values += value

    limited_target_classes['Other'] = other_values
    print(limited_target_classes)
    plt.pie([v for v in limited_target_classes.values()], labels=[k for k in limited_target_classes.keys()],
            autopct=None)
    plt.show()


all_drugs, all_targets = FL.load_from_files(config.DRUGS_LIST, config.TARGETS_LIST)
target_classes_sorted = get_all_target_classes_sorted(all_targets)
# writeClassesToFile(config.TARGET_CLASSES_GENERAL_GROUPS, target_classes_sorted)
print_target_classes(target_classes_sorted)
