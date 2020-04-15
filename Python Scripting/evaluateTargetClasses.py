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


def display_target_classes(target_classes):
    counter = 0
    other_values = 0
    no_class_found = 0
    limited_target_classes = {}
    if isinstance(target_classes, dict):
        for key, value in target_classes.items():
            if key != '' and 30 <= value:
                limited_target_classes[key] = value
                counter += 1
            else:
                if key == "":
                    no_class_found += value
                else:
                    other_values += value
    elif isinstance(target_classes, list):
        for key, value in target_classes:
            if key != '' and 30 <= value:
                limited_target_classes[key] = value
                counter += 1
            else:
                if key == "":
                    no_class_found += value
                else:
                    other_values += value

    print("For " + str(no_class_found) + " Proteins were no classes found")
    print(str(other_values) + " Proteins had less than 3 Targets")
    limited_target_classes['Unknown'] = no_class_found
    limited_target_classes['Other'] = other_values
    plt.pie([v for v in limited_target_classes.values()], labels=[k for k in limited_target_classes.keys()],
            autopct=None)
    plt.show()
    for key, value in limited_target_classes.items():
        print(str(key) + " " + str(value))


def add_target_to_dict(dict, name, value):
    if name in dict:
        dict[name] += value
    else:
        dict[name] = value


def is_enzyme(target_class):
    split_classes = target_class.split(" ")
    for word in split_classes:
        if len(word) > 4 and (word[-3:] == "ASE" or word[-4:] == "ASES"):
            return True
    return False


def terms_in_line(list_of_terms, line):
    for term in list_of_terms:
        if term in line:
            return True
    return False


def generalize_target_classes(target_classes):
    generalized_target_classes = {}
    for target_class, occurrence in target_classes:
        target_class_upper = target_class.upper()
        if terms_in_line(["KINASE"], target_class_upper):
            add_target_to_dict(generalized_target_classes, "Kinase", occurrence)
        elif terms_in_line(["G-PROTEIN"], target_class_upper):
            add_target_to_dict(generalized_target_classes, "GPCR", occurrence)
        elif terms_in_line(["TRANSPORTER", "CHANNEL", "CARRIER", "TRANSPORT", "FACILITATOR", "SYMPORTER"],
                           target_class_upper):
            add_target_to_dict(generalized_target_classes, "Transporter", occurrence)
        elif is_enzyme(target_class_upper) or terms_in_line(["ENZYME", "P450"], target_class_upper):
            add_target_to_dict(generalized_target_classes, "Enzyme", occurrence)
        elif terms_in_line(["NUCLEAR HORMONE RECEPTOR"], target_class_upper):
            add_target_to_dict(generalized_target_classes, "Nuclear Hormone Receptor", occurrence)
        else:
            generalized_target_classes[target_class] = occurrence

    return generalized_target_classes


all_drugs, all_targets = FL.load_from_files(config.DRUGS_LIST, config.TARGETS_LIST)
target_classes_sorted = get_all_target_classes_sorted(all_targets)
print("There are " + str(len(target_classes_sorted)) + " classes before generalizing")
target_classes_sorted_generalized = generalize_target_classes(target_classes_sorted)
print("There are " + str(len(target_classes_sorted_generalized)) + " classes after generalizing")
# writeClassesToFile(config.TARGET_CLASSES_GENERAL_GROUPS, target_classes_sorted)
display_target_classes(target_classes_sorted_generalized)
