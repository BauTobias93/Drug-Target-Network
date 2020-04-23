from Drug import Drug
from Protein import Protein
import config
import fileLoader as FL
import sys

GLOABALCOUNTER = 0


# global GLOABALCOUNTER
# change GLOABALCOUNTER

def parse_protein_info(matches):
    split_matches = matches.split(";")
    proteins = {}
    # preprocess target groups
    for i in range(0, len(split_matches)):
        if "," in split_matches[i]:
            temp_targets_group = []
            split = split_matches[i].split(":")
            for entry in split:
                if "," in entry:
                    part_entry = entry.split(",")
                    if len(temp_targets_group) != 0:
                        # making sure that formats match
                        if len(part_entry) == len(temp_targets_group):
                            for j in range(0, len(part_entry)):
                                temp_targets_group[j] += ":" + str(part_entry[j])
                        else:
                            full_string = ""
                            for j in range(0, len(part_entry)):
                                full_string += part_entry[j] + ","
                            full_string = full_string[:-1]
                            for j in range(0, len(temp_targets_group)):
                                temp_targets_group[j] += ":" + full_string
                    else:
                        for j in range(0, len(part_entry)):
                            temp_targets_group.append(str(part_entry[j]))
                else:
                    if len(temp_targets_group) != 0:
                        for j in range(0, len(temp_targets_group)):
                            temp_targets_group[j] += ":" + str(entry)
                    else:
                        temp_targets_group.append(str(entry))
            # add all splited up targets to the matches, the first one overwrites the old one
            for j in range(0, len(temp_targets_group)):
                if j == 0:
                    split_matches[i] = str(temp_targets_group[j])
                else:
                    split_matches.append(str(temp_targets_group[j]))

    # start putting the targets in
    for match in split_matches:
        has_value = False
        match = match.split(":")
        protein = None
        if match[0] != '':
            protein = Protein(match[0])
            if len(match) > 1:
                if match[1] != '':
                    protein.mutation = match[1]
                if match[2] != '':
                    protein.species = match[2]
                if match[3] != '':
                    protein.type = match[3]
                if match[4] != '':
                    if '<' in match[4]:
                        match[4] = match[4][1:]
                        protein.less_greater_than_sign = -1
                    elif '>' in match[4]:
                        match[4] = match[4][1:]
                        protein.less_greater_than_sign = 1

                    protein.value = float(match[4])
                    has_value = True
                if match[5] != '':
                    protein.db_source = match[5]

        if protein is not None and not has_value:
            protein.value = 0
            protein.experimental_value = True
            has_value = True

        if protein is not None and has_value:  # proteins with no activity value will not be saved
            protein_name = protein.name.upper() + "_" + protein.species.upper()
            proteins[protein_name] = protein
    return proteins


def find_atc_for_drug_in_file(path, drug):
    file = open(path, "r")
    for line in file:
        line_split = line.split("$")
        if len(line_split) == 5 and drug.name.lower() == str(line_split[0]).lower() and line_split[3] != "":
            drug.atc = line_split[3].split(",")
            file.close()
            return 0

    file.close()
    return 1


def find_atc_for_drug(atc_path, drug, drugs_without_atc, atc_additional_path=None):
    count = 0
    count += find_atc_for_drug_in_file(atc_path, drug)
    if atc_additional_path is not None:
        count += find_atc_for_drug_in_file(atc_additional_path, drug)
    if count == 2 or (count == 1 and atc_additional_path is None):
        drugs_without_atc.append(drug.name)
        return 1
    else:
        return 0


def append_drug_name(all_drugs, drug_name):
    if drug_name.upper() not in all_drugs:
        all_drugs.append(drug_name.upper())


def append_target_names(all_targets, target_dict):
    for target in target_dict:
        if target not in all_targets:
            all_targets.append(target)


def add_target_values_to_list(values_list, targets):
    for key, target in targets.items():
        if target.value != 0:
            values_list.append(target.value)


def parse_drugs_info(drugs_path, atc_path, atc_additional_path=None):
    drugs = {}
    all_drugs = []
    all_targets = []
    pact_value_median = 0
    try:
        drug_file = open(drugs_path, "r")
        # counter for statistics
        count = 0
        no_protein_information_count = 0
        error_format_count = 0
        parsing_error_count = 0
        with_target_count = 0
        without_atc_code_count = 0
        drugs_without_atc = []
        pact_value_list = []
        drug_file.readline()  # read the first line, it is unnecessary header
        for line in drug_file:
            drug = None
            line = line.split("|")
            if len(line) != 6:  # check for correct format
                error_format_count += 1
                break
            if line[0] != '' and line[1] != '':  # check for name and status
                drug = Drug(line[0].upper(), line[1].upper())
                try:
                    if line[2] != '':  # check for target
                        tar = parse_protein_info(line[2])
                        drug.tar = tar
                        add_target_values_to_list(pact_value_list, tar)
                    if line[3] != '':  # check for enzyme
                        enz = parse_protein_info(line[3])
                        drug.enz = enz
                    if line[4] != '':  # check for transporter
                        tra = parse_protein_info(line[4])
                        drug.tra = tra
                    if line[5] != '\n':  # check for carrier
                        car = parse_protein_info(line[5].rstrip())
                        drug.car = car
                except:
                    parsing_error_count += 1
                    print(line)
                    print("Unexpected error:", sys.exc_info()[0])

                if not drug.tar and not drug.enz and not drug.tra and not drug.car:  # check if dicts are empty
                    no_protein_information_count += 1
            else:
                error_format_count += 1

            if drug is not None and drug.has_target_and_activity_value():
                without_atc_code_count += find_atc_for_drug(atc_path, drug, drugs_without_atc, atc_additional_path)
                drugs[drug.name.upper()] = drug
                append_drug_name(all_drugs, drug.name)
                append_target_names(all_targets, drug.tar)
                with_target_count += 1

            count += 1
        pact_value_list.sort()
        pact_value_median = pact_value_list[int(len(pact_value_list) / 2)]
        print_statistics(count, error_format_count, no_protein_information_count, parsing_error_count,
                         with_target_count,
                         without_atc_code_count, drugs_without_atc)

        drug_file.close()
    except:
        print("There was a problem processing the file: " + str(drugs_path) + " and/or " + str(
            atc_path) + " and/or " + str(atc_additional_path))
        print("Unexpected error:", sys.exc_info()[0])

    return drugs, all_drugs, all_targets, pact_value_median


def print_statistics(count, error_format_count, no_protein_information_count, parsing_error_count, with_target_count,
                     without_atc_code, drugs_without_atc):
    print("Drugs processed: " + str(count))
    print("Drugs processed with target: " + str(with_target_count))
    print("Drugs with parsing errors: " + str(parsing_error_count))
    print("Drugs without further information: " + str(no_protein_information_count))
    print("Lines with format errors: " + str(error_format_count))
    print("Drugs without atc code: " + str(without_atc_code))
    print("Drugs: " + str(drugs_without_atc))


drugs, all_drugs, all_targets, pact_value_median = parse_drugs_info(config.DRUGBANK_FILE_PATH,
                                                                    config.ATC_CODE_FILE_PATH,
                                                                    config.ATC_CODE_ADDITION_FILE_PATH)

print("Calculated median of pact values: " + str(pact_value_median))
print("I found " + str(len(all_drugs)) + " drugs")
print("I found " + str(len(all_targets)) + " targets")
FL.save_to_file(config.PACT_MEDIAN, pact_value_median)
FL.save_drugs_to_files(config.DRUG_FILE_PATH, drugs)
FL.save_to_file(config.DRUGS_LIST, all_drugs)
FL.save_to_file(config.TARGETS_LIST, all_targets)

# print("My gloabal counter found: " + str(GLOABALCOUNTER) + " appearances")
