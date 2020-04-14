import config
import fileLoader as FL
import operator
import matplotlib.pyplot as plt


def print_all_atc(drugs):
    categories = {}
    counter_no_atc = 0
    for d in drugs:
        if not drugs[d].atc:
            counter_no_atc += 1
        else:
            for atc in drugs[d].atc:
                if atc[0:1] not in categories:
                    categories[atc[0:1]] = 1
                else:
                    categories[atc[0:1]] += 1

    # categories = sorted(categories.items(), key=operator.itemgetter(1),reverse=True)
    print(categories)
    plt.pie([v for v in categories.values()], labels=[k for k in categories.keys()], autopct=None)
    plt.show()


drugs = FL.get_all_drugs(config.DRUG_FILE_PATH)
all_drugs, all_targets = FL.load_from_files(config.DRUGS_LIST, config.TARGETS_LIST)
print_all_atc(drugs)
