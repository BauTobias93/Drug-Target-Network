import config


def pre_process_file(filePath, newFilePath):
    file = open(filePath, "r")
    new_file = open(newFilePath, "w")
    for line in file:
        new_line = ""
        line_split = line.split("$")
        line_split[0] = line_split[0].replace(" ", "_")
        line_split[0] = line_split[0].replace("-", "_")
        for entry in line_split:
            new_line += str(entry) + str("$")

        new_line = new_line[:-1]  # delete last character
        new_file.write(new_line)

    file.close()
    new_file.close()


pre_process_file(config.ATC_CODE_FILE_PATH_OLD, config.ATC_CODE_FILE_PATH)
