input_folder = os.path.abspath(config["input_folder"])
output_folder = os.path.abspath(config["output_folder"])
log_folder = os.path.join(output_folder, "logs")

if not os.path.isdir(output_folder):
    os.mkdir(output_folder)

if not os.path.isdir(log_folder):
    os.mkdir(log_folder)