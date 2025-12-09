from .helpers import create_folders, read_json

print("> Initiating required folders: ")
create_folders()

print("> Loading parameters from `config.json` accessible as `config` variable. ")
congif = read_json()