from .helpers import create_folders, read_json

print("> Initiating required folders: ")
create_folders()

print("> Loading parameters from `congif.json` accessible as `congif` variable. ")
congif = read_json()