import os, shutil

ROOT = "/home/jtdinsmo/Dropbox (MIT)/GCE UROP/"

f = open("fig-sources.txt", 'r')
for line in f.readlines():
    if line == '':
        continue
    args = line[:-1].split(' ')
    if len(args) != 2:
        continue
    name, path = args
    if not os.path.isfile(ROOT + path):
        print("Could not find " + path)
        continue
    if path[-4:] not in [".eps", ".pdf"]:
        print("WARNING: .eps or .pdf file formats are preferred")
    shutil.copyfile(ROOT+path, "./"+name)

f.close()
