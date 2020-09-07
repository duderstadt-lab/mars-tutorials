//Get number of molecules with tags printed in the script window

#@ LogService logService
#@ MoleculeArchive archive
#@ String(label="Tags (comma separated list)") Tags

String[] tagList = Tags.split(",")

for (int i=0; i<tagList.length; i++) {
   tagList[i] = tagList[i].trim();
}

int number = archive.getMoleculeUIDs().stream().filter{UID ->
                            int tagCount = 0
                       for (int i=0; i<tagList.length; i++) {
                            for (String tag : archive.get(UID).getTags()) {
                               if (tagList[i].equals(tag)) {
                                  tagCount++;
                               }
                            }
                         }
                    if (tagCount == tagList.length)
                                             return true;
                                        else
                                             return false;

                    }.count()

println("There are " + number + " molecules with tags " + Tags + " in archive " + archive.getName())
