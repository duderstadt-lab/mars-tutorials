//The following script create a completely empty SingleMoleculeArchive

#@output MoleculeArchive(label="archive.yama") archive
import de.mpg.biochem.mars.molecule.*

archive = new SingleMoleculeArchive("archive.yama")

// A SingleMoleculeArchive with one SdmmImageMetadata Item and one SingleMolecule record based on an input table is created in the following script.
#@ MarsTable table
#@output MoleculeArchive(label="archive.yama") archive

import de.mpg.biochem.mars.molecule.*
import de.mpg.biochem.mars.table.*
import de.mpg.biochem.mars.util.*

//Create a new molecule archive
archive = new SingleMoleculeArchive("archive.yama")

//Create ImageMetaDataTable this table must have all slices that molecule records can have.
//created with 1 column and 0 rows.
MarsTable metaTable = new MarsTable(1,0)
metaTable.setColumnHeader(0, "slice")

for (int slice=1;slice<=100;slice++) {
 //Since we initialised the table with no rows
 //We have to increase the table size by one row
 //before we can add more values
 metaTable.appendRow();

//Row index starts at 0
 metaTable.setValue("slice", slice - 1, slice);
}

//Create empty ImageMetaData with random UID
metaUID = MarsMath.getUUID58().substring(0, 10)
SdmmImageMetadata meta = new SdmmImageMetadata(metaUID, metaTable)

archive.putImageMetadata(meta);

//Create a new molecule record using table input
SingleMolecule molecule = new SingleMolecule(MarsMath.getUUID58(), table)
molecule.setImageMetadataUID(metaUID);

archive.put(molecule);
