//This script adds a molecule table to the Archive

#@ MoleculeArchive archive
#@ MarsTable table

import de.mpg.biochem.mars.molecule.*
import de.mpg.biochem.mars.util.*

//The lock and unlock sequence will ensure the window is refreshed to reflect the changes.
archive.lock()

//Get the ImageMetadataUID for the ImageMetadata item at index 0
//This basically assumes you are working with an archive that only has
//one ImageMetadata item
metaUID = archive.getImageMetadata(0).getUID()

//Create a new molecule record using a new randomly generated UID key and the input table
//You can either create a molecule directly or use the archive.createMolecule method
//Since archives can contain different types of molecules using .createMolecule will always ensure
//the molecule you create is the correct type.
//molecule = new SingleMolecule(MarsMath.getUUID58(), table)
molecule = archive.createMolecule(MarsMath.getUUID58(), table)
molecule.setImageMetadataUID(metaUID);

//Add the molecule record to the archive given
archive.put(molecule)

archive.unlock()
