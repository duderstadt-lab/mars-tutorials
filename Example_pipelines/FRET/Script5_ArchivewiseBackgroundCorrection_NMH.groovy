// Import the required dependencies
#@ MoleculeArchive archive

import de.mpg.biochem.mars.molecule.*
import de.mpg.biochem.mars.table.*
import de.mpg.biochem.mars.util.*
import org.scijava.table.*
import groovy.lang.*

// Create a list containing the names of all videos present in the Archive
// make this a unique list

def names_list = []

archive.metadata().forEach{metadata ->
	names_list.add(metadata.getImage(0).getName())
}
names_list.unique()

// Iterate through all videos
// For each video, create a list with the UIDs of the three corresponding archives [FRET, AO, DO]
// Calculate the background values from the FRET archive and correct all archives with this value

for (def name : names_list) {

	// Create the [FRET, AO, DO] UID list for the video
	String[] metadataUIDs = new String[3]
	for (def mUID1 : archive.getMetadataUIDs()) {
		def metadata = archive.getMetadata(mUID1)
		if (metadata.getImage(0).getName().equals(name)) {
				if (metadata.hasTag("FRET"))
					metadataUIDs[0] = mUID1
				else if (metadata.hasTag("AO"))
					metadataUIDs[1] = mUID1
				else if (metadata.hasTag("DO"))
					metadataUIDs[2] = mUID1
		}
	}

	// Calculate the correction factors from the FRET archive
	if (metadataUIDs[0] != null) {
		def metadata = archive.getMetadata(metadataUIDs[0]) //select the FRET archive
 		def background_aemaex = []
		def background_aemdex = []
		def background_demdex = []
		archive.molecules()
		.filter{molecule -> molecule.getMetadataUID().equals(metadataUIDs[0])}
		.filter{molecule -> molecule.hasTag("Active_single")}	//only use molecules tagged 'Active_single'
		.forEach{molecule ->		//loop through all molecules
			if (molecule.hasSegmentsTable("T","0","active")) {
	  			MarsTable table = molecule.getSegmentsTable("T","0","active")
	  			double Ibackground = table.min("A")
	  			background_aemaex.add(Ibackground)
			}
	  		if (molecule.hasSegmentsTable("T","1 Green","active")) {
	  			MarsTable table2 = molecule.getSegmentsTable("T","1 Green","active")
	  			double Ibackground2 = table2.min("A")
	  			background_demdex.add(Ibackground2)
			}
	  		if (molecule.hasSegmentsTable("T","1 Red","active")) {
	  			MarsTable table3 = molecule.getSegmentsTable("T","1 Red","active")
	  			double Ibackground3 = table3.min("A")
	  			background_aemdex.add(Ibackground3)
	  	 	}
  		}
 
  	// Calculate the background correction factors and store in the metadata
  	double mean_aemaex = background_aemaex.sum() / background_aemaex.size()	//calculate the mean background values
	double mean_aemdex = background_aemdex.sum() / background_aemdex.size()
	double mean_demdex = background_demdex.sum() / background_demdex.size()	
	metadata.setParameter("background_aemaex", mean_aemaex) 
	metadata.setParameter("background_aemdex", mean_aemdex)
	metadata.setParameter("background_demdex", mean_demdex)

	// Add the background correction factors to the AO and DO archive as well
	if (metadataUIDs[0] != null && metadataUIDs[1] != null) {
		def metadataAO = archive.getMetadata(metadataUIDs[1]) //select the AO archive
		metadataAO.setParameter("background_aemaex", mean_aemaex) 
		metadataAO.setParameter("background_aemdex", mean_aemdex)
		metadataAO.setParameter("background_demdex", mean_demdex)
	}

	if (metadataUIDs[0] != null && metadataUIDs[2] != null) {
		def metadataDO = archive.getMetadata(metadataUIDs[2]) //select the DO archive
		metadataDO.setParameter("background_aemaex", mean_aemaex) 
		metadataDO.setParameter("background_aemdex", mean_aemdex)
		metadataDO.setParameter("background_demdex", mean_demdex)
	}	
	}
}

// Correct all I values with the background correction factor stored in the metadata record

for (def mUID : archive.getMetadataUIDs()) { //iterate through metadata records
	def metadata = archive.getMetadata(mUID)
	double background_Iaemaex = metadata.getParameter("background_aemaex")
	double background_Iaemdex = metadata.getParameter("background_aemdex")
	double background_Idemdex = metadata.getParameter("background_demdex")
	println(background_Iaemaex)

	archive.molecules()
	.filter{molecule -> molecule.getMetadataUID().equals(mUID)}
	.forEach{molecule -> //iterate through molecule records
		double Iaemdex = molecule.getParameter("Iaemdex")
  		double Idemdex = molecule.getParameter("Idemdex")
 		double Iaemaex = molecule.getParameter("Iaemaex")
 		double iiIaemdex = Iaemdex - background_Iaemdex
  		double iiIdemdex = Idemdex - background_Idemdex
  		double iiIaemaex = Iaemaex - background_Iaemaex
  		molecule.setParameter("iiIaemdex",iiIaemdex)
  		molecule.setParameter("iiIdemdex",iiIdemdex)
  		molecule.setParameter("iiIaemaex",iiIaemaex)
	}
	
}








