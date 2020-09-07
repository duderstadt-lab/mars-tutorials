//Converts the units from pixels into actual length (f.e. base pairs)

#@ MoleculeArchive archive

import de.mpg.biochem.mars.molecule.*
import de.mpg.biochem.mars.table.*
import org.scijava.table.*

double converter = 500 //set here how many pixels correspond to how many base pairs
String fromColumn = "y"
String toColumn = "bps"

//This will ensure you don't over write changes in the window
//Should not cause a problem in juptyer notebook but is also not needed..
archive.lock()

archive.forEach{ molecule ->
   MarsTable table = molecule.getDataTable()
   //Check if the column exists already
   //If not add it
   //Otherwise we overwrite the values
   if (!table.hasColumn(toColumn))
      table.appendColumn(toColumn)

   for (int row=0;row<table.getRowCount();row++)
      table.setValue(toColumn, row, table.getValue(fromColumn, row)*converter)

   archive.put(molecule)
}

archive.unlock()
